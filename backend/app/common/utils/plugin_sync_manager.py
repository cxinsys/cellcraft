"""
Plugin Synchronization Manager for CellCraft
Manages synchronization between cellcraft-plugin repository, database, and Docker images
"""

import os
import subprocess
import re
import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import pandas as pd
from sqlalchemy.orm import Session

from app.database import models
from app.database.conn import SessionLocal, initialize_plugins_from_csv
from app.database.crud import crud_plugin

logger = logging.getLogger(__name__)


class PluginSyncManager:
    """Manages plugin repository synchronization and version control"""
    
    def __init__(self):
        self.official_plugin_path = Path("./plugin/official")
        self.plugins_csv_path = self.official_plugin_path / "plugins.csv"
        
    def get_current_branch(self) -> str:
        """
        Get current branch of cellcraft-plugin submodule
        
        Returns:
            str: Current branch name
        """
        try:
            # Change to submodule directory
            original_dir = os.getcwd()
            os.chdir(self.official_plugin_path)
            
            # Get current branch
            result = subprocess.run(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                capture_output=True,
                text=True,
                check=True
            )
            
            branch = result.stdout.strip()
            logger.info(f"Current plugin branch: {branch}")
            return branch
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to get current branch: {e}")
            raise
        finally:
            os.chdir(original_dir)
    
    def get_available_branches(self) -> List[str]:
        """
        Get list of available release branches from remote
        
        Returns:
            List[str]: List of branch names matching release/plugins-v* pattern
        """
        try:
            original_dir = os.getcwd()
            os.chdir(self.official_plugin_path)
            
            # Fetch latest from remote
            subprocess.run(["git", "fetch", "origin"], check=True, capture_output=True)
            
            # Get remote branches
            result = subprocess.run(
                ["git", "branch", "-r"],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Filter for release branches
            branches = []
            for line in result.stdout.splitlines():
                line = line.strip()
                if "origin/release/plugins-v" in line:
                    # Extract branch name without 'origin/'
                    branch = line.replace("origin/", "").strip()
                    branches.append(branch)
            
            logger.info(f"Available release branches: {branches}")
            return sorted(branches)
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to get available branches: {e}")
            raise
        finally:
            os.chdir(original_dir)
    
    def switch_branch(self, branch: str) -> bool:
        """
        Switch cellcraft-plugin submodule to specified branch
        
        Args:
            branch: Target branch name (e.g., "release/plugins-v1.0")
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            original_dir = os.getcwd()
            os.chdir(self.official_plugin_path)
            
            # Fetch latest changes
            subprocess.run(["git", "fetch", "origin"], check=True, capture_output=True)
            
            # Checkout branch
            subprocess.run(
                ["git", "checkout", branch],
                check=True,
                capture_output=True
            )
            
            # Pull latest changes
            subprocess.run(
                ["git", "pull", "origin", branch],
                check=True,
                capture_output=True
            )
            
            logger.info(f"Successfully switched to branch: {branch}")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to switch branch: {e}")
            return False
        finally:
            os.chdir(original_dir)
    
    def extract_version_from_branch(self, branch: str) -> str:
        """
        Extract version number from branch name
        
        Args:
            branch: Branch name (e.g., "release/plugins-v1.0")
            
        Returns:
            str: Version string (e.g., "1.0")
        """
        # Match pattern release/plugins-v{version}
        match = re.match(r"release/plugins-v(\d+\.\d+)", branch)
        if match:
            return match.group(1)
        else:
            logger.warning(f"Could not extract version from branch: {branch}")
            return "latest"
    
    def sync_plugins_to_database(self) -> Dict[str, any]:
        """
        Synchronize plugins.csv to database
        
        Returns:
            Dict containing sync results
        """
        try:
            # Check if plugins.csv exists
            if not self.plugins_csv_path.exists():
                raise FileNotFoundError(f"plugins.csv not found at {self.plugins_csv_path}")
            
            # Use existing function to initialize plugins
            initialize_plugins_from_csv(str(self.plugins_csv_path))
            
            # Get current branch and version
            current_branch = self.get_current_branch()
            version = self.extract_version_from_branch(current_branch)
            
            # Update version for all official plugins
            self.update_plugin_versions(version)
            
            return {
                "success": True,
                "branch": current_branch,
                "version": version,
                "message": f"Successfully synced plugins from branch {current_branch}"
            }
            
        except Exception as e:
            logger.error(f"Failed to sync plugins: {e}")
            return {
                "success": False,
                "error": str(e)
            }
    
    def update_plugin_versions(self, version: str) -> None:
        """
        Update version field for all official plugins in database
        
        Args:
            version: Version string to set for all official plugins
        """
        db = SessionLocal()
        try:
            # Update all official plugins
            official_plugins = db.query(models.Plugin).filter_by(source="official").all()
            
            for plugin in official_plugins:
                plugin.version = version
                logger.info(f"Updated {plugin.name} to version {version}")
            
            db.commit()
            logger.info(f"Updated {len(official_plugins)} official plugins to version {version}")
            
        except Exception as e:
            db.rollback()
            logger.error(f"Failed to update plugin versions: {e}")
            raise
        finally:
            db.close()
    
    def get_sync_status(self) -> Dict[str, any]:
        """
        Get current synchronization status
        
        Returns:
            Dict containing current status information
        """
        try:
            current_branch = self.get_current_branch()
            version = self.extract_version_from_branch(current_branch)
            
            # Check database plugin versions
            db = SessionLocal()
            try:
                official_plugins = db.query(models.Plugin).filter_by(source="official").all()
                db_versions = {p.name: p.version for p in official_plugins}
            finally:
                db.close()
            
            return {
                "repository_branch": current_branch,
                "repository_version": version,
                "database_plugin_count": len(db_versions),
                "database_versions": db_versions,
                "plugins_csv_exists": self.plugins_csv_path.exists(),
                "plugins_csv_modified": self.plugins_csv_path.stat().st_mtime if self.plugins_csv_path.exists() else None
            }
            
        except Exception as e:
            logger.error(f"Failed to get sync status: {e}")
            return {
                "error": str(e)
            }