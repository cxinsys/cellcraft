"""
Plugin Synchronization Manager for CellCraft
Manages synchronization between cellcraft-plugin repository, database, and Docker images
"""

import os
import json
import re
import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import pandas as pd
from sqlalchemy.orm import Session
from datetime import datetime

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
        Get current branch from version.json file
        
        Returns:
            str: Current branch name
        """
        try:
            version_file = self.official_plugin_path / "version.json"
            
            if version_file.exists():
                with open(version_file, 'r') as f:
                    version_info = json.load(f)
                    branch = version_info.get("branch", "unknown")
                    logger.info(f"Current plugin branch from version.json: {branch}")
                    return branch
            else:
                # Default branch if version.json doesn't exist
                logger.warning("version.json not found, using default branch")
                print("version.json not found, using default branch")
                return "release/plugins-v1.0"
                
        except Exception as e:
            logger.error(f"Failed to get current branch: {e}")
            # Return a default value instead of raising
            return "release/plugins-v1.0"
    
    def get_available_branches(self) -> List[str]:
        """
        Get list of available release branches (fixed list for now)
        
        Returns:
            List[str]: List of branch names matching release/plugins-v* pattern
        """
        # Since we can't query git in runtime, return a fixed list
        # This could be updated via configuration file or environment variables
        branches = [
            "release/plugins-v1.0",
            "release/plugins-v1.1",
            "release/plugins-v2.0"
        ]
        logger.info(f"Available release branches (configured): {branches}")
        return branches
    
    def switch_branch(self, branch: str) -> bool:
        """
        Branch switching is not supported in runtime.
        This should be done during Docker image build.
        
        Args:
            branch: Target branch name (e.g., "release/plugins-v1.0")
            
        Returns:
            bool: Always returns False as runtime switching is not supported
        """
        logger.warning(f"Branch switching to {branch} is not supported in runtime. " +
                      "Please rebuild the Docker image with the desired branch.")
        return False
    
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
            
            # Get version info from file
            version_info = {}
            version_file = self.official_plugin_path / "version.json"
            if version_file.exists():
                with open(version_file, 'r') as f:
                    version_info = json.load(f)
            
            current_branch = version_info.get("branch", "release/plugins-v1.0")
            version = self.extract_version_from_branch(current_branch)
            
            # Update version for all official plugins
            self.update_plugin_versions(version)
            
            return {
                "success": True,
                "branch": current_branch,
                "version": version,
                "version_info": version_info,
                "message": f"Successfully synced plugins with version {version}"
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
            # Get version info from file
            version_info = {}
            version_file = self.official_plugin_path / "version.json"
            if version_file.exists():
                with open(version_file, 'r') as f:
                    version_info = json.load(f)
            
            current_branch = version_info.get("branch", "unknown")
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
                "version_info": version_info,
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