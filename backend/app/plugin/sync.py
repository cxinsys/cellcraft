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

from app import models
from app.db.session import SessionLocal, get_db_session
from app.plugin import crud as crud_plugin
from app.core.security import get_password_hash

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

        Note:
            Plugin versions are now managed per-plugin via plugins.csv.
            Each plugin has its own version (e.g., GENIE3: 1.0, FastTENET: 1.1).
            The version.json "version" field represents the plugin bundle version,
            NOT individual plugin versions.
        """
        try:
            # Check if plugins.csv exists
            if not self.plugins_csv_path.exists():
                raise FileNotFoundError(f"plugins.csv not found at {self.plugins_csv_path}")

            # Use existing function to initialize plugins
            # This properly reads individual plugin versions from plugins.csv
            initialize_plugins_from_csv(str(self.plugins_csv_path))

            # Get version info from file (for bundle version tracking only)
            version_info = {}
            version_file = self.official_plugin_path / "version.json"
            if version_file.exists():
                with open(version_file, 'r') as f:
                    version_info = json.load(f)

            current_branch = version_info.get("branch", "release/plugins-v1.0")
            bundle_version = self.extract_version_from_branch(current_branch)

            # NOTE: Do NOT call update_plugin_versions() here!
            # Individual plugin versions are already set correctly from plugins.csv
            # update_plugin_versions() would incorrectly overwrite all plugins with the bundle version

            return {
                "success": True,
                "branch": current_branch,
                "bundle_version": bundle_version,
                "version_info": version_info,
                "message": f"Successfully synced plugins from CSV (bundle: {bundle_version})"
            }
            
        except Exception as e:
            logger.error(f"Failed to sync plugins: {e}")
            return {
                "success": False,
                "error": str(e)
            }
    
    def update_plugin_versions(self, version: str) -> None:
        """
        DEPRECATED: Do not use this method.

        This method incorrectly overwrites individual plugin versions with a single
        bundle version. Each plugin has its own version defined in plugins.csv
        (e.g., GENIE3: 1.0, FastTENET: 1.1), and this method would overwrite
        all of them with the bundle version from version.json.

        Plugin versions should ONLY be updated via:
        1. Editing plugins.csv with correct individual versions
        2. Running initialize_plugins_from_csv() which respects per-plugin versions

        Args:
            version: Version string (IGNORED - method is deprecated)

        Warning:
            This method is kept for backward compatibility but logs a warning
            and does nothing. It will be removed in a future version.
        """
        logger.warning(
            "DEPRECATED: update_plugin_versions() is deprecated and does nothing. "
            "Plugin versions are now managed per-plugin via plugins.csv. "
            "Do not call this method - it would incorrectly overwrite individual versions."
        )
        # DO NOT update versions - each plugin has its own version from plugins.csv
        # The following code is intentionally disabled:
        # db = SessionLocal()
        # try:
        #     official_plugins = db.query(models.Plugin).filter_by(source="official").all()
        #     for plugin in official_plugins:
        #         plugin.version = version  # BAD: overwrites individual versions!
        #     db.commit()
        # finally:
        #     db.close()
    
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

# Moved from app.db.session in PR-4 (phase-2-domain-skeleton).
def initialize_plugins_from_csv(csv_file_path: str = None):
    """
    Initialize plugins from CSV file. 
    If no path is provided, reads from official plugins CSV.
    
    Args:
        csv_file_path (str, optional): Path to CSV file. If None, uses official plugins CSV.
        
    Returns:
        int: Number of plugins successfully initialized
    """
    # Default to official plugins CSV if no path provided
    if csv_file_path is None:
        csv_file_path = "./plugin/official/plugins.csv"
    
    # Check if CSV file exists
    if not os.path.exists(csv_file_path):
        print(f"Warning: Plugin CSV file not found at {csv_file_path}")
        return 0
    
    # CSV 파일 읽기
    try:
        df = pd.read_csv(csv_file_path)
    except Exception as e:
        print(f"Error reading CSV file {csv_file_path}: {e}")
        return 0

    # 세션 시작 - get_db_session() context manager 사용
    try:
        with get_db_session() as session:
            # 관리자 사용자 추가
            existing_user = session.query(models.User).filter_by(username="admin").first()
            if not existing_user:
                hashed_password = get_password_hash("cellcraft2024!")
                user = models.User(
                    username="admin",
                    email="cellcraft@cellcraft.com",
                    hashed_password=hashed_password,
                    is_active=True,
                    is_superuser=True
                )
                session.add(user)
                print("Created admin user")

            # 데이터 추가
            plugins_added = 0
            for index, row in df.iterrows():
                # 플러그인이 이미 존재하는지 확인 (name과 source로 확인)
                source = "official" if "official" in csv_file_path else "local"
                existing_plugin = session.query(models.Plugin).filter_by(
                    name=row['name'],
                    source=source
                ).first()

                # JSON 필드 파싱 (공통 로직)
                dependencies = {}
                if pd.notna(row['dependencies']) and str(row['dependencies']).strip():
                    try:
                        dependencies = json.loads(str(row['dependencies']))
                    except json.JSONDecodeError:
                        print(f"Warning: Invalid JSON in dependencies for plugin {row['name']}, using empty dict")
                        dependencies = {}

                drawflow = {}
                if pd.notna(row['drawflow']) and str(row['drawflow']).strip():
                    try:
                        drawflow = json.loads(str(row['drawflow']))
                    except json.JSONDecodeError:
                        print(f"Warning: Invalid JSON in drawflow for plugin {row['name']}, using empty dict")
                        drawflow = {}

                rules = {}
                if pd.notna(row['rules']) and str(row['rules']).strip():
                    try:
                        rules = json.loads(str(row['rules']))
                    except json.JSONDecodeError:
                        print(f"Warning: Invalid JSON in rules for plugin {row['name']}, using empty dict")
                        rules = {}

                # 추가 필드 파싱 (공통 로직)
                plugin_type = None
                if 'plugin_type' in row and pd.notna(row['plugin_type']) and str(row['plugin_type']).strip():
                    plugin_type = str(row['plugin_type']).strip()

                use_gpu = False
                if 'use_gpu' in row and pd.notna(row['use_gpu']):
                    if isinstance(row['use_gpu'], str):
                        use_gpu = row['use_gpu'].lower().strip() in ('true', '1', 'yes')
                    else:
                        use_gpu = bool(row['use_gpu'])

                if not existing_plugin:
                    try:
                        # Determine plugin path and attributes based on source
                        if source == "official":
                            plugin_path = f"./plugin/official/{row['name']}"
                            is_editable = False
                            # Extract version and submodule path if available
                            version = row.get('version', None) if 'version' in row else None
                            submodule_path = row.get('submodule_path', None) if 'submodule_path' in row else None
                        else:
                            plugin_path = f"./plugin/local/{row['name']}"
                            is_editable = True
                            version = None
                            submodule_path = None

                        # Handle plugin_path from CSV if provided
                        if 'plugin_path' in row and pd.notna(row['plugin_path']) and str(row['plugin_path']).strip():
                            plugin_path = str(row['plugin_path']).strip()

                        # Handle is_editable from CSV if provided (for consistency)
                        if 'is_editable' in row and pd.notna(row['is_editable']):
                            if isinstance(row['is_editable'], str):
                                is_editable = row['is_editable'].lower().strip() in ('true', '1', 'yes')
                            else:
                                is_editable = bool(row['is_editable'])

                        # Handle author field properly
                        author = str(row['author']) if pd.notna(row['author']) else 'Unknown'

                        plugin = models.Plugin(
                            name=str(row['name']),
                            description=str(row['description']) if pd.notna(row['description']) else '',
                            author=author,
                            plugin_path=plugin_path,
                            plugin_type=plugin_type,
                            dependencies=dependencies,
                            drawflow=drawflow,
                            rules=rules,
                            use_gpu=use_gpu,
                            source=source,
                            is_editable=is_editable,
                            version=version,
                            submodule_path=submodule_path
                        )
                        session.add(plugin)
                        plugins_added += 1
                        print(f"Added {source} plugin: {row['name']} (type: {plugin_type}, GPU: {use_gpu})")
                    except Exception as e:
                        print(f"Error creating plugin {row['name']}: {e}")
                        continue
                else:
                    # Update existing plugin with CSV data
                    try:
                        updated = False

                        # Update fields that might have changed
                        if existing_plugin.description != str(row['description']):
                            existing_plugin.description = str(row['description']) if pd.notna(row['description']) else ''
                            updated = True

                        if existing_plugin.plugin_type != plugin_type:
                            existing_plugin.plugin_type = plugin_type
                            updated = True

                        if existing_plugin.use_gpu != use_gpu:
                            existing_plugin.use_gpu = use_gpu
                            updated = True

                        # Update JSON fields if they differ
                        if existing_plugin.dependencies != dependencies:
                            existing_plugin.dependencies = dependencies
                            updated = True

                        if existing_plugin.drawflow != drawflow:
                            existing_plugin.drawflow = drawflow
                            updated = True

                        if existing_plugin.rules != rules:
                            existing_plugin.rules = rules
                            updated = True

                        # Update version info for official plugins
                        if source == "official":
                            version = row.get('version', None) if 'version' in row else None
                            submodule_path = row.get('submodule_path', None) if 'submodule_path' in row else None

                            if existing_plugin.version != version:
                                existing_plugin.version = version
                                updated = True

                            if existing_plugin.submodule_path != submodule_path:
                                existing_plugin.submodule_path = submodule_path
                                updated = True

                        if updated:
                            print(f"Updated {source} plugin: {row['name']} (type: {plugin_type}, GPU: {use_gpu})")
                        else:
                            print(f"Plugin {row['name']} from {source} is up to date")

                    except Exception as e:
                        print(f"Error updating plugin {row['name']}: {e}")
                        continue

            # 커밋은 context manager가 자동 처리
            print(f"Successfully initialized {plugins_added} plugins from {csv_file_path}")
            return plugins_added
    except Exception as e:
        # rollback/close는 context manager가 자동 처리
        print(f"Error during plugin initialization: {e}")
        return 0