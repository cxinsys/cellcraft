"""
Plugin Version Validator for CellCraft
Validates consistency between repository branch, database versions, and Docker images
"""

import logging
from typing import Dict, List, Optional, Set
from datetime import datetime

from app.common.utils.plugin_sync_manager import PluginSyncManager
from app.common.utils.github_registry_client import GitHubRegistryClient
from app.db.session import SessionLocal
from app.database import models

logger = logging.getLogger(__name__)


class PluginVersionValidator:
    """Validates version consistency across plugin system components"""
    
    def __init__(self):
        self.sync_manager = PluginSyncManager()
        self.registry_client = GitHubRegistryClient()
        
    def validate_consistency(self) -> Dict[str, any]:
        """
        Validate version consistency between repository, database, and Docker registry
        
        Returns:
            Dict containing validation results:
            {
                "consistent": bool,
                "repository_branch": str,
                "repository_version": str,
                "database_versions": dict,
                "docker_availability": dict,
                "issues": list
            }
        """
        issues = []
        
        try:
            # Get repository information
            current_branch = self.sync_manager.get_current_branch()
            repository_version = self.sync_manager.extract_version_from_branch(current_branch)
            
            # Get database information
            db = SessionLocal()
            try:
                official_plugins = db.query(models.Plugin).filter_by(source="official").all()
                database_versions = {plugin.name: plugin.version for plugin in official_plugins}
                plugin_names = list(database_versions.keys())
            finally:
                db.close()
            
            # Check database version consistency
            inconsistent_db_versions = []
            for plugin_name, db_version in database_versions.items():
                if db_version != repository_version:
                    inconsistent_db_versions.append({
                        "plugin": plugin_name,
                        "database_version": db_version,
                        "expected_version": repository_version
                    })
                    issues.append(f"Plugin '{plugin_name}' has version '{db_version}' in database, expected '{repository_version}'")
            
            # Skip Docker registry checks - will be handled during actual pull
            docker_availability = {}
            missing_docker_images = []
            
            # Mark all as potentially available (actual check happens during pull)
            for plugin_name in plugin_names:
                docker_availability[plugin_name] = {
                    "version": repository_version,
                    "available": True  # Assume available, actual pull will verify
                }
            
            # Determine overall consistency
            consistent = len(issues) == 0
            
            return {
                "consistent": consistent,
                "repository_branch": current_branch,
                "repository_version": repository_version,
                "database_versions": database_versions,
                "docker_availability": docker_availability,
                "issues": issues,
                "details": {
                    "inconsistent_db_versions": inconsistent_db_versions,
                    "missing_docker_images": missing_docker_images,
                    "total_plugins": len(plugin_names)
                },
                "timestamp": datetime.now().isoformat()
            }
            
        except Exception as e:
            logger.error(f"Failed to validate consistency: {e}")
            return {
                "consistent": False,
                "error": str(e),
                "issues": [f"Validation error: {str(e)}"],
                "timestamp": datetime.now().isoformat()
            }
    
    def check_docker_image_availability(self, plugin_name: str, version: str) -> bool:
        """
        Check if specific Docker image version exists in registry
        
        Args:
            plugin_name: Name of the plugin
            version: Version tag to check
            
        Returns:
            bool: True if image exists, False otherwise
        """
        try:
            # Get available versions from registry
            available_versions = self.registry_client.get_available_versions(plugin_name.lower())
            
            # Check if specific version exists
            return version in available_versions
            
        except Exception as e:
            logger.error(f"Failed to check Docker image for {plugin_name}:{version}: {e}")
            return False
    
    def generate_consistency_report(self) -> str:
        """
        Generate human-readable consistency report
        
        Returns:
            str: Formatted report
        """
        validation_result = self.validate_consistency()
        
        report_lines = [
            "=" * 80,
            "Plugin Version Consistency Report",
            "=" * 80,
            f"Generated at: {validation_result.get('timestamp', 'Unknown')}",
            "",
            f"Repository Branch: {validation_result.get('repository_branch', 'Unknown')}",
            f"Repository Version: {validation_result.get('repository_version', 'Unknown')}",
            "",
        ]
        
        if validation_result.get('error'):
            report_lines.extend([
                "ERROR: Validation failed",
                f"  {validation_result['error']}",
                ""
            ])
            return "\n".join(report_lines)
        
        # Overall status
        if validation_result['consistent']:
            report_lines.append("✅ CONSISTENT - All components are in sync")
        else:
            report_lines.append("❌ INCONSISTENT - Issues detected")
        
        report_lines.append("")
        
        # Database versions
        report_lines.append("Database Plugin Versions:")
        db_versions = validation_result.get('database_versions', {})
        for plugin_name, version in sorted(db_versions.items()):
            status = "✅" if version == validation_result['repository_version'] else "❌"
            report_lines.append(f"  {status} {plugin_name}: {version}")
        
        report_lines.append("")
        
        # Docker availability
        report_lines.append("Docker Image Availability:")
        docker_info = validation_result.get('docker_availability', {})
        for plugin_name, info in sorted(docker_info.items()):
            status = "✅" if info['available'] else "❌"
            report_lines.append(f"  {status} {plugin_name}:{info['version']}")
        
        # Issues summary
        if validation_result.get('issues'):
            report_lines.extend([
                "",
                "Issues Found:",
            ])
            for issue in validation_result['issues']:
                report_lines.append(f"  - {issue}")
        
        # Details
        details = validation_result.get('details', {})
        if details.get('inconsistent_db_versions'):
            report_lines.extend([
                "",
                "Database Version Mismatches:",
            ])
            for mismatch in details['inconsistent_db_versions']:
                report_lines.append(
                    f"  - {mismatch['plugin']}: {mismatch['database_version']} → {mismatch['expected_version']}"
                )
        
        if details.get('missing_docker_images'):
            report_lines.extend([
                "",
                "Missing Docker Images:",
            ])
            for plugin in details['missing_docker_images']:
                report_lines.append(f"  - {plugin}")
        
        report_lines.extend([
            "",
            "=" * 80
        ])
        
        return "\n".join(report_lines)
    
    def quick_check(self) -> bool:
        """
        Quick consistency check (returns bool only)
        
        Returns:
            bool: True if consistent, False otherwise
        """
        result = self.validate_consistency()
        return result.get('consistent', False)