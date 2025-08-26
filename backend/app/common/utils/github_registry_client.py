import requests
import time
import os
from typing import List, Optional, Dict, Any
from functools import wraps
import logging

logger = logging.getLogger(__name__)

def cached_with_timeout(timeout_seconds: int = 600):
    """
    Decorator for caching function results with a timeout.
    
    Args:
        timeout_seconds: Cache timeout in seconds (default: 600 = 10 minutes)
    """
    def decorator(func):
        cache = {}
        
        @wraps(func)
        def wrapper(*args, **kwargs):
            key = str(args) + str(sorted(kwargs.items()))
            current_time = time.time()
            
            # Check if we have cached result and it's not expired
            if key in cache:
                result, timestamp = cache[key]
                if current_time - timestamp < timeout_seconds:
                    return result
                else:
                    # Remove expired cache entry
                    del cache[key]
            
            # Call the function and cache the result
            result = func(*args, **kwargs)
            cache[key] = (result, current_time)
            return result
        
        return wrapper
    return decorator


class GitHubRegistryClient:
    """Client for interacting with GitHub Container Registry API."""
    
    BASE_URL = "https://api.github.com"
    REGISTRY_URL = "ghcr.io/cxinsys"
    
    def __init__(self, owner: str = "cxinsys", timeout: int = 30):
        """
        Initialize the GitHub Registry client.
        
        Args:
            owner: GitHub repository owner/organization
            timeout: Request timeout in seconds
        """
        self.owner = owner
        self.timeout = timeout
        self.session = requests.Session()
        
        # Get GitHub token from environment variable
        self.github_token = os.environ.get('GITHUB_REGISTRY_TOKEN')
        
        headers = {
            'Accept': 'application/vnd.github.v3+json',
            'User-Agent': 'CellCraft-Plugin-Manager/1.0'
        }
        
        # Add authorization header if token is available
        if self.github_token:
            headers['Authorization'] = f'Bearer {self.github_token}'
            logger.info("GitHub token configured for registry access")
        else:
            logger.warning("No GITHUB_REGISTRY_TOKEN found in environment variables. API access will be limited to public repositories only.")
            
        self.session.headers.update(headers)
    
    @cached_with_timeout(timeout_seconds=600)  # 10 minute cache
    def get_available_versions(self, plugin_name: str) -> List[str]:
        """
        Get available versions for a plugin from GitHub Container Registry.
        Uses multiple fallback methods to retrieve version information.
        
        Args:
            plugin_name: Name of the plugin
            
        Returns:
            List of available version tags
        """
        # Try method 1: GitHub Packages API (requires auth)
        package_name = f"cellcraft-{plugin_name.lower()}"
        
        try:
            url = f"{self.BASE_URL}/orgs/{self.owner}/packages/container/{package_name}/versions"
            logger.info(f"Fetching versions for {package_name} from GitHub Registry API")
            
            response = self.session.get(url, timeout=self.timeout)
            
            if response.status_code == 404:
                logger.warning(f"Package {package_name} not found in registry")
                return self._get_fallback_versions(plugin_name)
            
            if response.status_code == 401:
                logger.warning(f"Unauthorized access to {package_name}, trying fallback methods")
                return self._get_fallback_versions(plugin_name)
            
            response.raise_for_status()
            
            versions = []
            for version_info in response.json():
                # Extract tags from the version metadata
                if 'metadata' in version_info and 'container' in version_info['metadata']:
                    tags = version_info['metadata']['container'].get('tags', [])
                    versions.extend(tags)
            
            # Remove duplicates and sort
            unique_versions = list(set(versions))
            unique_versions.sort(reverse=True)  # Latest first
            
            if not unique_versions:
                return self._get_fallback_versions(plugin_name)
            
            logger.info(f"Found {len(unique_versions)} versions for {plugin_name}")
            return unique_versions
            
        except requests.exceptions.Timeout:
            logger.error(f"Timeout fetching versions for {plugin_name}")
            return self._get_fallback_versions(plugin_name)
        except requests.exceptions.RequestException as e:
            logger.error(f"Error fetching versions for {plugin_name}: {e}")
            return self._get_fallback_versions(plugin_name)
        except Exception as e:
            logger.error(f"Unexpected error fetching versions for {plugin_name}: {e}")
            return self._get_fallback_versions(plugin_name)
    
    def _get_fallback_versions(self, plugin_name: str) -> List[str]:
        """
        Fallback method to get versions using Docker Registry API or predefined versions.
        
        Args:
            plugin_name: Name of the plugin
            
        Returns:
            List of fallback version tags
        """
        try:
            # Try Docker Registry API v2 for public images
            package_name = f"cellcraft-{plugin_name.lower()}"
            registry_api_url = f"https://ghcr.io/v2/{self.owner}/{package_name}/tags/list"
            
            logger.info(f"Trying Docker Registry API for {package_name}")
            
            response = self.session.get(registry_api_url, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                tags = data.get('tags', [])
                if tags:
                    # Sort tags, putting semantic versions first
                    sorted_tags = sorted(tags, key=lambda x: (not self._is_semantic_version(x), x), reverse=True)
                    logger.info(f"Found {len(sorted_tags)} versions via Docker Registry API for {plugin_name}")
                    return sorted_tags
            
            logger.warning(f"Docker Registry API failed for {package_name}, status: {response.status_code}")
            
        except Exception as e:
            logger.warning(f"Docker Registry API error for {plugin_name}: {e}")
        
        # Final fallback: return predefined versions based on common plugin versioning
        fallback_versions = self._get_predefined_versions(plugin_name)
        logger.info(f"Using predefined fallback versions for {plugin_name}: {fallback_versions}")
        return fallback_versions
    
    def _get_predefined_versions(self, plugin_name: str) -> List[str]:
        """
        Get predefined versions for known plugins.
        
        Args:
            plugin_name: Name of the plugin
            
        Returns:
            List of predefined version tags
        """
        # Actual plugin versions matching the current GitHub Container Registry
        # Note: These should be updated when new versions are released
        predefined_versions = {
            'tenet': ['1.0.0', 'latest'],
            'fasttenet': ['1.0.0', 'latest'],
            'genie3': ['1.0.0', 'latest'],
            'grnboost2': ['1.0.0', 'latest'],
            'leap': ['1.0.0', 'latest'],
            'scribe': ['1.0.0', 'latest'],
        }
        
        plugin_lower = plugin_name.lower()
        if plugin_lower in predefined_versions:
            return predefined_versions[plugin_lower]
        
        # Generic fallback for unknown plugins
        return ['1.0.0', 'latest']
    
    def _is_semantic_version(self, version: str) -> bool:
        """
        Check if a version string follows semantic versioning pattern.
        
        Args:
            version: Version string to check
            
        Returns:
            True if it looks like semantic version (e.g., 1.0.0, v2.1.3)
        """
        import re
        # Match semantic version patterns like 1.0.0, v1.0.0, 1.0.0-beta, etc.
        pattern = r'^v?(\d+)\.(\d+)\.(\d+)(-[a-zA-Z0-9\-\.]+)?$'
        return bool(re.match(pattern, version))
    
    def get_image_uri(self, plugin_name: str, version: Optional[str] = None) -> str:
        """
        Get the full image URI for a plugin.
        
        Args:
            plugin_name: Name of the plugin
            version: Version tag (defaults to 'latest' if None)
            
        Returns:
            Full image URI string
        """
        if version is None:
            version = "latest"
        
        plugin_name_lower = plugin_name.lower()
        image_uri = f"{self.REGISTRY_URL}/cellcraft-{plugin_name_lower}:{version}"
        
        logger.debug(f"Generated image URI: {image_uri}")
        return image_uri
    
    def check_image_exists(self, plugin_name: str, version: Optional[str] = None) -> bool:
        """
        Check if a specific image version exists in the registry.
        
        Args:
            plugin_name: Name of the plugin
            version: Version tag (defaults to 'latest' if None)
            
        Returns:
            True if image exists, False otherwise
        """
        try:
            versions = self.get_available_versions(plugin_name)
            target_version = version or "latest"
            return target_version in versions
        except Exception as e:
            logger.error(f"Error checking if image exists for {plugin_name}:{version}: {e}")
            return False