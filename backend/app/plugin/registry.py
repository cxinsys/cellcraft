import requests
import time
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
        
        # Configure headers for public repository access
        headers = {
            'Accept': 'application/vnd.github.v3+json',
            'User-Agent': 'CellCraft-Plugin-Manager/1.0'
        }
        
        self.session.headers.update(headers)
        logger.info("GitHub registry client configured for public repository access")
    
    @cached_with_timeout(timeout_seconds=600)  # 10 minute cache
    def get_available_versions(self, plugin_name: str) -> List[str]:
        """
        Get available versions for a plugin from GitHub Container Registry.
        Uses Docker Registry v2 API directly to avoid authentication issues.
        
        Args:
            plugin_name: Name of the plugin
            
        Returns:
            List of available version tags
        """
        # Primary method: Docker Registry v2 API (public access)
        package_name = f"cellcraft-{plugin_name.lower()}"
        
        try:
            # Try Docker Registry API v2 for public images first
            registry_api_url = f"https://ghcr.io/v2/{self.owner}/{package_name}/tags/list"
            
            logger.info(f"Fetching versions for {package_name} from Docker Registry v2 API")
            
            response = self.session.get(registry_api_url, timeout=self.timeout)
            
            if response.status_code == 200:
                data = response.json()
                tags = data.get('tags', [])
                if tags:
                    # Sort tags, putting semantic versions first
                    sorted_tags = sorted(tags, key=lambda x: (not self._is_semantic_version(x), x), reverse=True)
                    logger.info(f"Found {len(sorted_tags)} versions via Docker Registry API for {plugin_name}")
                    return sorted_tags
                else:
                    logger.warning(f"No tags found for {package_name}")
            elif response.status_code == 404:
                logger.warning(f"Package {package_name} not found in registry")
            else:
                logger.warning(f"Docker Registry API failed for {package_name}, status: {response.status_code}")
            
        except requests.exceptions.Timeout:
            logger.error(f"Timeout fetching versions for {plugin_name}")
        except requests.exceptions.RequestException as e:
            logger.error(f"Docker Registry API error for {plugin_name}: {e}")
        except Exception as e:
            logger.error(f"Unexpected error fetching versions for {plugin_name}: {e}")
        
        # Fallback to predefined versions
        return self._get_fallback_versions(plugin_name)
    
    def _get_fallback_versions(self, plugin_name: str) -> List[str]:
        """
        Fallback method to return predefined versions for known plugins.
        
        Args:
            plugin_name: Name of the plugin
            
        Returns:
            List of fallback version tags
        """
        # Use predefined versions for CPU-only branch plugins
        fallback_versions = self._get_predefined_versions(plugin_name)
        logger.info(f"Using predefined fallback versions for {plugin_name}: {fallback_versions}")
        return fallback_versions
    
    def _get_predefined_versions(self, plugin_name: str) -> List[str]:
        """
        Get predefined versions for known plugins in CPU-only mode.
        
        Args:
            plugin_name: Name of the plugin
            
        Returns:
            List of predefined version tags
        """
        # CPU-only compatible plugin versions for release/plugins-v1.0-cpu branch
        # Note: FastSCODE and FastTENET are GPU-only and excluded in CPU mode
        cpu_compatible_versions = {
            'tenet': ['1.0', 'latest'],
            'genie3': ['1.0', 'latest'],
            'grnboost2': ['1.0', 'latest'],
            'grnviz': ['1.0', 'latest'],
            'leap': ['1.0', 'latest'],
            'scribe': ['1.0', 'latest'],
        }
        
        plugin_lower = plugin_name.lower()
        if plugin_lower in cpu_compatible_versions:
            return cpu_compatible_versions[plugin_lower]
        
        # For GPU-only plugins, return empty list to indicate unavailability
        gpu_only_plugins = {'fastscode', 'fasttenet'}
        if plugin_lower in gpu_only_plugins:
            logger.warning(f"Plugin {plugin_name} is GPU-only and not available in CPU mode")
            return []
        
        # Generic fallback for unknown plugins
        return ['1.0', 'latest']
    
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