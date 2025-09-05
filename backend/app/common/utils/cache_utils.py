"""
Visualization Caching Utilities

Provides centralized caching system for visualization results with 7-day expiry policy.
Uses symbolic links to maintain compatibility with existing file structure while 
optimizing storage through cache deduplication.
"""

import os
import json
import hashlib
import shutil
import random
import logging
from datetime import datetime, timezone, timedelta
from typing import Dict, List, Tuple, Optional, Any

logger = logging.getLogger(__name__)


def generate_cache_key(plugin_name: str, script_name: str, params: Dict[str, Any], input_files: List[str]) -> str:
    """
    Generate a unique cache key based on plugin, script, parameters, and input files.
    
    Args:
        plugin_name: Name of the visualization plugin
        script_name: Name of the visualization script
        params: Dictionary of visualization parameters
        input_files: List of input file paths
        
    Returns:
        16-character hexadecimal cache key
    """
    try:
        # Normalize parameters for consistent hashing
        param_str = json.dumps(params, sort_keys=True, separators=(',', ':'))
        
        # Use only basenames of files for consistency
        files_str = ''.join(sorted([os.path.basename(f) for f in input_files]))
        
        # Combine all components
        combined = f"{plugin_name}_{script_name}_{param_str}_{files_str}"
        
        # Generate MD5 hash and take first 16 characters
        cache_key = hashlib.md5(combined.encode('utf-8')).hexdigest()[:16]
        
        logger.debug(f"Generated cache key: {cache_key} for plugin={plugin_name}, script={script_name}")
        return cache_key
        
    except Exception as e:
        logger.error(f"Failed to generate cache key: {e}")
        # Fallback to timestamp-based key if hashing fails
        return f"fallback_{int(datetime.now().timestamp())}"


def load_cache_metadata(metadata_file: str) -> Dict[str, Any]:
    """
    Load cache metadata from JSON file. Creates default structure if file doesn't exist.
    
    Args:
        metadata_file: Path to metadata JSON file
        
    Returns:
        Dictionary containing cache metadata
    """
    default_metadata = {
        "cache_entries": {},
        "settings": {
            "max_age_days": 7,
            "last_cleanup": datetime.now(timezone.utc).isoformat()
        }
    }
    
    try:
        if os.path.exists(metadata_file):
            with open(metadata_file, 'r', encoding='utf-8') as f:
                metadata = json.load(f)
                # Ensure settings exist
                if 'settings' not in metadata:
                    metadata['settings'] = default_metadata['settings']
                return metadata
        else:
            # Create metadata file with default structure
            os.makedirs(os.path.dirname(metadata_file), exist_ok=True)
            save_cache_metadata(metadata_file, default_metadata)
            return default_metadata
            
    except (json.JSONDecodeError, IOError) as e:
        logger.warning(f"Failed to load cache metadata from {metadata_file}: {e}")
        logger.info("Creating new metadata file with default structure")
        save_cache_metadata(metadata_file, default_metadata)
        return default_metadata


def save_cache_metadata(metadata_file: str, metadata: Dict[str, Any]) -> bool:
    """
    Save cache metadata to JSON file.
    
    Args:
        metadata_file: Path to metadata JSON file
        metadata: Metadata dictionary to save
        
    Returns:
        True if successful, False otherwise
    """
    try:
        os.makedirs(os.path.dirname(metadata_file), exist_ok=True)
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        return True
        
    except IOError as e:
        logger.error(f"Failed to save cache metadata to {metadata_file}: {e}")
        return False


def get_metadata_file_path(user_path: str) -> str:
    """Get the path to the cache metadata file for a user."""
    return os.path.join(user_path, '.cache', 'metadata', 'viz_cache.json')


def get_cache_dir_path(user_path: str) -> str:
    """Get the path to the cache directory for a user."""
    return os.path.join(user_path, '.cache', 'viz')


def check_cache_with_expiry(user_path: str, cache_key: str, max_age_days: int = 7) -> Tuple[bool, Optional[str]]:
    """
    Check if cached result exists and is not expired.
    
    Args:
        user_path: User's base directory path
        cache_key: Cache key to check
        max_age_days: Maximum age in days before cache expires
        
    Returns:
        Tuple of (cache_hit: bool, cache_file_path: Optional[str])
    """
    metadata_file = get_metadata_file_path(user_path)
    cache_dir = get_cache_dir_path(user_path)
    
    try:
        metadata = load_cache_metadata(metadata_file)
        
        # Check if cache entry exists
        if cache_key not in metadata['cache_entries']:
            logger.debug(f"Cache miss: key {cache_key} not found")
            return False, None
            
        entry = metadata['cache_entries'][cache_key]
        created_at_str = entry.get('created_at')
        
        if not created_at_str:
            logger.warning(f"Cache entry {cache_key} missing created_at timestamp")
            return False, None
            
        # Parse creation time
        try:
            created_at = datetime.fromisoformat(created_at_str.replace('Z', '+00:00'))
        except ValueError as e:
            logger.warning(f"Invalid timestamp format for cache key {cache_key}: {e}")
            return False, None
            
        # Check if expired
        current_time = datetime.now(timezone.utc)
        age = current_time - created_at
        
        if age > timedelta(days=max_age_days):
            logger.info(f"Cache expired: key {cache_key}, age {age.days} days")
            # Remove expired cache entry
            remove_expired_cache_entry(user_path, cache_key, entry)
            return False, None
            
        # Check if cache file actually exists
        cache_file = os.path.join(cache_dir, entry['file_path'])
        if not os.path.exists(cache_file):
            logger.warning(f"Cache file missing: {cache_file}")
            # Clean up metadata entry for missing file
            remove_cache_entry_from_metadata(user_path, cache_key)
            return False, None
            
        # Update access information
        entry['last_accessed'] = current_time.isoformat()
        entry['access_count'] = entry.get('access_count', 0) + 1
        save_cache_metadata(metadata_file, metadata)
        
        logger.info(f"Cache hit: key {cache_key}, file {cache_file}")
        return True, cache_file
        
    except Exception as e:
        logger.error(f"Error checking cache for key {cache_key}: {e}")
        return False, None


def remove_expired_cache_entry(user_path: str, cache_key: str, entry: Dict[str, Any]) -> None:
    """
    Remove expired cache entry including file and all symbolic links.
    
    Args:
        user_path: User's base directory path
        cache_key: Cache key to remove
        entry: Cache entry metadata
    """
    try:
        cache_dir = get_cache_dir_path(user_path)
        cache_file = os.path.join(cache_dir, entry['file_path'])
        
        # Remove actual cache file
        if os.path.exists(cache_file):
            os.remove(cache_file)
            logger.info(f"Removed expired cache file: {cache_file}")
            
        # Remove all symbolic links
        for link_path in entry.get('linked_locations', []):
            full_link_path = os.path.join(user_path, link_path)
            if os.path.islink(full_link_path):
                os.remove(full_link_path)
                logger.debug(f"Removed symbolic link: {full_link_path}")
            elif os.path.exists(full_link_path):
                # If it's not a symbolic link but exists, it might be a regular file
                # Remove it as well to clean up
                os.remove(full_link_path)
                logger.debug(f"Removed regular file: {full_link_path}")
                
        # Remove from metadata
        remove_cache_entry_from_metadata(user_path, cache_key)
        
    except Exception as e:
        logger.error(f"Error removing expired cache entry {cache_key}: {e}")


def remove_cache_entry_from_metadata(user_path: str, cache_key: str) -> None:
    """Remove a cache entry from metadata file."""
    try:
        metadata_file = get_metadata_file_path(user_path)
        metadata = load_cache_metadata(metadata_file)
        
        if cache_key in metadata['cache_entries']:
            del metadata['cache_entries'][cache_key]
            save_cache_metadata(metadata_file, metadata)
            logger.debug(f"Removed cache key {cache_key} from metadata")
            
    except Exception as e:
        logger.error(f"Error removing cache key {cache_key} from metadata: {e}")


def create_symbolic_link(cache_file_path: str, result_link_path: str) -> bool:
    """
    Create symbolic link from result location to cached file.
    
    Args:
        cache_file_path: Path to actual cached file
        result_link_path: Path where symbolic link should be created
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Create result directory if it doesn't exist
        os.makedirs(os.path.dirname(result_link_path), exist_ok=True)
        
        # Remove existing file/link if it exists
        if os.path.exists(result_link_path) or os.path.islink(result_link_path):
            os.remove(result_link_path)
            
        # Calculate relative path for the symbolic link
        rel_path = os.path.relpath(cache_file_path, os.path.dirname(result_link_path))
        
        # Create symbolic link
        os.symlink(rel_path, result_link_path)
        logger.debug(f"Created symbolic link: {result_link_path} -> {cache_file_path}")
        return True
        
    except OSError as e:
        logger.error(f"Failed to create symbolic link {result_link_path}: {e}")
        # Fallback: copy file instead of creating link
        try:
            shutil.copy2(cache_file_path, result_link_path)
            logger.info(f"Fallback: copied file to {result_link_path}")
            return True
        except Exception as copy_e:
            logger.error(f"Fallback copy also failed: {copy_e}")
            return False


def save_result_to_cache(result_file_path: str, user_path: str, cache_key: str, 
                        plugin_name: str, script_name: str, 
                        linked_location: str) -> bool:
    """
    Save visualization result to cache and create symbolic link.
    
    Args:
        result_file_path: Path to the result file to cache
        user_path: User's base directory path
        cache_key: Cache key for this result
        plugin_name: Name of the visualization plugin
        script_name: Name of the visualization script
        linked_location: Relative path where symbolic link will be created
        
    Returns:
        True if successful, False otherwise
    """
    try:
        cache_dir = get_cache_dir_path(user_path)
        metadata_file = get_metadata_file_path(user_path)
        
        # Ensure cache directory exists
        os.makedirs(cache_dir, exist_ok=True)
        
        # Determine cache file extension
        _, ext = os.path.splitext(result_file_path)
        cache_filename = f"{cache_key}{ext}"
        cache_file_path = os.path.join(cache_dir, cache_filename)
        
        # Move result file to cache (or copy if move fails)
        if os.path.exists(result_file_path):
            try:
                shutil.move(result_file_path, cache_file_path)
                logger.info(f"Moved result to cache: {cache_file_path}")
            except Exception as e:
                logger.warning(f"Failed to move file, copying instead: {e}")
                shutil.copy2(result_file_path, cache_file_path)
                os.remove(result_file_path)  # Remove original after copying
                logger.info(f"Copied result to cache: {cache_file_path}")
        else:
            logger.error(f"Result file not found: {result_file_path}")
            return False
            
        # Create symbolic link at original location
        if not create_symbolic_link(cache_file_path, result_file_path):
            logger.error(f"Failed to create symbolic link for {cache_key}")
            return False
            
        # Update metadata
        metadata = load_cache_metadata(metadata_file)
        
        # Get file size for statistics
        file_size = os.path.getsize(cache_file_path) if os.path.exists(cache_file_path) else 0
        
        metadata['cache_entries'][cache_key] = {
            'created_at': datetime.now(timezone.utc).isoformat(),
            'last_accessed': datetime.now(timezone.utc).isoformat(),
            'file_path': cache_filename,
            'file_size': file_size,
            'access_count': 1,
            'plugin_name': plugin_name,
            'script_name': script_name,
            'linked_locations': [linked_location]
        }
        
        if save_cache_metadata(metadata_file, metadata):
            logger.info(f"Saved cache entry for key {cache_key}")
            return True
        else:
            logger.error(f"Failed to save metadata for cache key {cache_key}")
            return False
            
    except Exception as e:
        logger.error(f"Error saving result to cache for key {cache_key}: {e}")
        return False


def update_cache_link_location(user_path: str, cache_key: str, linked_location: str) -> bool:
    """
    Add a new linked location to an existing cache entry.
    
    Args:
        user_path: User's base directory path
        cache_key: Cache key
        linked_location: New relative path to add to linked locations
        
    Returns:
        True if successful, False otherwise
    """
    try:
        metadata_file = get_metadata_file_path(user_path)
        metadata = load_cache_metadata(metadata_file)
        
        if cache_key in metadata['cache_entries']:
            entry = metadata['cache_entries'][cache_key]
            linked_locations = entry.get('linked_locations', [])
            
            if linked_location not in linked_locations:
                linked_locations.append(linked_location)
                entry['linked_locations'] = linked_locations
                return save_cache_metadata(metadata_file, metadata)
                
        return True  # Already exists or entry not found
        
    except Exception as e:
        logger.error(f"Error updating cache link location for key {cache_key}: {e}")
        return False


def cleanup_expired_cache(user_path: str, max_age_days: int = 7) -> int:
    """
    Clean up all expired cache entries for a user.
    
    Args:
        user_path: User's base directory path
        max_age_days: Maximum age in days before cache expires
        
    Returns:
        Number of expired entries cleaned up
    """
    try:
        metadata_file = get_metadata_file_path(user_path)
        metadata = load_cache_metadata(metadata_file)
        
        expired_keys = []
        current_time = datetime.now(timezone.utc)
        
        # Find expired entries
        for cache_key, entry in metadata['cache_entries'].items():
            created_at_str = entry.get('created_at')
            if not created_at_str:
                expired_keys.append(cache_key)
                continue
                
            try:
                created_at = datetime.fromisoformat(created_at_str.replace('Z', '+00:00'))
                age = current_time - created_at
                
                if age > timedelta(days=max_age_days):
                    expired_keys.append(cache_key)
                    
            except ValueError:
                # Invalid timestamp, mark for removal
                expired_keys.append(cache_key)
                
        # Remove expired entries
        for cache_key in expired_keys:
            entry = metadata['cache_entries'].get(cache_key, {})
            remove_expired_cache_entry(user_path, cache_key, entry)
            
        # Update last cleanup time
        metadata = load_cache_metadata(metadata_file)  # Reload after cleanup
        metadata['settings']['last_cleanup'] = current_time.isoformat()
        save_cache_metadata(metadata_file, metadata)
        
        if expired_keys:
            logger.info(f"Cleaned up {len(expired_keys)} expired cache entries for user")
            
        return len(expired_keys)
        
    except Exception as e:
        logger.error(f"Error during cache cleanup: {e}")
        return 0


def maybe_cleanup_cache(user_path: str, cleanup_probability: float = 0.1) -> int:
    """
    Probabilistically run cache cleanup.
    
    Args:
        user_path: User's base directory path
        cleanup_probability: Probability (0.0-1.0) of running cleanup
        
    Returns:
        Number of expired entries cleaned up (0 if cleanup not triggered)
    """
    if random.random() < cleanup_probability:
        logger.debug(f"Triggered probabilistic cache cleanup (p={cleanup_probability})")
        return cleanup_expired_cache(user_path)
    return 0


def get_cache_statistics(user_path: str) -> Dict[str, Any]:
    """
    Get cache statistics for a user.
    
    Args:
        user_path: User's base directory path
        
    Returns:
        Dictionary containing cache statistics
    """
    try:
        metadata_file = get_metadata_file_path(user_path)
        metadata = load_cache_metadata(metadata_file)
        
        entries = metadata['cache_entries']
        total_size = sum(entry.get('file_size', 0) for entry in entries.values())
        total_entries = len(entries)
        total_accesses = sum(entry.get('access_count', 0) for entry in entries.values())
        
        # Count by plugin
        plugin_stats = {}
        for entry in entries.values():
            plugin = entry.get('plugin_name', 'unknown')
            plugin_stats[plugin] = plugin_stats.get(plugin, 0) + 1
            
        return {
            "total_entries": total_entries,
            "total_size_bytes": total_size,
            "total_size_mb": round(total_size / (1024 * 1024), 2),
            "total_accesses": total_accesses,
            "avg_accesses_per_entry": round(total_accesses / max(total_entries, 1), 2),
            "plugin_distribution": plugin_stats,
            "last_cleanup": metadata['settings'].get('last_cleanup', 'never'),
            "max_age_days": metadata['settings'].get('max_age_days', 7)
        }
        
    except Exception as e:
        logger.error(f"Error getting cache statistics: {e}")
        return {"error": str(e)}


def clear_all_cache(user_path: str) -> bool:
    """
    Clear all cache entries for a user. Use with caution.
    
    Args:
        user_path: User's base directory path
        
    Returns:
        True if successful, False otherwise
    """
    try:
        metadata_file = get_metadata_file_path(user_path)
        metadata = load_cache_metadata(metadata_file)
        
        # Remove all cache entries
        for cache_key, entry in metadata['cache_entries'].items():
            remove_expired_cache_entry(user_path, cache_key, entry)
            
        logger.info(f"Cleared all cache entries for user")
        return True
        
    except Exception as e:
        logger.error(f"Error clearing all cache: {e}")
        return False