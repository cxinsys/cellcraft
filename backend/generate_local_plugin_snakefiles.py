#!/usr/bin/env python3
"""
Generate Snakefiles for all plugins in the local plugin directory.
"""
import os
import json
import sys
sys.path.append('/app')
from app.common.utils.plugin_utils import generate_snakemake_code

LOCAL_PLUGINS_DIR = '/app/plugin/local'

def main():
    """Generate Snakefiles for all local plugins."""
    # Get all plugin directories in the local plugin folder
    if not os.path.exists(LOCAL_PLUGINS_DIR):
        print(f'Local plugins directory not found: {LOCAL_PLUGINS_DIR}')
        return
    
    plugin_dirs = [
        d for d in os.listdir(LOCAL_PLUGINS_DIR)
        if os.path.isdir(os.path.join(LOCAL_PLUGINS_DIR, d)) and not d.startswith('.')
    ]
    
    if not plugin_dirs:
        print('No plugins found in the local plugins directory.')
        return
    
    print(f'Found {len(plugin_dirs)} plugins in local directory:')
    for plugin_name in sorted(plugin_dirs):
        print(f'  - {plugin_name}')
    
    print('\nGenerating Snakefiles...\n')
    
    success_count = 0
    error_count = 0
    
    for plugin_name in sorted(plugin_dirs):
        plugin_path = os.path.join(LOCAL_PLUGINS_DIR, plugin_name)
        metadata_path = os.path.join(plugin_path, 'metadata.json')
        
        print(f'Processing plugin: {plugin_name}')
        
        # Check if metadata.json exists
        if not os.path.exists(metadata_path):
            print(f'  ❌ Error: metadata.json not found in {plugin_path}')
            error_count += 1
            continue
        
        try:
            # Load metadata.json
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
            
            # Check if rules exist in metadata
            if 'rules' not in metadata:
                print(f'  ❌ Error: rules key not found in metadata.json')
                error_count += 1
                continue
            
            rules_data = metadata['rules']
            
            # Generate Snakefile
            generate_snakemake_code(
                rules_data=rules_data,
                output_folder_path=plugin_path,
                plugin_name=plugin_name
            )
            
            # Verify Snakefile was created
            snakefile_path = os.path.join(plugin_path, 'Snakefile')
            visualization_snakefile_path = os.path.join(plugin_path, 'visualization_Snakefile')
            
            files_created = []
            if os.path.exists(snakefile_path):
                files_created.append('Snakefile')
            if os.path.exists(visualization_snakefile_path):
                files_created.append('visualization_Snakefile')
            
            if files_created:
                print(f"  ✅ Success: Created {', '.join(files_created)}")
                success_count += 1
            else:
                print(f'  ⚠️  Warning: No Snakefiles were created')
                error_count += 1
            
        except json.JSONDecodeError as e:
            print(f'  ❌ Error: Failed to parse metadata.json - {str(e)}')
            error_count += 1
        except Exception as e:
            print(f'  ❌ Error: {str(e)}')
            error_count += 1
        
        print()  # Empty line for readability
    
    # Summary
    print('='*50)
    print(f'Summary: {success_count} successful, {error_count} errors')
    print('='*50)

if __name__ == "__main__":
    main()