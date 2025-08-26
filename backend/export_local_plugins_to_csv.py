#!/usr/bin/env python3
"""
Export local plugins from database to CSV format for official plugin repository.
This script converts local plugins to the format expected by the plugin synchronization system.
Output is saved to /app/plugin/official/plugins.csv (replaces legacy plugin_initialization.csv).
"""
import os
import sys
import csv
import json
import subprocess
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Add the backend directory to Python path
sys.path.append('/app')

from app.common.config import settings

def export_local_plugins_to_csv(output_path: str):
    """
    Export local plugins from database to CSV file with all required columns.
    
    Args:
        output_path (str): Path to save the CSV file
    """
    # Create database connection
    engine = create_engine(settings.SQLALCHEMY_DATABASE_URI, pool_pre_ping=True)
    SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    session = SessionLocal()
    
    # Import models after creating the connection to avoid circular import
    from app.database import models
    
    try:
        # Query all local plugins
        local_plugins = session.query(models.Plugin).filter_by(source='local').order_by(models.Plugin.name).all()
        
        if not local_plugins:
            print("No local plugins found in database")
            return
        
        print(f"Found {len(local_plugins)} local plugins")
        
        # Open CSV file for writing
        with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write header with all required columns
            writer.writerow([
                'name', 'description', 'author', 'plugin_path', 'plugin_type',
                'dependencies', 'drawflow', 'rules', 'use_gpu', 'source',
                'is_editable', 'version', 'submodule_path'
            ])
            
            # Write plugin data
            for plugin in local_plugins:
                # Modify plugin path from local to official format
                # From: ./plugin/local/PluginName/ to ./plugin/official/PluginName/
                modified_path = plugin.plugin_path.replace('/local/', '/official/')
                
                # Convert JSONB fields to JSON strings
                dependencies_json = json.dumps(plugin.dependencies) if plugin.dependencies else '{}'
                drawflow_json = json.dumps(plugin.drawflow) if plugin.drawflow else '{}'
                rules_json = json.dumps(plugin.rules) if plugin.rules else '{}'
                
                # Convert plugin_type enum to string
                plugin_type_str = plugin.plugin_type.value if plugin.plugin_type else ''
                
                # Generate submodule_path in format: cellcraft-plugin/PluginName
                submodule_path = f"cellcraft-plugin/{plugin.name}"
                
                writer.writerow([
                    plugin.name,
                    plugin.description,
                    "admin",  # Change author to admin
                    modified_path,
                    plugin_type_str,
                    dependencies_json,
                    drawflow_json,
                    rules_json,
                    plugin.use_gpu,
                    "official",  # Change source to official
                    False,  # is_editable = False for official plugins
                    "1.0",  # Set version to 1.0
                    submodule_path
                ])
                
                print(f"  - Exported: {plugin.name}")
        
        print(f"\nSuccessfully exported {len(local_plugins)} plugins to {output_path}")
        print("\nNote: The following modifications were made:")
        print("  - Plugin paths: changed from '/local/' to '/official/'")
        print("  - Author: set to 'admin'")
        print("  - Source: set to 'official'")
        print("  - is_editable: set to False")
        print("  - Version: set to '1.0'")
        print("  - submodule_path: set to 'cellcraft-plugin/PluginName'")
        print(f"\nThis file replaces the legacy plugin_initialization.csv system")
        print(f"and can be used directly with the plugin synchronization system.")
        
    except Exception as e:
        print(f"Error exporting plugins: {e}")
        raise
    finally:
        session.close()

def main():
    # Output path - now uses official plugins location
    output_path = "/app/plugin/official/plugins.csv"
    
    # Check if we're running inside the backend container
    if not os.path.exists('/app'):
        print("Error: This script must be run inside the backend container")
        print("Use: docker compose -f docker-compose.dev.yml exec backend python /app/export_local_plugins_to_csv.py")
        return 1
    
    # Create backup of existing file if it exists
    if os.path.exists(output_path):
        backup_path = output_path + '.backup'
        print(f"Backing up existing file to {backup_path}")
        import shutil
        shutil.copy2(output_path, backup_path)
    
    # Export plugins
    export_local_plugins_to_csv(output_path)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())