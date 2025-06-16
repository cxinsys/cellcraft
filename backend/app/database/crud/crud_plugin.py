from sqlalchemy.orm import Session
from sqlalchemy.exc import SQLAlchemyError
from app.database import models
from app.database.schemas import plugin
from fastapi import HTTPException

def get_plugin_by_id(db: Session, plugin_id: int):
    """
    Get a plugin by its ID.
    
    Args:
        db (Session): Database session
        plugin_id (int): Plugin ID
    
    Returns:
        Plugin: Plugin object if found
    
    Raises:
        HTTPException: If plugin not found
    """
    try:
        plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()
        if not plugin:
            raise HTTPException(status_code=404, detail=f"Plugin with id {plugin_id} not found")
        return plugin
    except SQLAlchemyError as e:
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")

def get_plugin_by_name(db: Session, name: str):
    """
    Get a plugin by its name.
    
    Args:
        db (Session): Database session
        name (str): Plugin name
    
    Returns:
        Plugin: Plugin object if found
    
    Raises:
        HTTPException: If plugin not found
    """
    try:
        plugin = db.query(models.Plugin).filter(models.Plugin.name == name).first()
        if not plugin:
            raise HTTPException(status_code=404, detail=f"Plugin with name {name} not found")
        return plugin
    except SQLAlchemyError as e:
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")

def get_plugins(db: Session, skip: int = 0, limit: int = 100):
    """
    Get all plugins with pagination.
    
    Args:
        db (Session): Database session
        skip (int): Number of records to skip
        limit (int): Maximum number of records to return
    
    Returns:
        List[Plugin]: List of plugins
    
    Raises:
        HTTPException: If database error occurs
    """
    try:
        return db.query(models.Plugin).offset(skip).limit(limit).all()
    except SQLAlchemyError as e:
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")

def create_plugin(db: Session, plugin: plugin.PluginUpdate):
    """
    Create a new plugin.
    
    Args:
        db (Session): Database session
        plugin (PluginUpdate): Plugin data
    
    Returns:
        Plugin: Created plugin object
    
    Raises:
        HTTPException: If plugin creation fails
    """
    try:
        # Check if plugin with same name already exists
        existing_plugin = db.query(models.Plugin).filter(models.Plugin.name == plugin.name).first()
        if existing_plugin:
            raise HTTPException(status_code=400, detail=f"Plugin with name {plugin.name} already exists")

        plugin_data = plugin.dict(exclude={"reference_folders"})
        db_plugin = models.Plugin(**plugin_data)
        db.add(db_plugin)
        db.commit()
        db.refresh(db_plugin)
        return db_plugin
    except SQLAlchemyError as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Error creating plugin: {str(e)}")

def update_plugin(db: Session, plugin: plugin.PluginUpdate, plugin_id: int):
    """
    Update an existing plugin.
    
    Args:
        db (Session): Database session
        plugin (PluginUpdate): Updated plugin data
        plugin_id (int): ID of plugin to update
    
    Returns:
        Plugin: Updated plugin object
    
    Raises:
        HTTPException: If plugin update fails
    """
    try:
        db_plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()
        if not db_plugin:
            raise HTTPException(status_code=404, detail=f"Plugin with id {plugin_id} not found")

        # Check if name is being changed and if new name already exists
        if plugin.name and plugin.name != db_plugin.name:
            existing_plugin = db.query(models.Plugin).filter(models.Plugin.name == plugin.name).first()
            if existing_plugin:
                raise HTTPException(status_code=400, detail=f"Plugin with name {plugin.name} already exists")

        update_data = plugin.dict(exclude_unset=True)
        for key, value in update_data.items():
            setattr(db_plugin, key, value)

        db.commit()
        db.refresh(db_plugin)
        return db_plugin
    except SQLAlchemyError as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Error updating plugin: {str(e)}")

def delete_plugin(db: Session, plugin_id: int):
    """
    Delete a plugin.
    
    Args:
        db (Session): Database session
        plugin_id (int): ID of plugin to delete
    
    Returns:
        Plugin: Deleted plugin object
    
    Raises:
        HTTPException: If plugin deletion fails
    """
    try:
        db_plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()
        if not db_plugin:
            raise HTTPException(status_code=404, detail=f"Plugin with id {plugin_id} not found")

        db.delete(db_plugin)
        db.commit()
        return db_plugin
    except SQLAlchemyError as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Error deleting plugin: {str(e)}")

def associate_user_plugin(db: Session, user_id: int, plugin_id: int):
    """
    Associate a user with a plugin.
    
    Args:
        db (Session): Database session
        user_id (int): User ID
        plugin_id (int): Plugin ID
    
    Returns:
        dict: Success message
    
    Raises:
        HTTPException: If association fails
    """
    try:
        user = db.query(models.User).filter(models.User.id == user_id).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()
        if not plugin:
            raise HTTPException(status_code=404, detail="Plugin not found")
        
        if plugin in user.plugins:
            raise HTTPException(status_code=400, detail="Plugin already associated with user")

        user.plugins.append(plugin)
        db.commit()
        return {"message": "Plugin associated with user successfully"}
    except SQLAlchemyError as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=str(e))

def dissociate_user_plugin(db: Session, user_id: int, plugin_id: int):
    """
    Dissociate a user from a plugin.
    
    Args:
        db (Session): Database session
        user_id (int): User ID
        plugin_id (int): Plugin ID
    
    Returns:
        dict: Success message
    
    Raises:
        HTTPException: If dissociation fails
    """
    try:
        user = db.query(models.User).filter(models.User.id == user_id).first()
        if not user:
            raise HTTPException(status_code=404, detail="User not found")

        plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()
        if not plugin:
            raise HTTPException(status_code=404, detail="Plugin not found")
        
        if plugin not in user.plugins:
            raise HTTPException(status_code=400, detail="Plugin not associated with user")

        user.plugins.remove(plugin)
        db.commit()
        return {"message": "Plugin dissociated from user successfully"}
    except SQLAlchemyError as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=str(e))
