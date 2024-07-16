from sqlalchemy.orm import Session
from app.database import models
from app.database.schemas import plugin
from fastapi import HTTPException

def get_plugin_by_id(db: Session, plugin_id: int):
    return db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()

def get_plugins(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Plugin).offset(skip).limit(limit).all()

def create_plugin(db: Session, plugin: plugin.PluginCreate):
    db_plugin = models.Plugin(**plugin.dict())
    db.add(db_plugin)
    db.commit()
    db.refresh(db_plugin)
    return db_plugin

def update_plugin(db: Session, plugin: plugin.PluginUpdate, plugin_id: int):
    db_plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()
    for key, value in plugin.dict().items():
        setattr(db_plugin, key, value)
    db.commit()
    db.refresh(db_plugin)
    return db_plugin

def delete_plugin(db: Session, plugin_id: int):
    db_plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()
    db.delete(db_plugin)
    db.commit()
    return db_plugin

def associate_user_plugin(db: Session, user_id: int, plugin_id: int):
    user = db.query(models.User).filter(models.User.id == user_id).first()
    plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()

    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    if not plugin:
        raise HTTPException(status_code=404, detail="Plugin not found")
    
    if plugin not in user.plugins:
        user.plugins.append(plugin)
        db.commit()
        return {"message": "Plugin associated with user successfully"}
    else:
        raise HTTPException(status_code=400, detail="Plugin already associated with user")

def dissociate_user_plugin(db: Session, user_id: int, plugin_id: int):
    user = db.query(models.User).filter(models.User.id == user_id).first()
    plugin = db.query(models.Plugin).filter(models.Plugin.id == plugin_id).first()
    
    if not user:
        raise HTTPException(status_code=404, detail="User not found")
    if not plugin:
        raise HTTPException(status_code=404, detail="Plugin not found")
    
    if plugin in user.plugins:
        user.plugins.remove(plugin)
        db.commit()
        return {"message": "Plugin dissociated from user successfully"}
    else:
        raise HTTPException(status_code=400, detail="Plugin not associated with user")
