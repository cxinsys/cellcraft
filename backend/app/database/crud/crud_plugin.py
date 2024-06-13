from sqlalchemy.orm import Session
from app.database import models
from app.database.schemas import plugin

def get_plugin(db: Session, plugin_id: int):
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