from enum import Enum


class PluginType(str, Enum):
    ANALYSIS = "analysis"
    VISUALIZATION = "visualization"