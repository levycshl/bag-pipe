"""Configuration handling for bagpipe"""

import os
import yaml
from typing import Dict, Any

def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file"""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def get_default_config() -> Dict[str, Any]:
    """Get default configuration"""
    return {
        "output_dir": "./output",
        "threads": 4,
        "samples": [],
        "reference": {
            "refflat": None
        },
        "parameters": {
            "min_mapping_quality": 10,
            "min_read_length": 30
        }
    }