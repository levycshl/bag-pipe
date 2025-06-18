"""Stats module for bagpipe"""

from collections import Counter
from typing import Dict, List, Tuple

def count_varietal_tags(tags: List[str]) -> Dict[str, int]:
    """Count occurrences of varietal tags"""
    return dict(Counter(tags))

def calculate_mapping_stats(mapping_qualities: List[int]) -> Dict[str, float]:
    """Calculate basic statistics for mapping qualities"""
    if not mapping_qualities:
        return {"mean": 0, "median": 0, "min": 0, "max": 0}
    
    return {
        "mean": sum(mapping_qualities) / len(mapping_qualities),
        "median": sorted(mapping_qualities)[len(mapping_qualities) // 2],
        "min": min(mapping_qualities),
        "max": max(mapping_qualities)
    }