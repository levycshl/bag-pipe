from bag_pipe.stats import calculate_mapping_stats  # stats impl :contentReference[oaicite:4]{index=4}

def test_mapping_stats_values():
    vals = [10, 20, 30]
    stats = calculate_mapping_stats(vals)
    assert stats == {"mean": 20, "median": 20, "min": 10, "max": 30}

def test_mapping_stats_empty():
    assert calculate_mapping_stats([]) == {"mean": 0, "median": 0, "min": 0, "max": 0}
