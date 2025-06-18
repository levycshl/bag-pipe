from bag_pipe.tags import hamming, split_tag, rollup_tags  # tags impl :contentReference[oaicite:3]{index=3}

def test_hamming_basic():
    assert hamming("AAA", "AAA") == 0
    assert hamming("AAA", "AAT") == 1

def test_split_tag_integrity():
    parts = split_tag("ABCDEFGHIJ", k=2)  # 3 blocks expected
    assert len(parts) == 3
    assert "".join(parts) == "ABCDEFGHIJ"

def test_rollup_simple_case():
    keep, mapping, rejects = rollup_tags(
        [("AAA", 100), ("AAB", 5), ("CCC", 90)], k=1, threshold=0.2
    )
    assert keep == {"AAA", "CCC"}
    assert mapping == {"AAB": "AAA"}
    assert ("AAB", "AAA") in rejects
