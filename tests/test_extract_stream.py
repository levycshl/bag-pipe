# tests/test_extract_stream.py
from bag_pipe.io import write_feature_tsv
from pathlib import Path

def test_write_feature_tsv(tmp_path):
    sam = Path(__file__).parent.parent / "data" / "test_data" / "hisat2pe_sorted.sam"
    out = tmp_path / "tiny.tsv"
    rows = write_feature_tsv(sam, out)
    assert out.exists() and rows > 0
    assert sum(1 for _ in out.open()) == rows + 1   # +1 for header
