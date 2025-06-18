from bag_pipe.io import parse_cigar
from bag_pipe.io import TagMapRow

def test_parse_cigar_spliced():
    """CIGAR with an intron returns two exons + one junction."""
    s, e, js, je = parse_cigar(1000, "76M225N74M")
    print(s, e, js, je)
    assert s == [1000, 1301]
    assert e == [1075, 1374]
    assert js == [1075] and je == [1301]

def test_parse_cigar_unmapped():
    """Unmapped read (`*`) yields empty lists."""
    assert parse_cigar(42, "*") == ([], [], [], [])

def test_tagmap_splice_keys(tagmap_row):
    """`splice_keys` must format chr.end.start correctly."""
    assert tagmap_row.splice_keys() == {"chr1.1301.1075"}
