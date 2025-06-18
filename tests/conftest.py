"""
Shared fixtures for bag-pipe tests.
Creates a TagMapRow with one splice junction so splice_key logic can be reused.
"""
import pytest
from bag_pipe.io import parse_cigar, TagMapRow

@pytest.fixture
def tagmap_row():
    start = 1000
    cigar = "76M225N74M"
    ex_starts, ex_ends, sj_starts, sj_ends = parse_cigar(start, cigar)
    return TagMapRow(
        read_id="r1",
        barcode="BC01",
        vt="VT1",
        chrom="chr1",
        start=ex_starts[0],
        end=ex_ends[-1],
        flag=0,
        mapq=60,
        cigar=cigar,
        block_starts=ex_starts,
        block_ends=ex_ends,
        sj_starts=sj_starts,
        sj_ends=sj_ends,
    )
