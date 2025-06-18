"""
End-to-end test: run ReadParser on a miniature FASTQ pair and verify that
clipped outputs are produced and have the expected read count.
"""
from pathlib import Path
import gzip
import yaml

from bag_pipe.config import load_config
from read_break.io import FastqReader, FastqWriter
from read_break.parser import ReadParser, read_clip_and_write


def _count_fastq_reads(path: Path) -> int:
    """Return number of reads in a (possibly gzipped) FASTQ file."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def test_read_parser_end_to_end(tmp_path):
    # --- Locate bundled test data ------------------------------------------
    data_dir = Path(__file__).parent.parent / "data"
    confs_dir = Path(__file__).parent.parent / "configs"
    test_data_dir = data_dir / "test_data"

    bp_cfg = load_config(confs_dir / "bagpipe_config.yaml")
    rb_cfg = yaml.safe_load((confs_dir / "bagpipe_2N.yaml").read_text())

    sample = bp_cfg["samples"][0]
    r1 = sample["fastq_r1"]
    r2 = sample["fastq_r2"]

    # --- Instantiate parser & I/O ------------------------------------------
    parser = ReadParser(
        pipeline_cfg={"pipeline": rb_cfg["pipeline"]},
        globals_cfg=rb_cfg["params"],
        globals_namespace="params",
        base_dir=str(data_dir),
    )

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    reader = FastqReader(str(r1), str(r2), trim_tail=True)
    writer = FastqWriter(str(out_dir), sample["id"])

    # --- Run clipping -------------------------------------------------------
    read_clip_and_write(
        reader,
        parser,
        writer,
        default_start_r1=parser.globals.get("start_r1", 0),
        default_start_r2=parser.globals.get("start_r2", 0),
    )
    writer.close()

    # --- Assertions ---------------------------------------------------------
    out_r1 = out_dir / f"{sample['id']}.R1.fastq.gz"
    out_r2 = out_dir / f"{sample['id']}.R2.fastq.gz"

    assert out_r1.exists() and out_r2.exists(), "Output FASTQs not created"
    r1 = Path(r1)
    # Make sure we still have at least one read per input read (trimmed counts equal is OK)
    in_reads = _count_fastq_reads(r1)
    out_reads = _count_fastq_reads(out_r1)
    assert 0 < out_reads <= in_reads, "Unexpected clipped-read count"
