"""
Command-line entry points for Bag-Pipe.
Run “bag_pipe --help” to see the list of sub-commands.
"""
from pathlib import Path
import click
import yaml
from itertools import groupby
from read_break.io import FastqReader, FastqWriter
from read_break.parser import ReadParser, read_clip_and_write
from bag_pipe.io import write_feature_tsv
import gzip
from bag_pipe.io          import tag_and_map_stream, load_refflat
from bag_pipe.gene_models import GeneModel, summarise_vt
from bag_pipe.tags        import rollup_tags


# ─────────────────────────────────────────────────────────────
#  CLI root
# ─────────────────────────────────────────────────────────────
@click.group()
def main():
    """Bag-Pipe command-line interface."""
    pass


# ─────────────────────────────────────────────────────────────
#  bag_pipe trim
# ─────────────────────────────────────────────────────────────
@main.command("trim")
@click.option("--r1", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Path to R1 FASTQ (gz or plain)")
@click.option("--r2", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Path to R2 FASTQ (gz or plain)")
@click.option("--parse-config", "-c", required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="Read-Breaker YAML (e.g. bagpipe_2N.yaml)")
@click.option("--sample-id", "-s", required=True,
              help="Sample ID used to name output files")
@click.option("--outdir", "-o", required=True, type=click.Path(file_okay=False),
              help="Directory to write clipped FASTQs")
@click.option("--trim-tail/--no-trim-tail", default=True, show_default=True,
              help="Pass trim_tail flag to FastqReader")
def trim_cmd(r1, r2, parse_config, sample_id, outdir, trim_tail):
    """Clip + barcode-tag one paired-end sample."""
    cfg_path = Path(parse_config)
    rb_cfg = yaml.safe_load(cfg_path.read_text())

    parser = ReadParser(
        pipeline_cfg={"pipeline": rb_cfg["pipeline"]},
        globals_cfg=rb_cfg["params"],
        globals_namespace="params",
        base_dir=str(cfg_path.parent),
    )

    out_root = Path(outdir)
    out_root.mkdir(parents=True, exist_ok=True)

    reader = FastqReader(r1, r2, trim_tail=trim_tail)
    writer = FastqWriter(out_root, sample_id)

    read_clip_and_write(
        reader, parser, writer,
        default_start_r1=parser.globals.get("start_r1", 0),
        default_start_r2=parser.globals.get("start_r2", 0),
    )
    writer.close()
    click.echo(f"✅  wrote {sample_id}.R1.fastq.gz / .R2.fastq.gz to {out_root}")


# ─────────────────────────────────────────────────────────────
#  bag_pipe extract
# ─────────────────────────────────────────────────────────────
@main.command("extract")
@click.argument("sam", type=click.Path(exists=True, dir_okay=False))
@click.argument("tsv", type=click.Path(dir_okay=False))
@click.option("--min-mapq", default=1, show_default=True, type=int,
              help="Drop alignments with MAPQ < MIN_MAPQ")
def extract_cmd(sam, tsv, min_mapq):
    """Stream SAM/BAM → TSV feature table."""
    n = write_feature_tsv(sam, tsv, min_mapq=min_mapq)
    click.echo(f"✅  wrote {n:,} rows → {Path(tsv).name}")



# ─────────────────────────────────────────────────────────────
#  bag_pipe summarise  – stream VT summaries per barcode
# ─────────────────────────────────────────────────────────────

@main.command("summarise")
@click.option("--chrom",        "-c", required=True,
              help="Chromosome name (e.g. chr21) used to filter refFlat")
@click.option("--tsv-in",       "-i", "tsv_in", required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="Input TSV (output of bag_pipe extract). Can be .gz")
@click.option("--ref-flat",     "-r", "ref_flat", required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="refFlat file for gene-model metrics")
@click.option("--tsv-out",      "-o", "tsv_out", required=True,
              type=click.Path(dir_okay=False),
              help="Output TSV (gz if filename ends with .gz)")
@click.option("--k",            default=1, show_default=True, type=int,
              help="Maximum Hamming distance for VT roll-up")
@click.option("--threshold",    default=0.10, show_default=True, type=float,
              help="Relative-abundance threshold for roll-up")
@click.option("--header/--no-header", default=True, show_default=True,
              help="Write header line to output")
def summarise_cmd(chrom, tsv_in, ref_flat, tsv_out,
                  k, threshold, header):
    """
    Stream barcode-level VT summaries into a TSV.
    Runs constant-memory, writing each VTSummary row as soon as it is ready.
    """
    # ---------- set up iterators & models ----------
    rows_iter = tag_and_map_stream(tsv_in)
    model     = GeneModel(load_refflat(ref_flat, chrom=chrom))

    # auto gzip if user gave *.gz
    out_path  = Path(tsv_out)
    opener    = gzip.open if out_path.suffix == ".gz" else open

    with opener(out_path, "wt") as fout:
        header_written = False

        # group rows by barcode (requires TSV already sorted by barcode)
        for bc, bc_iter in groupby(rows_iter, key=lambda r: r.barcode):
            # bucket by VT
            vt_buckets = {vt: list(vt_iter)
                          for vt, vt_iter in groupby(bc_iter, key=lambda r: r.vt)}

            # roll-up
            tag_counts = [(vt, len(rows)) for vt, rows in vt_buckets.items()]
            keepers, _, _ = rollup_tags(tag_counts, k=k, threshold=threshold)

            # summarise each surviving VT
            for vt in keepers:
                summary = summarise_vt(vt_buckets[vt], model)   # VTSummary dataclass
                # first time: write header
                if header and not header_written:
                    fout.write("\t".join(summary.header()) + "\n")
                    header_written = True

                fout.write(str(summary) + "\n")

    click.echo(f"✅  Wrote summaries → {out_path}")
