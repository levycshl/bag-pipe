"""IO module for bagpipe"""
import re
import os
import pandas as pd
from typing import Generator, List, Set, TextIO, Optional, NamedTuple, Tuple, Iterator
import sys
import gzip
import pathlib
from dataclasses import dataclass

REFCOLS = (
    "gene", "transcript", "chrom", "strand",
    "tx_start", "tx_end", "cds_start", "cds_end",
    "exon_count", "exon_starts", "exon_ends"
)


def load_refflat(path: str, chrom: str | None = None) -> pd.DataFrame:
    """
    Read UCSC refFlat into a tidy DataFrame.

    * Converts 0-based starts ➜ 1-based closed.
    * Splits comma lists into Python lists of `int`.
    * Optionally filters to a single chromosome.
    """
    df = pd.read_csv(
        path, sep="\t", header=None, names=REFCOLS, usecols=range(len(REFCOLS)),
        dtype={"chrom": "category", "strand": "category"}
    )
    if chrom:
        df = df.query("chrom == @chrom").copy()

    # 1-based coordinate fix
    df["tx_start"]  += 1
    df["cds_start"] += 1

    def _split0(s: str, adjust: bool = False) -> list[int]:
        items = [int(x) for x in s.rstrip(",").split(",")]
        if adjust:  # exon_starts need +1
            items = [i + 1 for i in items]
        return items

    df["exon_starts"] = df["exon_starts"].map(lambda s: _split0(s, adjust=True))
    df["exon_ends"]   = df["exon_ends"].map(_split0)

    return df.reset_index(drop=True)

def open_fastq(filepath: str) -> TextIO:
    """Open a FASTQ file (gzipped or not)"""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')




# ──────────────────────────────────────────────────────────────────────
#  1. Lightweight alignment record
# ──────────────────────────────────────────────────────────────────────
class SamAlignment(NamedTuple):
    """
    Minimal view of a SAM alignment.

    Only the fields required downstream are kept; loading more would just
    waste memory.

    Attributes
    ----------
    qname : str
        Query / read name (field 0).
    flag : int
        Bit-wise FLAG (field 1).
    rname : str
        Reference sequence name (field 2).
    pos : int
        1-based leftmost reference coordinate (field 3).
    mapq : int
        Mapping quality (field 4).
    cigar : str
        Raw CIGAR string (field 5).
    """
    qname: str
    flag: int
    rname: str
    pos: int
    mapq: int
    cigar: str

# ──────────────────────────────────────────────────────────────────────
#  2. SAMReader – text-only, streaming, Windows-safe
# ──────────────────────────────────────────────────────────────────────
class SAMReader:
    """
    Iterate over SAM text and yield :class:`SamAlignment` objects.

    The constructor accepts:
      • a *file-like* object (already open)  – or –
      • a *path / "-"* ("-" = stdin)

    Examples
    --------
    > import sys
    > for aln in SAMReader(sys.stdin):
      ...     print(aln.qname, aln.cigar)
    """
    def __init__(self, handle: TextIO | str):
        if isinstance(handle, os.PathLike):
            handle = os.fspath(handle)
        if isinstance(handle, str):
            self._fh = sys.stdin if handle == "-" else open(handle, "rt")
            self._close_when_done = handle != "-"
        else:
            self._fh = handle
            self._close_when_done = False

    def __iter__(self) -> Iterator[SamAlignment]:
        """
        Skip @header lines and yield `SamAlignment` tuples.

        Using split is ~3× faster than csv.reader for this fixed-width format
        and avoids pulling in an extra dependency.
        """
        for line in self._fh:
            if line.startswith("@"):
                continue                       # header
            fields = line.rstrip("\n").split("\t")

            # Only the six fields we care about – everything else is ignored
            yield SamAlignment(
                qname=fields[0],
                flag=int(fields[1]),
                rname=fields[2],
                pos=int(fields[3]),
                mapq=int(fields[4]),
                cigar=fields[5],
            )

    # context-manager sugar so callers can `with SAMReader(path) as rdr:`
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._close_when_done:
            self._fh.close()


# ──────────────────────────────────────────────────────────────────────
#  3. parse_cigar – reference-space bookkeeping
# ──────────────────────────────────────────────────────────────────────
_CIGAR_SPLIT = re.compile(r"(\d+)")   # captures numbers and keeps operators

def parse_cigar(
    ref_start: int,
    cigar: str,
) -> Tuple[List[int], List[int], List[int], List[int]]:
    """
    Expand a CIGAR string into exon and junction coordinates.

    Parameters
    ----------
    ref_start
        1-based reference start (SAM POS).
    cigar
        Raw CIGAR string (e.g. ``76M225N74M``).

    Returns
    -------
    exon_starts, exon_ends, junc_starts, junc_ends
        • Every list is 1-based and inclusive.
        • Empty lists mean the read is unmapped (``*`` CIGAR) or contains no
          skipped regions (introns).

    Notes
    -----
    The algorithm walks the reference coordinate system once:

    * **M / = / X / D**  – consume *n* reference bases.
    * **N**              – end current exon, record an intron, then skip *n*.
    * **I / S / H / P**  – consume query only → no change in reference.
    """
    if cigar == "*":                       # unmapped read
        return [], [], [], []
    ## odd for lengths, even for operations
    lengths = list(map(int, _CIGAR_SPLIT.split(cigar)[1::2]))
    ops     = _CIGAR_SPLIT.split(cigar)[2::2]
    print("CIGAR:", cigar, "→", lengths, ops)
    ## called exons, but really "intervals of coverage"
    exon_starts: List[int] = []
    exon_ends:   List[int] = []
    junc_starts: List[int] = []
    junc_ends:   List[int] = []

    ref_cursor = ref_start                # leftmost base (1-based)
    exon_start = ref_cursor
    ## SOFT CLIP is handled correctly by the CIGAR parser
    for length, op in zip(lengths, ops):
        if op in ("M", "=", "X", "D"):
            ref_cursor += length
        elif op == "N":                   # intron: close exon, open next
            exon_starts.append(exon_start)
            exon_ends.append(ref_cursor - 1)

            junc_starts.append(ref_cursor - 1)
            ref_cursor += length          # jump over intron
            junc_ends.append(ref_cursor)

            exon_start = ref_cursor       # new exon begins here
        # I / S / H / P: reference coordinate stays unchanged

    # add final exon
    exon_starts.append(exon_start)
    exon_ends.append(ref_cursor - 1)

    return exon_starts, exon_ends, junc_starts, junc_ends



# ──────────────────────────────────────────────────────────────────────
#  4. extract_alignment_features – per-alignment summary
# ──────────────────────────────────────────────────────────────────────
def extract_alignment_features(
    aln: SamAlignment,
    min_mapq: int = 1,
) -> Optional[Tuple[str, ...]]:
    """
    Convert one :class:`SamAlignment` into the 13-column summary row.

    A ``None`` return means *skip this alignment*.

    Filtering rules
    ---------------
    * ``mapq < min_mapq``        → drop
    * ``flag > 255`` (secondary/supplementary) → drop
    * malformed read name        → drop
    * ``*`` CIGAR (unmapped)     → drop
    """
    # ── 1. quick filters ──────────────────────────────────────────────
    if aln.mapq < min_mapq or aln.flag > 255:
        return None

    # ── 2. parse read name → bead tag & varietal tag ─────────────────
    #     Expected pattern:  <root>/<pair>_<bead>_<vt>
    #                        e.g.  A00327:12:H7Y5VDRXX/1_BC01_VT17
    try:
        read_root, bead_tag, varietal_tag = aln.qname.split("_")
        read_root = read_root.split("/")[0]     # strip /1 or /2
    except ValueError:
        # Unexpected naming scheme → safer to skip than mis-annotate
        return None

    # ── 3. derive exon & junction coordinates ────────────────────────
    exon_starts, exon_ends, j_start, j_end = parse_cigar(aln.pos, aln.cigar)
    if not exon_starts:                      # unmapped or weird
        return None

    # ── 4. assemble final tuple (as strings) ─────────────────────────
    return (
        read_root,
        bead_tag,
        varietal_tag,
        aln.rname,
        str(exon_starts[0]),
        str(exon_ends[-1]),
        str(aln.flag),
        str(aln.mapq),
        aln.cigar,
        ",".join(map(str, exon_starts)),
        ",".join(map(str, exon_ends)),
        ",".join(map(str, j_start)),
        ",".join(map(str, j_end)),
    )





TSV_COLUMNS = (
    "read_id", "barcode", "vt", "chrom", "start", "end",
    "flag", "mapq", "cigar",
    "cigar_block_starts", "cigar_block_ends",
    "splice_starts", "splice_ends"
)

def _open(path: str | pathlib.Path | None) -> TextIO:
    """
    Return a text handle for 'path'.
    * '-' means stdin.
    * '.gz' ends open in text mode via gzip.
    """
    if path in (None, '-', ''):
        return sys.stdin
    path = pathlib.Path(path)
    if path.suffix == '.gz':
        return gzip.open(path, 'rt')
    return path.open()




# ── simple text-handle opener ─────────────────────────────────────────
def _open(path: str | pathlib.Path | None) -> TextIO:
    if path in (None, '-', ''):
        return sys.stdin
    path = pathlib.Path(path)
    if path.suffix == '.gz':
        return gzip.open(path, 'rt')
    return path.open()


# ── dataclass to describe one alignment-summary row ──────────────────
@dataclass(slots=True)
class TagMapRow:
    read_id: str
    barcode: str
    vt: str
    chrom: str
    start: int
    end: int
    flag: int
    mapq: int
    cigar: str
    block_starts: List[int]
    block_ends: List[int]
    sj_starts: List[int]
    sj_ends: List[int]

    # convenience: union of CIGAR blocks as list[(start,end)]
    @property
    def cigar_intervals(self) -> List[tuple[int, int]]:
        return list(zip(self.block_starts, self.block_ends))

    # convenience: splice-junction keys  "chr.end.start"
    def splice_keys(self) -> Set[str]:
        if not self.sj_starts:
            return set()
        chrom = self.chrom
        return {
            f"{chrom}.{e}.{s}"
            for e, s in zip(self.sj_ends, self.sj_starts)
        }


# ── parsing helpers ──────────────────────────────────────────────────
def _split_ints(cell: str) -> List[int]:
    if cell == '' or cell == ',':
        return []
    return [int(x) for x in cell.rstrip(',').split(',')]



def tag_and_map_stream(
    tsv_path: str | pathlib.Path | None,
    *,
    sample_whitelist: Set[str] | None = None
) -> Generator[TagMapRow, None, None]:
    """
    Stream TagMapRow objects from a hisat2-sorted TSV.

    Filters applied (same as legacy script)
      • barcode not in sample_whitelist
      • FLAG > 255   (secondary / supplementary)
      • chrom starts with 'chrM', contains '_', or empty
    """
    with _open(tsv_path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')

            bc    = f[1]
            flag  = int(f[6])
            chrom = f[3]

            if sample_whitelist and bc not in sample_whitelist:
                continue
            if flag > 255:
                continue
            if chrom.startswith('chrM') or '_' in chrom or not chrom:
                continue

            # ---- type conversion ------------------------------------------------
            row = TagMapRow(
                read_id       = f[0],
                barcode       = bc,
                vt            = f[2],
                chrom         = chrom,
                start         = int(f[4]),
                end           = int(f[5]),
                flag          = flag,
                mapq          = int(f[7]),
                cigar         = f[8],
                block_starts  = _split_ints(f[9]),
                block_ends    = _split_ints(f[10]),
                sj_starts     = _split_ints(f[11]) if len(f) > 11 else [],
                sj_ends       = _split_ints(f[12]) if len(f) > 12 else [],
            )
            yield row



# ──────────────────────────────────────────────────────────────
#  write_feature_tsv – stream SAM/BAM ➜ 13-column TSV
# ──────────────────────────────────────────────────────────────
def write_feature_tsv(
    sam_path: str | pathlib.Path,
    tsv_path: str | pathlib.Path,
    *,
    min_mapq: int = 1,
    progress_every: int | None = 100_000,
) -> int:
    """
    Convert *sam_path* to the 13-column feature table defined by
    :func:`extract_alignment_features`, writing one line at a time.

    Parameters
    ----------
    sam_path
        Path to **coordinate-sorted** SAM or BAM (BAM requires `samtools view -h` pipe).
    tsv_path
        Output TSV (``.tsv`` or ``.tsv.gz``).
    min_mapq
        Minimum MAPQ to keep (see :func:`extract_alignment_features`).
    progress_every
        If set, log every *n* passing alignments to stderr.

    Returns
    -------
    int
        Number of rows written.
    """
    # prepare output handle (gzip if requested)
    tsv_path = pathlib.Path(tsv_path)
    opener = gzip.open if tsv_path.suffix == ".gz" else open
    rows_written = 0

    with SAMReader(sam_path) as rdr, opener(tsv_path, "wt") as out:
        out.write("\t".join(TSV_COLUMNS) + "\n")          # header once

        for aln in rdr:
            row = extract_alignment_features(aln, min_mapq=min_mapq)
            if row:
                out.write("\t".join(row) + "\n")
                rows_written += 1
                if progress_every and rows_written % progress_every == 0:
                    print(f"[write_feature_tsv] {rows_written:,} rows", file=sys.stderr)

    return rows_written

