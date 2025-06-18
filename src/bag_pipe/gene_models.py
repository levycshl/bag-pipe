# bag_pipe/gene_models.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Iterable
import pandas as pd

from bisect import bisect_left

from bag_pipe.io import TagMapRow

from typing import List, Set, Tuple



Interval = Tuple[int, int]                 # (start, end) inclusive
FeatureMap = Dict[int, List[Interval]]     # gene_id → intervals

@dataclass(slots=True)
class GeneModel:
    """
    GeneModel encapsulates **one chromosome** worth of refFlat annotation.

    Parameters
    ----------
    ref_df : DataFrame
        Output of `load_refflat(path, chrom=...)` (1-chrom, 1-based, list cols).
    """

    ref_df: pd.DataFrame
    _chrom: str = field(init=False)

    # private caches  (populated lazily)
    _exons: FeatureMap            = field(default_factory=dict, init=False)
    _cds: FeatureMap              = field(default_factory=dict, init=False)
    _utr5: FeatureMap             = field(default_factory=dict, init=False)
    _utr3: FeatureMap             = field(default_factory=dict, init=False)
    _splice_sites: FeatureMap     = field(default_factory=dict, init=False)
    _splice_junctions: Dict[int, Dict[str, int]] = field(default_factory=dict, init=False)
    # for transcript bounds and quick overlap checks
    _tx_bounds: List[Tuple[int, int, int]] = field(default_factory=list, init=False)
    _starts: List[int] = field(default_factory=list, init=False)  # convenience for bisect

    def __post_init__(self):
        # inside GeneModel.__post_init__
        self._tx_bounds = []  # list[(start, end, gid)]
        for gid, row in enumerate(self.ref_df.itertuples(index=False), 1):
            self._tx_bounds.append((row.tx_start, row.tx_end, gid))
        # sort by tx_start once
        self._tx_bounds.sort(key=lambda x: x[0])
        # convenience list of starts for bisect
        self._starts = [s for s, _, _ in self._tx_bounds]
        if self.ref_df["chrom"].nunique() != 1:
            raise ValueError("GeneModel expects a single-chromosome DataFrame")
        self._chrom = self.ref_df["chrom"].iat[0]

    def overlapping_gids(self, q_start: int, q_end: int) -> List[int]:
        """
        Return gene IDs whose transcript span overlaps [q_start, q_end].
        """
        idx = bisect_left(self._starts, q_end + 1)  # first start > queryEnd
        candidates = []
        for s, e, gid in self._tx_bounds[:idx]:
            if e < q_start:  # ends before query → no overlap
                continue
            candidates.append(gid)
        return candidates

    # ------------------------------------------------------------------
    # Public – simple accessors
    # ------------------------------------------------------------------
    def exons(self, gid: int) -> List[Interval]:
        return self._build_if_missing(gid, "_exons")

    def cds_exons(self, gid: int) -> List[Interval]:
        return self._build_if_missing(gid, "_cds")

    def utr5(self, gid: int) -> List[Interval]:
        return self._build_if_missing(gid, "_utr5")

    def utr3(self, gid: int) -> List[Interval]:
        return self._build_if_missing(gid, "_utr3")

    def splice_sites(self, gid: int) -> List[Interval]:
        return self._build_if_missing(gid, "_splice_sites")

    def splice_junctions(self, gid: int) -> Dict[str, int]:
        return self._build_if_missing(gid, "_splice_junctions")

    # ------------------------------------------------------------------
    # Public – composite helper used by roll-up scoring
    # ------------------------------------------------------------------
    def metrics(
        self,
        gid: int,
        union: List[Interval],
        vt_junc_keys: Iterable[str],
    ) -> Dict[str, int]:
        """Calculate the exact seven numbers used in the tie-breaker cascade."""
        return {
            "basesTx"  : _covered_bp(union, [(self.tx_start(gid), self.tx_end(gid))]),
            "basesEx"  : _covered_bp(union, self.exons(gid)),
            "nSJ"      : sum(1 for k in vt_junc_keys if k in self.splice_junctions(gid)),
            "basesCds" : _covered_bp(union, self.cds_exons(gid)),
            "bases5"   : _covered_bp(union, self.utr5(gid)),
            "bases3"   : _covered_bp(union, self.utr3(gid)),
            "minSS"    : _covered_splice_sites(self.splice_sites(gid), union),
        }

    # ------------------------------------------------------------------
    # Convenience wrappers for tx bounds
    # ------------------------------------------------------------------
    def tx_start(self, gid: int) -> int:
        return self.ref_df.iloc[gid - 1].tx_start     # enumerate() is 1-based
    def tx_end(self, gid: int) -> int:
        return self.ref_df.iloc[gid - 1].tx_end

    # ------------------------------------------------------------------
    # Internal – lazy builder
    # ------------------------------------------------------------------
    def _build_if_missing(self, gid: int, attr: str):
        store = getattr(self, attr)
        if gid not in store:
            self._populate(gid)
        return store[gid]

    def _populate(self, gid: int):
        """Fill every cached feature for the given gene_id in one shot."""
        row = self.ref_df.iloc[gid - 1]           # 0-based index in DataFrame
        xs, xe = row.exon_starts, row.exon_ends
        ex_iv  = list(zip(xs, xe))
        store  = lambda name, val: getattr(self, name).__setitem__(gid, val)

        # exons
        store("_exons", ex_iv)

        # CDS ∩ exon
        cds = [
            (max(s, row.cds_start), min(e, row.cds_end))
            for s, e in ex_iv if max(s, row.cds_start) <= min(e, row.cds_end)
        ]
        store("_cds", cds)

        # 5'/3' UTR
        if row.strand == "+":
            utr5_reg = (row.tx_start, max(row.cds_start - 1, row.tx_start))
            utr3_reg = (min(row.cds_end + 1, row.tx_end), row.tx_end)
        else:
            utr5_reg = (min(row.cds_end + 1, row.tx_end), row.tx_end)
            utr3_reg = (row.tx_start, max(row.cds_start - 1, row.tx_start))

        store("_utr5", _intersect_region(utr5_reg, ex_iv))
        store("_utr3", _intersect_region(utr3_reg, ex_iv))

        # splice sites (2-bp windows)
        ss: List[Interval] = []
        if row.exon_count > 1:
            for i in range(row.exon_count):
                if i == 0:
                    ss.append((xe[i], xe[i] + 1))
                elif i == row.exon_count - 1:
                    ss.append((xs[i] - 1, xs[i]))
                else:
                    ss.extend([(xs[i] - 1, xs[i]), (xe[i], xe[i] + 1)])
        store("_splice_sites", ss)

        # splice-junction dict
        sj: Dict[str, int] = {}
        for i in range(1, row.exon_count):
            sj[f"{self._chrom}.{xe[i-1]}.{xs[i]}"] = 0
        store("_splice_junctions", sj)


# ────────────────────────────────────────────────────────────────────
# Helper utilities (module-internal)
# ────────────────────────────────────────────────────────────────────
def _intersect_region(region: Interval, exons: List[Interval]) -> List[Interval]:
    s, e = region
    if s > e:
        return []
    return [(max(s, xs), min(e, xe)) for xs, xe in exons if max(s, xs) <= min(e, xe)]

def _covered_bp(query: List[Interval], target: List[Interval]) -> int:
    """Return total bases of **query** that overlap any interval in target."""
    if not query or not target:
        return 0
    tot = 0
    for qs, qe in query:
        for ts, te in target:
            s = max(qs, ts)
            e = min(qe, te)
            if s <= e:
                tot += e - s + 1
    return tot

def _covered_splice_sites(ss: List[Interval], query: List[Interval]) -> int:
    """Count splice-site windows touched by query."""
    return sum(
        1
        for s, e in ss
        for qs, qe in query
        if max(s, qs) <= min(e, qe)
    )


@dataclass(slots=True)
class VTSummary:
    barcode:       str
    vt:            str
    chrom:         str
    min_start:     int
    max_end:       int
    union_starts:  List[int]
    union_ends:    List[int]
    span_bp:       int
    read_count:    int
    gene_count:    int
    genes:         List[str]      # list of gene symbols, empty -> NOGENE
    bases_tx:      int
    bases_ex:      int
    n_sj:          int
    bases_cds:     int
    bases_5utr:    int
    bases_3utr:    int
    min_ss:        int

    def header(self):
        return ["barcode", "vt", "chrom","min_start", "max_end", "union_starts", "union_ends",
                "span_bp", "read_count", "gene_count", "genes", "bases_tx", "bases_ex",
                "n_sj", "bases_cds", "bases_5utr", "bases_3utr", "min_ss"]

    # ------------------------------------------------------------
    # serialization helpers
    # ------------------------------------------------------------
    def to_row(self) -> List[str]:
        """Return values in legacy order, ready for '\t'.join(row)."""
        return [
            self.barcode,
            self.vt,
            self.chrom,
            str(self.min_start), str(self.max_end),
            ",".join(map(str, self.union_starts)),
            ",".join(map(str, self.union_ends)),
            str(self.span_bp),
            str(self.read_count),
            str(self.gene_count),
            ",".join(self.genes) if self.genes else "NOGENE",
            str(self.bases_tx), str(self.bases_ex), str(self.n_sj),
            str(self.bases_cds), str(self.bases_5utr), str(self.bases_3utr),
            str(self.min_ss)
        ]

    def __str__(self) -> str:
        return "\t".join(self.to_row())


def union_intervals(rows: List[TagMapRow]) -> List[Interval]:
    """
    Merge all CIGAR blocks of the VT into a non-overlapping list[ (start,end) ].
    Assumes every row comes from the same chromosome.
    """
    intervals: List[Interval] = []
    for r in rows:
        intervals.extend(r.cigar_intervals)          # property on TagMapRow

    if not intervals:
        return []

    intervals.sort(key=lambda x: x[0])               # sort by start
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        if s > merged[-1][1] + 1:             # gap → new segment
            merged.append((s, e))
        else:                                 # overlap/adjacent → extend
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
    return merged


def splice_keys(rows: List[TagMapRow]) -> Set[str]:
    """
    Collect splice-junction keys ("chr.end.start") observed in the VT’s reads.
    """
    keys: Set[str] = set()
    for r in rows:
        keys.update(r.splice_keys())          # method on TagMapRow
    return keys


def composite_key(m: dict) -> tuple:
    """
    Create a tuple of metrics for lexicographic comparison.
    This tuple is used to determine the "best" transcript based on
    the tie-breaker cascade.
    The tuple is ordered such that a *larger* tuple is considered "better".
    The order of the tuple elements is important and follows the cascade:
    basesTx, basesEx, nSJ, basesCds, bases5, bases3, -minSS
    This means that:
    - more bases in the transcript is better
    - more bases in the exons is better
    - more splice junctions is better
    - more bases in the CDS is better
    - more bases in the 5' UTR is better
    - more bases in the 3' UTR is better
    - a *smaller* minSS is better (hence the negation)
    """
    # negate minSS so that a *smaller* minSS makes the tuple larger lexicographically
    return (m["basesTx"], m["basesEx"], m["nSJ"],
            m["basesCds"], m["bases5"], m["bases3"],
            -m["minSS"])


def summarise_vt(rows: list[TagMapRow], model: GeneModel) -> VTSummary:
    """
    Summarise a VT (variant transcript) from a list of TagMapRow objects.
    This function computes the union of all intervals in the VT,
    determines the best transcript based on the tie-breaker cascade,
    and returns a VTSummary object containing the relevant metrics.
    The summary includes:
    - barcode
    - vt (variant transcript identifier)
    - chromosome
    - minimum start position
    - maximum end position
    - union of intervals (start and end positions)
    - span in base pairs
    - read count
    - gene count
    - list of genes
    - bases in transcript
    - bases in exons
    - number of splice junctions
    - bases in CDS
    - bases in 5' UTR
    - bases in 3' UTR
    - minimum splice site score
    :param rows:
    :param model:
    """
    union = union_intervals(rows)
    sj_keys = splice_keys(rows)

    min_start = min(r.start for r in rows)
    max_end   = max(r.end   for r in rows)

    # --- pick best transcript --------------------------------------------
    def metrics_for_gid(gid: int):
        return model.metrics(gid, union, sj_keys)

    cands = model.overlapping_gids(min_start, max_end)
    metric_items = [(gid, metrics_for_gid(gid))
                    for gid in cands]

    best_key = None
    best_genes = set()
    for gid, m in metric_items:
        if (m["basesTx"] == 0 and m["basesEx"] == 0 and
        m["nSJ"] == 0 and m["basesCds"] == 0):
            continue
        key = composite_key(m)
        if best_key is None or key > best_key:
            best_key = key
            best_genes = {model.ref_df.iloc[gid-1].gene}
        elif key == best_key:
            best_genes.add(model.ref_df.iloc[gid-1].gene)

    if best_key is None:
        best_key = (0, 0, 0, 0, 0, 0, 0)  # no valid transcript found
    basesTx,basesEx,nSJ,basesCds,bases5,bases3,neg_minSS = best_key
    minSS = -neg_minSS
    genes_list = sorted(best_genes)
    gene_count = len(genes_list)
    span_bp    = sum(e-s+1 for s,e in union)

    r0 = rows[0]
    return VTSummary(
        barcode      = r0.barcode,
        vt           = r0.vt,
        chrom        = r0.chrom,
        min_start    = min_start,
        max_end      = max_end,
        union_starts = [s for s, _ in union],
        union_ends   = [e for _, e in union],
        span_bp      = span_bp,
        read_count   = len(rows),
        gene_count   = gene_count,
        genes        = genes_list,
        bases_tx     = basesTx,
        bases_ex     = basesEx,
        n_sj         = nSJ,
        bases_cds    = basesCds,
        bases_5utr   = bases5,
        bases_3utr   = bases3,
        min_ss       = minSS,
    )
