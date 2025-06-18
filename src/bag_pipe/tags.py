from __future__ import annotations
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple


# ---------------------------------------------------------------------- #
# helpers
# ---------------------------------------------------------------------- #
def hamming(a: str, b: str) -> int:
    """Hamming distance for equal-length strings."""
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(aa != bb for aa, bb in zip(a, b))


def split_tag(tag: str, k: int) -> List[str]:
    """
    Split the tag into (k+1) non-overlapping substrings.
    Any two tags with Hamming ≤ k are guaranteed to share at least
    one identical block.
    """
    n_blocks = k + 1
    base_len, extra = divmod(len(tag), n_blocks)
    parts, start = [], 0
    for i in range(n_blocks):
        end = start + base_len + (1 if i < extra else 0)
        parts.append(tag[start:end])
        start = end
    return parts


# ---------------------------------------------------------------------- #
# main routine
# ---------------------------------------------------------------------- #
def rollup_tags(
    tags: List[Tuple[str, int]],
    k: int = 1,
    threshold: float = 0.10,
) -> Tuple[Set[str], Dict[str, str], List[Tuple[str, str]]]:
    """
    Collapse low-abundance tags onto higher-abundance parents.

    Parameters
    ----------
    tags       : list of (tag, read_count)
    k          : max Hamming distance
    threshold  : child_count < threshold * parent_count  (default 0.10)

    Returns
    -------
    keepers    : set[str]                 – tags that survive roll-up
    mapping    : dict[child] -> parent    – only for rolled tags
    rejects    : list[(child, parent)]    – same info, order preserved
    """
    # -- 1. sort by descending abundance ----------------------------------
    tags_sorted = sorted(tags, key=lambda t: t[1], reverse=True)
    counts: Dict[str, int] = {t: c for t, c in tags_sorted}
    tag_to_rank: Dict[str, int] = {t: i for i, (t, _) in enumerate(tags_sorted)}

    # -- 2. substring index & Union-Find ----------------------------------
    substring_index: Dict[str, Set[str]] = defaultdict(set)
    parent: Dict[str, str] = {}

    substr_cache: Dict[str, List[str]] = {}
    for tag, _ in tags_sorted:
        parts = split_tag(tag, k)
        substr_cache[tag] = parts
        for p in parts:
            substring_index[p].add(tag)
        parent[tag] = tag                      # each starts as its own root

    def find(tag: str) -> str:
        while parent[tag] != tag:
            parent[tag] = parent[parent[tag]]  # path compression
            tag = parent[tag]
        return tag

    # -- 3. main roll-up pass --------------------------------------------
    rejects: List[Tuple[str, str]] = []

    for tag, _ in tags_sorted:                # descending abundance
        root = find(tag)
        if root != tag:                       # already merged
            continue

        for sub in substr_cache[tag]:
            for cand in substring_index[sub]:
                if cand == tag:
                    continue
                cand_root = find(cand)
                # keep only higher-abundance candidates
                if tag_to_rank[cand_root] >= tag_to_rank[tag]:
                    continue
                # distance check
                if hamming(tag, cand_root) > k:
                    continue
                # abundance cutoff
                if counts[tag] >= threshold * counts[cand_root]:
                    continue

                # ---- roll up ------------------------------------------
                parent[tag] = cand_root
                rejects.append((tag, cand_root))
                root = cand_root             # not needed further, but clear
                break
            if parent[tag] != tag:
                break

    # -- 4. final keeper / mapping sets -----------------------------------
    keepers: Set[str] = {t for t in parent if parent[t] == t}
    mapping: Dict[str, str] = {t: p for t, p in parent.items() if t != p}

    return keepers, mapping, rejects
