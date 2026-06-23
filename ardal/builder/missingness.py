"""
Missingness parsing and matrix-column projection helpers.
"""

import bisect
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from .vcf import open_maybe_gzip

def position_key(
    chrom : str,
    pos : str | int,
) -> str:
    """Compose a string key for positional lookups."""
    if ( chrom is None or chrom == "" ):
        return f"{pos}"
    
    return f"{chrom}.{pos}"

def normalize_missing_entries(
    entries : Iterable[object],
) -> List[str]:
    """Normalize missing site markers to strings.

    Supports string site keys as well as [chrom, pos] and (chrom, pos) tuples.
    """
    normalized = []
    for entry in entries or []:
        if isinstance(entry, (list, tuple)):
            if len(entry) == 2:
                chrom, pos = entry
                entry_str = position_key(str(chrom), pos)
            elif len(entry) == 1:
                entry_str = str(entry[0])
            else:
                entry_str = str(entry)
        else:
            entry_str = str(entry)

        if ( not _ALLELE_ID_HAS_CHR ):
            if ( "." in entry_str ):
                _, tail = entry_str.rsplit(".", 1)
                if ( tail.isdigit() ):
                    entry_str = tail
        normalized.append(entry_str)

    return normalized

def missing_blocks_from_entries(
    entries : Iterable[object],
) -> List[Tuple[str, int, int]]:
    """Parse v2 missing blocks as (chrom, start, end), preserving only valid numeric intervals."""
    blocks: List[Tuple[str, int, int]] = []
    for entry in entries or []:
        if not isinstance(entry, (list, tuple)) or len(entry) != 3:
            continue
        chrom, start, end = entry
        try:
            start_i = int(start)
            end_i = int(end)
        except (TypeError, ValueError):
            continue
        if end_i < start_i:
            start_i, end_i = end_i, start_i
        blocks.append((str(chrom), start_i, end_i))

    return blocks

def missing_entry_contains_position(entry, chrom: str, pos: str | int) -> bool:
    """Return whether a legacy missing key or v2 missing block contains a site."""
    try:
        pos_i = int(pos)
    except (TypeError, ValueError):
        pos_i = None

    if isinstance(entry, (list, tuple)) and len(entry) == 3:
        block_chrom, start, end = entry
        if str(block_chrom) != str(chrom) or pos_i is None:
            return False
        try:
            return int(start) <= pos_i <= int(end)
        except (TypeError, ValueError):
            return False

    return str(entry) == position_key(str(chrom), pos)

def any_missing_entry_contains_position(entries, chrom: str, pos: str | int) -> bool:
    return any(missing_entry_contains_position(entry, chrom, pos) for entry in entries or [])

def compact_index_ranges(indices: Iterable[int]) -> List[List[int]]:
    """Compact sorted integer indices into inclusive [start, end] ranges."""
    ranges: List[List[int]] = []
    sorted_indices = sorted(set(int(i) for i in indices))
    if not sorted_indices:
        return ranges
    start = prev = sorted_indices[0]
    for idx in sorted_indices[1:]:
        if idx == prev + 1:
            prev = idx
            continue
        ranges.append([start, prev])
        start = prev = idx
    ranges.append([start, prev])
    return ranges

def expand_index_ranges(ranges: Iterable[Sequence[int]]) -> set:
    expanded = set()
    for item in ranges or []:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            continue
        start, end = item
        expanded.update(range(int(start), int(end) + 1))
    return expanded

def build_columns_by_genomic_site(
    ordered_alleles : Sequence[str],
    allele_to_idx : Dict[str, int],
    allele_genomic_sites : Dict[str, Tuple[str, int]],
) -> Dict[str, Tuple[List[int], List[List[int]]]]:
    """Index matrix columns by genomic coordinate for fast missing-block projection."""
    print("[Build missing-site column index]", flush=True)
    grouped: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(list))
    for allele_id in ordered_alleles:
        genomic_site = allele_genomic_sites.get(allele_id)
        if genomic_site is None:
            continue
        chrom, pos = genomic_site
        grouped[str(chrom)][int(pos)].append(allele_to_idx[allele_id])

    indexed: Dict[str, Tuple[List[int], List[List[int]]]] = {}
    for chrom, pos_to_cols in grouped.items():
        positions = sorted(pos_to_cols)
        columns = [pos_to_cols[pos] for pos in positions]
        indexed[chrom] = (positions, columns)
    print(f"[Build missing-site column index] {len(indexed)} chromosomes", flush=True)

    return indexed

def column_masks_from_missing_blocks(
    sample_to_idx,
    sample_missing_positions : Dict[str, list],
    columns_by_site : Dict[str, Tuple[List[int], List[List[int]]]],
) -> Dict[str, List[List[int]]]:
    """Project missing intervals onto existing matrix allele columns."""
    column_masks = {}
    total = len(sample_to_idx)
    for i, sample_id in enumerate(sample_to_idx.keys(), 1):
        if ( i == 1 or i % 25 == 0 or i == total ):
            print(f"\r[Missing masks] {i}/{total}: {sample_id}\x1b[K", end="", flush=True)
        masked_cols = []
        for entry in sample_missing_positions.get(sample_id, []):
            if isinstance(entry, (list, tuple)) and len(entry) == 3:
                chrom, start, end = entry
                site_index = columns_by_site.get(str(chrom))
                if site_index is None:
                    continue
                positions, columns = site_index
                left = bisect.bisect_left(positions, int(start))
                right = bisect.bisect_right(positions, int(end))
                for cols_at_pos in columns[left:right]:
                    masked_cols.extend(cols_at_pos)
            else:
                if ( "." not in str(entry) ):
                    continue
                chrom, pos = str(entry).rsplit(".", 1)
                site_index = columns_by_site.get(chrom)
                if site_index is None:
                    continue
                positions, columns = site_index
                try:
                    pos_i = int(pos)
                except ValueError:
                    continue
                idx = bisect.bisect_left(positions, pos_i)
                if idx < len(positions) and positions[idx] == pos_i:
                    masked_cols.extend(columns[idx])
        compacted = compact_index_ranges(masked_cols)
        if ( compacted ):
            column_masks[sample_id] = compacted
    print(flush=True)

    return column_masks

def append_column_range(
    ranges : List[List[int]],
    start : int,
    end : int,
):
    """Append a [start, end) column range, merging with the previous range when adjacent."""
    if ( end <= start ):
        return
    if ( ranges and start <= ranges[-1][1] ):
        if ( end > ranges[-1][1] ):
            ranges[-1][1] = end
        return
    ranges.append([int(start), int(end)])

def missing_column_ranges_from_entries(
    entries : Iterable[object],
    columns_by_site : Dict[str, Tuple[List[int], List[List[int]]]],
) -> List[List[int]]:
    """Project one sample's missing genomic entries to compact matrix-column [start, end) ranges."""
    ranges: List[List[int]] = []
    for entry in entries or []:
        if isinstance(entry, (list, tuple)) and len(entry) == 3:
            chrom, start, end = entry
            site_index = columns_by_site.get(str(chrom))
            if site_index is None:
                continue
            positions, columns = site_index
            left = bisect.bisect_left(positions, int(start))
            right = bisect.bisect_right(positions, int(end))
            for cols_at_pos in columns[left:right]:
                if ( not cols_at_pos ):
                    continue
                append_column_range(ranges, min(cols_at_pos), max(cols_at_pos) + 1)
            continue

        if ( "." not in str(entry) ):
            continue
        chrom, pos = str(entry).rsplit(".", 1)
        site_index = columns_by_site.get(chrom)
        if site_index is None:
            continue
        positions, columns = site_index
        try:
            pos_i = int(pos)
        except ValueError:
            continue
        idx = bisect.bisect_left(positions, pos_i)
        if idx < len(positions) and positions[idx] == pos_i:
            cols_at_pos = columns[idx]
            if ( cols_at_pos ):
                append_column_range(ranges, min(cols_at_pos), max(cols_at_pos) + 1)

    return ranges

def missing_blocks_from_bed(depth_path: Path, min_depth: int) -> List[List]:
    missing: List[List] = []
    current_chrom: str | None = None
    current_start: int | None = None
    current_end: int | None = None

    def flush_current() -> None:
        nonlocal current_chrom, current_start, current_end
        if current_chrom is not None and current_start is not None and current_end is not None:
            missing.append([current_chrom, current_start, current_end])
        current_chrom = None
        current_start = None
        current_end = None

    with open_maybe_gzip(depth_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, pos, depth_str = parts[:3]
            try:
                depth_val = int(float(depth_str))
            except ValueError:
                continue
            if depth_val >= min_depth:
                flush_current()
                continue
            try:
                pos_int = int(pos)
            except ValueError:
                continue
            if current_chrom == chrom and current_end is not None and pos_int == current_end + 1:
                current_end = pos_int
                continue
            flush_current()
            current_chrom = chrom
            current_start = pos_int
            current_end = pos_int
    flush_current()
    return missing
