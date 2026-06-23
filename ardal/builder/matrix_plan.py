"""
Matrix planning and site-statistics helpers.
"""

import sys
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Dict, Sequence, Tuple

from .alleles import parse_allele_key, site_identifier, site_sort_key
from .sample_json import iter_sample_projected_entries

def allele_genomic_sites_from_json(
    json_files : Sequence[Path],
    data_type : str,
    allele_id_format : str | None,
    nonsyn_allele_id_format : str | None,
) -> Dict[str, Tuple[str, int]]:
    """Map matrix allele IDs to the genomic positions used for missing-block masks."""
    allele_sites: Dict[str, Tuple[str, int]] = {}
    total = len(json_files)
    for i, jf in enumerate(json_files, 1):
        if ( i % 25 == 0 or i == total ):
            print(f"\r[Allele genomic sites] {i}/{total}: {jf}\x1b[K", end="", flush=True)
        for _sample_id, allele_entries, _missing in iter_sample_projected_entries(
            jf,
            data_type=data_type,
            allele_id_format=allele_id_format,
            nonsyn_allele_id_format=nonsyn_allele_id_format,
        ):
            for allele_id, chrom, pos in allele_entries:
                if allele_id not in allele_sites:
                    try:
                        allele_sites[allele_id] = (str(chrom), int(pos))
                    except (TypeError, ValueError):
                        pass
    print(flush=True)

    return allele_sites

def index_samples_and_alleles(
    json_files : Sequence[Path],
    data_type : str = "snps",
    allele_id_format : str | None = None,
    nonsyn_allele_id_format : str | None = None,
    collect_site_alt_samples : bool = False,
):
    """
    Description:
        Traverse all sample JSON artifacts to build lookup tables for samples, alleles, and site-level stats.

    Inputs:
        json_files: Sequence of JSON paths emitted by upstream allele callers.

    Outputs:
        Returns mapping dictionaries plus site/sample metadata required for downstream matrix construction.

    Exceptions:
        Propagates parsing errors raised while decoding malformed JSON.
    """
    sample_to_idx, allele_to_idx = {}, {}
    site_alt_counts = defaultdict(lambda: defaultdict(int))
    site_alt_samples = defaultdict(set) if collect_site_alt_samples else None
    site_alt_sample_counts = defaultdict(int)
    site_reference = {}
    site_coordinates = {}
    sample_missing_positions: Dict[str, list] = defaultdict(list)
    
    for ( i, jf ) in enumerate(json_files, 1):
        print(f"\r[Indexing] {i}/{len(json_files)}: {jf}\x1b[K", end="", flush=True)
        for sample_id, allele_entries, missing in iter_sample_projected_entries(
            jf,
            data_type=data_type,
            allele_id_format=allele_id_format,
            nonsyn_allele_id_format=nonsyn_allele_id_format,
        ):
            if ( sample_id not in sample_to_idx ):
                sample_to_idx[sample_id] = len(sample_to_idx)
            sample_idx = sample_to_idx[sample_id]
            if ( missing ):
                sample_missing_positions[sample_id].extend(missing)
            if ( not allele_entries ):
                continue
            seen = set()
            per_site_in_sample = set()
            for allele_id, _genomic_chrom, _genomic_pos in allele_entries:
                if ( allele_id in seen ):
                    continue
                seen.add(allele_id)
                chrom, pos, ref, _ = parse_allele_key(allele_id)
                site_id = site_identifier(chrom, pos, ref)
                site_reference[site_id] = ref
                site_coordinates[site_id] = {"chrom": chrom, "pos": pos}
                if ( allele_id not in allele_to_idx ):
                    allele_to_idx[allele_id] = len(allele_to_idx)
                site_alt_counts[site_id][allele_id] += 1
                per_site_in_sample.add(site_id)
            for site_id in per_site_in_sample:
                if ( site_alt_samples is not None ):
                    site_alt_samples[site_id].add(sample_idx)
                site_alt_sample_counts[site_id] += 1
    print()
    
    return (
        sample_to_idx,
        allele_to_idx,
        site_reference,
        site_coordinates,
        site_alt_counts,
        site_alt_samples if site_alt_samples is not None else site_alt_sample_counts,
        sample_missing_positions,
    )

def derive_site_statistics(
    sample_to_idx,
    allele_to_idx,
    site_reference,
    site_coordinates,
    site_alt_counts,
    site_alt_samples,
):
    """
    Description:
        Summarize the major/minor allele makeup for every site observed in the dataset.

    Inputs:
        sample_to_idx / allele_to_idx: Index lookups constructed earlier.
        site_reference / site_coordinates: Context for reference alleles and positions.
        site_alt_counts / site_alt_samples: Aggregated counts for alt alleles and carriers.

    Outputs:
        Tuple of (site_metadata dict, minor column labels, alt-to-idx map, ref-to-idx map).

    Exceptions:
        Propagates parsing errors arising from malformed allele identifiers.
    """
    n_samples = len(sample_to_idx)
    site_metadata = OrderedDict()
    minor_columns = []
    minor_alt_to_idx = {}
    minor_ref_to_idx = {}

    ## iterate sites in sorted order for deterministic output
    for site_id in sorted(site_reference.keys(), key=site_sort_key):
        ref_base = site_reference[site_id]
        coords = site_coordinates[site_id]
        alt_counts_map = site_alt_counts.get(site_id, {})
        ## convert to regular dict for stable iteration later
        alt_counts = dict(alt_counts_map)
        alt_sample_entry = site_alt_samples.get(site_id, set())
        alt_sample_count = (
            len(alt_sample_entry)
            if isinstance(alt_sample_entry, set)
            else int(alt_sample_entry)
        )
        ref_count = n_samples - alt_sample_count

        major_label = "REF"
        major_count = ref_count

        for allele_id, count in alt_counts.items():
            if ( count > major_count ):
                major_label = allele_id
                major_count = count
            elif ( count == major_count ):
                if ( major_label == "REF" ):
                    continue
                if ( major_label != "REF" and allele_id < major_label ):
                    major_label = allele_id

        if ( major_label == "REF" ):
            major_type = "REF"
            major_base = ref_base
            major_allele_id = None
        else:
            major_type = "ALT"
            _, _, _, alt_base = parse_allele_key(major_label)
            major_base = alt_base
            major_allele_id = major_label

        minor_alt_labels = []
        alt_items = sorted(alt_counts.items(), key=lambda kv: allele_to_idx.get(kv[0], sys.maxsize))
        for allele_id, _count in alt_items:
            if ( allele_id == major_label ):
                continue
            idx = len(minor_columns)
            minor_columns.append(allele_id)
            minor_alt_to_idx[allele_id] = idx
            minor_alt_labels.append(allele_id)

        ref_minor_label = None
        if ( major_label != "REF" ):
            ref_minor_label = f"{site_id}.REF"
            idx = len(minor_columns)
            minor_columns.append(ref_minor_label)
            minor_ref_to_idx[site_id] = idx

        alt_allele_details = []
        for allele_id, count in alt_items:
            _, _, _, alt_base = parse_allele_key(allele_id)
            alt_allele_details.append(
                {
                    "allele_id": allele_id,
                    "allele": alt_base,
                    "count": count,
                    "is_minor": allele_id in minor_alt_to_idx,
                }
            )

        site_metadata[site_id] = {
            "chrom": coords["chrom"],
            "pos": coords["pos"],
            "ref": ref_base,
            "total_samples": n_samples,
            "alt_sample_count": alt_sample_count,
            "major": {
                "type": major_type,
                "allele": major_base,
                "allele_id": major_allele_id,
                "count": major_count,
            },
            "alt_alleles": alt_allele_details,
            "minor_columns": {
                "ref": ref_minor_label,
                "alt": minor_alt_labels,
            },
        }
    
    return site_metadata, minor_columns, minor_alt_to_idx, minor_ref_to_idx
