"""
Allele identifier parsing, rendering, and ordering helpers.
"""

import re
from typing import Dict, List, Tuple

from .schemas import DEFAULT_ALLELE_ID_FORMAT

_ALLELE_ID_FORMAT: str | None = None
_ALLELE_ID_PATTERN: re.Pattern | None = None
_ALLELE_ID_POS_KEY: str | None = None
_ALLELE_ID_HAS_CHR: bool = True
_ALLELE_ID_HAS_GENE: bool = False

def allele_id_from_fields(
    *,
    chrom,
    start,
    end,
    ref,
    alt,
    gene = None,
    allele_id_format : str | None = None,
) -> str:
    """Render an allele record using the requested matrix-time identifier format."""
    fmt = allele_id_format or DEFAULT_ALLELE_ID_FORMAT
    pos = start
    values = {
        "chr": chrom,
        "start": start,
        "end": end,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "gene": gene,
    }
    try:
        return fmt.format(**values)
    except KeyError as exc:
        raise ValueError(f"Unsupported placeholder in allele_id_format: {exc}") from exc

def parse_allele_key(
    allele_id : str,
):
    """Split an allele identifier into chrom, position, reference, and alternate entries."""
    if ( _ALLELE_ID_PATTERN is None ):
        try:
            chrom, pos, ref, alt = allele_id.rsplit(".", 3)
        except ValueError as exc:
            raise ValueError(f"Unexpected allele identifier format: {allele_id}") from exc
        
        return chrom, pos, ref, alt

    match = _ALLELE_ID_PATTERN.match(allele_id)
    if ( not match ):
        raise ValueError(
            f"Allele ID '{allele_id}' does not match the format '{_ALLELE_ID_FORMAT}'."
        )

    parts = match.groupdict()
    chrom = parts.get("chr") if _ALLELE_ID_HAS_CHR else None
    pos = parts.get(_ALLELE_ID_POS_KEY)
    gene = parts.get("gene") if _ALLELE_ID_HAS_GENE else None
    ref = parts.get("ref")
    alt = parts.get("alt")

    if ( pos is None or ref is None or alt is None or ( _ALLELE_ID_HAS_CHR and chrom is None ) ):
        raise ValueError(
            "allele_id_format must include {pos}/{start}, {ref}, and {alt} placeholders."
        )

    if ( gene is not None ):
        pos = f"{gene}.{pos}"

    return chrom, pos, ref, alt

def compile_allele_id_format(
    allele_id_format : str,
) -> Tuple[re.Pattern, str, bool, bool]:
    """Compile an allele_id_format string into a regex pattern and return the position placeholder key."""
    if ( not allele_id_format ):
        raise ValueError("allele_id_format cannot be empty.")

    placeholders = re.findall(r"\{([^}]+)\}", allele_id_format)
    if ( not placeholders ):
        raise ValueError(
            "allele_id_format must include placeholders like {pos}, {ref}, {alt}."
        )

    allowed = {"chr", "gene", "pos", "start", "end", "ref", "alt"}
    invalid = sorted({p for p in placeholders if p not in allowed})
    if ( invalid ):
        raise ValueError(
            "Unsupported placeholders in allele_id_format: "
            + ", ".join(invalid)
            + f". Allowed: {', '.join(sorted(allowed))}."
        )

    seen = set()
    dupes = set()
    for p in placeholders:
        if p in seen:
            dupes.add(p)
        seen.add(p)
    if ( dupes ):
        raise ValueError(
            "allele_id_format contains duplicate placeholders: "
            + ", ".join(sorted(dupes))
        )

    if ( "pos" in placeholders and "start" in placeholders ):
        raise ValueError("allele_id_format must use either {pos} or {start}, not both.")

    pos_key = "pos" if "pos" in placeholders else "start"
    if ( pos_key not in placeholders ):
        raise ValueError("allele_id_format must include {pos} or {start}.")

    missing_required = []
    if ( "ref" not in placeholders ):
        missing_required.append("{ref}")
    if ( "alt" not in placeholders ):
        missing_required.append("{alt}")
    if ( missing_required ):
        raise ValueError(
            "allele_id_format missing required placeholders: "
            + ", ".join(missing_required)
        )

    pattern = re.escape(allele_id_format)
    replacements = {
        "ref": r"(?P<ref>.+)",
        "alt": r"(?P<alt>.+)",
        "chr": r"(?P<chr>.+)",
        "gene": r"(?P<gene>.+)",
        "start": r"(?P<start>\d+)",
        "pos": r"(?P<pos>\d+)",
        "end": r"(?P<end>\d+)",
    }
    for key, regex_pattern in replacements.items():
        escaped_placeholder = re.escape(f"{{{key}}}")
        if escaped_placeholder in pattern:
            pattern = pattern.replace(escaped_placeholder, regex_pattern)

    pattern = f"^{pattern}$"

    return re.compile(pattern), pos_key, ("chr" in placeholders), ("gene" in placeholders)

def configure_allele_id_format(
    allele_id_format : str | None,
):
    """Configure how allele identifiers are parsed within this module."""
    global _ALLELE_ID_FORMAT, _ALLELE_ID_PATTERN, _ALLELE_ID_POS_KEY, _ALLELE_ID_HAS_CHR, _ALLELE_ID_HAS_GENE

    if ( allele_id_format is None ):
        _ALLELE_ID_FORMAT = None
        _ALLELE_ID_PATTERN = None
        _ALLELE_ID_POS_KEY = None
        _ALLELE_ID_HAS_CHR = True
        _ALLELE_ID_HAS_GENE = False
        
        return

    pattern, pos_key, has_chr, has_gene = compile_allele_id_format(allele_id_format)
    _ALLELE_ID_FORMAT = allele_id_format
    _ALLELE_ID_PATTERN = pattern
    _ALLELE_ID_POS_KEY = pos_key
    _ALLELE_ID_HAS_CHR = has_chr
    _ALLELE_ID_HAS_GENE = has_gene

def allele_sort_key(
    allele_id : str,
) -> Tuple[str, int, object, str, str]:
    """Provide a deterministic ordering tuple for allele identifiers."""
    chrom, pos, ref, alt = parse_allele_key(allele_id)
    try:
        pos_val = int(pos)
        
        return chrom, 0, pos_val, ref, alt
    except ValueError:
        
        return chrom, 1, pos, ref, alt

def reorder_alleles_by_position(
    allele_to_idx : Dict[str, int],
) -> Tuple[Dict[str, int], List[str]]:
    """Provide a remapped allele index dict plus the ordered allele list."""
    ordered = sorted(allele_to_idx.keys(), key=allele_sort_key)
    remapped = {allele_id: idx for idx, allele_id in enumerate(ordered)}
    
    return remapped, ordered

def site_sort_key(
    site_id : str,
) -> Tuple[str, int, object, str]:
    """Return a deterministic ordering tuple for site identifiers."""
    parts = site_id.rsplit(".", 2)
    if ( len(parts) == 3 ):
        chrom, pos, ref = parts
    elif ( len(parts) == 2 ):
        chrom, pos, ref = "", parts[0], parts[1]
    else:
        raise ValueError(f"Unexpected site identifier format: {site_id}")
    try:
        pos_val = int(pos)
        
        return chrom, 0, pos_val, ref
    except ValueError:
        
        return chrom, 1, pos, ref

def site_identifier(
    chrom : str,
    pos : int | str,
    ref : str,
):
    """Construct the canonical identifier for a site."""
    if ( chrom is None or chrom == "" ):
        return f"{pos}.{ref}"
    return f"{chrom}.{pos}.{ref}"
