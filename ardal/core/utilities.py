""" utilities.py
This module provides utility functions for working with GUIDs and alleles in the Ardal framework.
"""
import importlib
from collections import defaultdict

from ..exceptions.exceptions import *


def require_package(package_name, import_as=None, attr: str = None, error_message=None):
    """
    Attempt to import a package. Raise informative error if not found.

    Args:
        package_name (str): Name of the package to import (e.g. "matplotlib").
        import_as (str, optional): Name to import as (e.g. "plt" for matplotlib.pyplot).
        error_message (str, optional): Custom error message to display if not found.

    Returns:
        module: The imported module or submodule.

    Raises:
        ImportError: If the package is not installed.
    """
    try:
        module = importlib.import_module(import_as or package_name)
        if attr:
            return getattr(module, attr)
        return module
    except (ImportError, AttributeError):
        raise ImportError(
            error_message or
            f"The package '{package_name}' is required but not installed or missing an attribute. "
            f"Install it with `pip install {package_name}`."
        )



def checkGUIDs( guids : list,
                headers : dict,
                filter : bool = False ) -> list:
    """ Check guids are present within the matrix and construct a list of present guids to proceed with
    if filter is True
    """
    ## check input
    if not isinstance(guids, list):
        raise InvalidTypeError("guid_list must be a list.")
    
    absent_guids = []
    present_guids = []
    for id in guids:
        if id not in headers["guids"]:
            absent_guids.append(id)
        else:
            present_guids.append(id)

    if filter:
        return present_guids

    if len(absent_guids) > 0:
        raise InvalidGUIDQueryError(f"guids {absent_guids} not found in allele matrix.")


def checkAlleles( alleles : list,
                  headers : dict,
                  filter : bool = False ) -> list:
    """ Check alleles are present within the matrix and construct a list of present alleles to proceed with
    if filter is True
    """
    ## check input
    if not isinstance(alleles, list):
        raise InvalidTypeError("allele_list must be a list.")
    
    absent_alleles = []
    present_alleles = []
    for allele_id in alleles:
        if allele_id not in headers["alleles"]:
            absent_alleles.append(allele_id)
        else:
            present_alleles.append(allele_id)

    if filter:
        return present_alleles
    
    if len(absent_alleles) > 0:
        raise InvalidAlleleQueryError(f"alleles {absent_alleles} not found in allele matrix.")


def encodeGuid( guid : str,
                headers : dict ) -> int:
    """ Encode a GUID to its corresponding index in the headers dictionary.
    """
    return headers["guids"].index(guid)


def decodeGuid( row_coord : int,
                headers : dict ) -> str:
    """ Decode a row coordinate to its corresponding GUID in the headers dictionary.
    """
    return headers["guids"][row_coord]


def encodeAllele( allele : str,
                 headers : dict ) -> int:
    """ Encode an allele to its corresponding index in the headers dictionary.
    """
    return headers["alleles"].index(allele)


def decodeAllele( col_coord : int,
                 headers : dict ) -> str:
    """ Decode a column coordinate to its corresponding allele in the headers dictionary.
    """
    return headers["alleles"][col_coord]


def decodeAlleleID( allele_id : str,
                    allele_id_format : str = "{chr}.{start}.{ref}.{alt}" ) -> tuple :
    """ Decodes an allele ID string into its constituent parts based on a format string.

    Args:
        allele_id (str): The allele ID string to decode (e.g., `chr1.100.200.A.T``).
        allele_id_format (str): The format string defining the structure of the allele ID,
                                 using placeholders like {ref}, {chr}, {start}, {end}, {alt}.
                                 (e.g., "{chr}.{start}.{end}.{ref}.{alt}").

    Returns:
        tuple: A tuple containing the decoded positional tuple (chr, start, end, ref, alt).

    Raises:
        ValueError: If the allele_id does not match the allele_id_format.
    """
    re = require_package("re", "re")
    
    accepted_placeholder_patterns = ["chr", "start", "end", "ref", "alt"]
    
    if not allele_id_format:
        raise ValueError("allele_id_format cannot be empty.")
    if not allele_id:
        raise ValueError("allele_id cannot be empty.")
    
    ## check placeholders are valid
    pattern_placeholders = re.findall(r'\{([^}]+)\}', allele_id_format)
    for p in pattern_placeholders:
        if p not in accepted_placeholder_patterns:
            raise AllelePatternError(f"{p} not a valid placeholder for allele_id_format pattern. Accepted placeholders : {accepted_placeholder_patterns}")
        
    pattern = re.escape(allele_id_format)
    
    if not pattern:
        raise ValueError("allele_id_format could not be parsed. Please check your format string and ensure it contains valid placeholders.")

    ## define placeholders and their regex pattern
    placeholders = {
        'ref': r'(?P<ref>.+)',
        'chr': r'(?P<chr>.+)',
        'start': r'(?P<start>\d+)',
        'end': r'(?P<end>\d+)',
        'alt': r'(?P<alt>.+)'
    }

    found_placeholders = 0
    ## eeplace escaped placeholders with their regex capture groups
    for key, regex_pattern in placeholders.items():
        escaped_placeholder = re.escape(f"{{{key}}}")
        if escaped_placeholder in pattern:
            pattern = pattern.replace(escaped_placeholder, regex_pattern)
            found_placeholders += 1

    if found_placeholders == 0:
        raise ValueError("allele_id_format does not contain any valid placeholders (e.g., {chr}, {start}).")

    ## anchor the pattern to match the entire string
    pattern = f"^{pattern}$"

    ## compile the regex and match against the allele ID
    match = re.compile(pattern).match(allele_id)

    if not match:
        raise ValueError(f"Allele ID '{allele_id}' does not match the format '{allele_id_format}'.")

    parts = match.groupdict()

    start = int(parts['start']) if 'start' in parts else None
    end = int(parts['end']) if 'end' in parts else None

    return (
        parts.get('chr'),
        start,
        end,
        parts.get('ref'),
        parts.get('alt')
    )


def getIntervalAlleles( intervals : list,
                        headers : dict,
                        allele_id_format : str,
                        coords_bed : str=None ) -> set:
    bisect = require_package("bisect", "bisect")
    
    ## get allele positions
    allele_positions = getAllelePositions( headers,
                                           allele_id_format,
                                           coords_bed )
    ## find alleles which fall within the given intervals
    interval_alleles = []
    for interval in intervals:
        try:
            if len(interval) == 3:
                chr_key = str(interval[0])
                start = int(interval[1])
                end = int(interval[2])
            elif len(interval) == 2:
                chr_key = None
                start = int(interval[0])
                end = int(interval[1])
            elif len(interval) == 1:
                chr_key = str(interval[0])
                start = None
                end = None
            else:
                raise IntervalError(f"Intervals must be either [chr_key, start, end], [start, end] or [chr_key] : {interval}")
        except (ValueError, TypeError, IndexError) as e:
            raise IntervalError(f"Failed to parse interval {interval} : {e}")

        if chr_key not in allele_positions:
            ardal_warn(f"Chromosome {chr_key} not found in allele positions. Skipping...")
            continue

        chr_pos_dict = allele_positions[chr_key]
        if not chr_pos_dict:
            ardal_warn(f"No alleles found for chromosome {chr_key}. Skipping...")
            continue
        
        ## if only chromosome was passed as an interval, return all alleles in that chromosome
        if start is None and end is None:
            chr_alleles = [item for sublist in list(chr_pos_dict.values()) for item in sublist]
            interval_alleles.extend(chr_alleles)
            continue
        
        sorted_positions = sorted(chr_pos_dict.keys())
        
        ## binary search the alleles out
        start_index = bisect.bisect_left(sorted_positions, start)
        for i in range(start_index, len(sorted_positions)):
            pos = sorted_positions[i]
            if pos > end:
                break  ## exceeded the end of the interval
            interval_alleles.extend(chr_pos_dict[pos])
                
    if len(interval_alleles) == 0:
        ardal_warn(f"No alleles found within the given intervals. Please check your intervals and that allele IDs are in the expected form defined by allele_id_format (currently {allele_id_format}).")
        
    return set(interval_alleles)
            

def getAllelePositions( headers : dict,
                        allele_id_format : str,
                        coords_bed : str=None ) -> dict:

    if coords_bed:
        coords = readCoordsBED(coords_bed)
    else:
        coords = None

    ## get the position of each allele
    allele_positions = defaultdict(lambda: defaultdict(list))
    for allele in headers["alleles"]:
        try:
            if coords:
                if allele not in coords:
                    continue
                chr_key, start, end = coords[allele]
                allele_positions[str(chr_key)][int(start)].append(allele)
            else:
                chr_key, start, end, ref, alt = decodeAlleleID(allele, allele_id_format)
                if start is not None:
                    allele_positions[chr_key][start].append(allele)
        except ValueError as e:
            ardal_warn(f"Could not parse allele ID '{allele}': {e}")
            continue
        
    return allele_positions
        
        
def readCoordsBED( coords_bed : str ) -> dict:
    """ Read the coords BED file containing position : allele_id mappings.
    Must be structures as:
    chr start end allele_id
    """
    with open(coords_bed, 'r') as fin:
        lines = fin.readlines()
    
    coords = defaultdict(list)
    for line in lines:
        chr, start, end, allele_id = line.strip().split()
        coords[allele_id] = [chr, start, end]
    
    return coords
