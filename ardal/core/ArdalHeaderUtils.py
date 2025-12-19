import os
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Union, List, Tuple, Set, Dict

from ..utils.misc import require_package
from ..utils.decorators import check_allele_id_format, check_bed_paths, check_alleles_list, check_guids_list, validate_type
from ..utils.exceptions import AllelePatternError, IntervalError, InvalidGUIDQueryError, InvalidAlleleQueryError, ParameterError
from ..utils.logger import get_logger

log = get_logger()


# core/ArdalHeaderUtils.py
class ArdalHeaderUtils:

    def __init__( self,
                  headers : dict,
                  meta : dict,
                  allele_coords_bed : Union[str, None] = None,
                  allele_id_format : Union[str, None] = None,
                  allele_positions : Union[Dict, None] = None,
                  missing_masks: Union[Dict, None] = None ):

        self.headers = headers
        self.meta = meta
        self.n_guids = len(self.headers['guids'])
        self.n_alleles = len(self.headers['alleles'])
        self.allele_positions = None
        self._allele_coords_bed = allele_coords_bed
        self._allele_id_format = allele_id_format
        self._decoded_headers = self._decode_headers()
        self._snv_lookup_cache = None
        self._snv_lookup_format = None
        self._site_lookup_cache = None
        self._missing_masks: Dict[str, List] = {}
        self._missing_mask_rows: Union[List[List[int]], None] = None
        self._allele_position_map_cache: Union[Dict, None] = None
        self._site_count_all_cache: Union[int, None] = None
        
        ## allele positions were not provided
        if not allele_positions:
            self.allele_positions = defaultdict(lambda: defaultdict(list))
            
            ## position cache switches
            self._allele_positions_from_bed = False   ## positions provided from bed
            self._allele_positions_from_ids = False   ## positions decoded from allele_ids
            
            ## if allele_id_format is provided then just decode positions now
            if self._allele_id_format :
                self.allele_positions = self.compute_allele_positions(allele_id_format = self._allele_id_format,
                                                                      allele_coords_bed = self._allele_coords_bed)
        
        ## allele positions were provided in a dict
        else:
            ## TODO: some checks on this
            self.allele_positions = allele_positions
            self._allele_positions_from_bed = False
            self._allele_positions_from_ids = True

        self._init_missing_masks(missing_masks)
            
        log.info(f"""Initialised HeaderUtils object with {len(self.headers['guids'])} GUIDs and {len(self.headers['alleles'])} alleles.
                                    allele_coords_bed = {self._allele_coords_bed}
                                    allele_id_format = {self._allele_id_format}
                                    allele_positions_from_bed = {self._allele_positions_from_bed}
                                    allele_positions_from_ids = {self._allele_positions_from_ids}""")


    def _decode_headers( self ) -> Dict:
        log.info(f"Decoding headers.")
        return  { "guids" : dict(zip(self.headers["guids"], range(len(self.headers["guids"])))),
                  "alleles" : dict(zip(self.headers["alleles"], range(len(self.headers["alleles"]))))
                }


    def get_snv_lookup( self,
                        allele_id_format : Union[str, None] = None,
                        *,
                        force: bool = False
                        ) -> Tuple[np.ndarray, np.ndarray, int]:
        """
        Return lookup arrays mapping alleles to SNV loci and nucleotide bins.

        Args:
            allele_id_format: Format string used to decode allele identifiers.
            force: When True, rebuild the lookup even if the cached format matches.
        """
        ## return cache if it already exists and not forced to rebuild
        if (
            not force
            and self._snv_lookup_cache is not None
            and self._snv_lookup_format == allele_id_format
        ):
            return self._snv_lookup_cache

        ## else build/rebuild
        allele_id_format = allele_id_format or self._allele_id_format
        if allele_id_format is None:
            raise AllelePatternError("Allele ID format required to build SNV lookup.")
        
        pattern = self._check_allele_format_grammar(allele_id_format=allele_id_format)
        base_map = {"A": 1, "C": 2, "G": 3, "T": 4}

        n_alleles = self.n_alleles
        allele_to_locus = np.full(n_alleles, np.uint32(np.iinfo(np.uint32).max), dtype=np.uint32)
        allele_to_base = np.zeros(n_alleles, dtype=np.uint8)

        locus_to_alts: Dict[Tuple[str, int], Set[str]] = defaultdict(set)
        parsed_entries: List[Tuple[int, Tuple[str, int], str]] = []
        bad_loci: Set[Tuple[str, int]] = set()

        for idx, allele_id in enumerate(self.headers["alleles"]):
            try:
                chr_key, start, _, ref, alt = self._decode_allele_position(
                    allele_id=allele_id,
                    pattern=pattern,
                    allele_id_format=allele_id_format,
                )
            except ValueError as exc:
                log.debug(f"Skipping allele '{allele_id}' for SNV lookup: {exc}")
                continue

            if start is None:
                bad_loci.add((str(chr_key), None))
                continue

            locus_key = (str(chr_key), int(start))
            ref_base = (ref or "").upper()
            alt_base = (alt or "").upper()

            if not alt_base or len(alt_base) != 1 or alt_base not in base_map:
                bad_loci.add(locus_key)
                continue
            if ref_base and len(ref_base) != 1:
                bad_loci.add(locus_key)
                continue

            locus_to_alts[locus_key].add(alt_base)
            parsed_entries.append((idx, locus_key, alt_base))

        for locus_key, alt_set in locus_to_alts.items():
            if len(alt_set) > 1:
                bad_loci.add(locus_key)

        valid_loci = [lk for lk in locus_to_alts.keys() if lk not in bad_loci]
        valid_loci.sort(key=lambda item: (item[0], item[1]))
        locus_id_map: Dict[Tuple[str, int], int] = {
            locus: idx for idx, locus in enumerate(valid_loci)
        }

        for allele_idx, locus_key, alt_base in parsed_entries:
            locus_id = locus_id_map.get(locus_key)
            if locus_id is None:
                continue
            allele_to_locus[allele_idx] = np.uint32(locus_id)
            allele_to_base[allele_idx] = base_map[alt_base]

        snv_loci = len(valid_loci)
        self._snv_lookup_cache = (allele_to_locus, allele_to_base, snv_loci)
        self._snv_lookup_format = allele_id_format
        log.info(
            f"Constructed SNV lookup with {snv_loci} loci. "
            f"Skipped {len(bad_loci)} loci due to non-SNP, multi-allelic, or ambiguous bases."
        )
        return self._snv_lookup_cache


    @check_allele_id_format
    @check_bed_paths
    def compute_allele_positions( self,
                                  allele_id_format : str = "{chr}.{start}.{ref}.{alt}",
                                  allele_coords_bed : Union[str, None] = None,
                                  recompute_positions : bool = False
                                  ) -> Dict:
        """
        Compute or retrieve allele positions from either BED file or allele ID grammar.

        BED file takes precedence: any allele defined in the BED will not be parsed from ID.
        Positions are cached separately for BED and ID-derived alleles.
        
        Returns:
            Dict[str, Dict[int, List[str]]]: Dictionary of chromosome keys to position-to-allele mappings.
        """
        log.info(f"""Request to compute allele positions with allele_id_format={allele_id_format}, bed={allele_coords_bed}, recompute={recompute_positions}
                                    allele_coords_bed = {self._allele_coords_bed}
                                    allele_id_format = {self._allele_id_format}
                                    allele_positions_from_bed = {self._allele_positions_from_bed}
                                    allele_positions_from_ids = {self._allele_positions_from_ids}""")

        use_bed = bool(allele_coords_bed)

        ## reuse cached results if possible and no recomputation is requested
        if not recompute_positions:
            if self._allele_positions_from_bed and self._allele_positions_from_ids:
                log.debug("Using fully cached allele positions.")
                return self.allele_positions
            elif self._allele_positions_from_bed and not use_bed:
                log.debug("Using cached BED-based positions only.")
                return self.allele_positions
            elif self._allele_positions_from_ids and not use_bed:
                log.debug("Using cached ID-based positions only.")
                return self.allele_positions

        ## initialize or extend allele_positions and cache flags
        if recompute_positions:
            allele_positions = defaultdict(lambda: defaultdict(list))
            self._allele_positions_from_bed = False
            self._allele_positions_from_ids = False
        else:
            allele_positions = self.allele_positions

        bed_defined_alleles = set()

        ## decode BED if requested or not already cached
        if use_bed and (recompute_positions or not self._allele_positions_from_bed):
            log.debug("Decoding allele positions from BED file.")
            
            ## get the positions of each allele to check whether they are being redefined
            cached_allele_positions = self.get_allele_positions()
                        
            coords = self._read_allele_coords_bed(allele_coords_bed)
            for allele in self.headers["alleles"]:
                if allele in coords:
                    chr_key, start, end = coords[allele]
                    
                    ## raise a warning if an allele position is being redefined                    
                    if allele in cached_allele_positions:
                        cached_allele_pos = cached_allele_positions[allele]
                        if cached_allele_pos != [chr_key, start]:
                            log.warning(f"Position for allele {allele} being overwritten from {cached_allele_pos} to {[chr_key, start]}")
                    
                    allele_positions[chr_key][start].append(allele)
                    bed_defined_alleles.add(allele)
                    
            ## switch the bed flag
            self._allele_positions_from_bed = True
            self._allele_coords_bed = allele_coords_bed

        ## compute the regex pattern
        pattern = self._check_allele_format_grammar(allele_id_format=allele_id_format)
        
        ## decode from allele ID format if requested or not already cached
        if recompute_positions or not self._allele_positions_from_ids:
            log.debug("Decoding allele positions from ID format.")
            for allele in self.headers["alleles"]:
                if allele in bed_defined_alleles:
                    continue  ## skip if already defined via BED
                try:
                    chr_key, start, end, ref, alt = self._decode_allele_position(
                        allele_id=allele,
                        pattern=pattern,
                        allele_id_format=allele_id_format
                    )
                    if start is not None:
                        allele_positions[str(chr_key)][int(start)].append(allele)
                except ValueError as e:
                    log.debug(f"Skipping allele ID '{allele}': {e}")
            self._allele_positions_from_ids = True
            self._allele_id_format = allele_id_format

        self.allele_positions = allele_positions
        self._allele_position_map_cache = None
        self._site_count_all_cache = None
        return self.allele_positions


    def _decode_allele_position( self,
                                 allele_id : str,
                                 pattern : str,
                                 allele_id_format : str
                                 ) -> Tuple[ Union[str, None], int, Union[int, None], Union[str, None], str]:
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
        
        if not allele_id:
            raise ValueError("allele_id cannot be empty.")

        ## compile the regex and match against the allele ID
        match = re.compile(pattern).match(allele_id)

        if not match:
            raise ValueError(f"Allele ID '{allele_id}' does not match the format '{allele_id_format}'.")

        parts = match.groupdict()

        start = int(parts['start'])   ## exists or will raise error in previous checks
        end = int(parts['end']) if 'end' in parts else None

        return (
            parts.get('chr'),
            start,
            end,
            parts.get('ref'),
            parts.get('alt')
        )
    

    def _check_allele_format_grammar( self,
                                      allele_id_format : str
                                      ) -> str:
        """ Checks the allele id format grammar and constructs a pattern for regex-ing.
        """
        
        re = require_package("re", "re")
        
        accepted_placeholder_patterns = ["chr", "start", "end", "ref", "alt"]
        
        if not allele_id_format:
            raise AllelePatternError("allele_id_format cannot be empty.")
        
        ## check placeholders are valid
        pattern_placeholders = re.findall(r'\{([^}]+)\}', allele_id_format)
        for p in pattern_placeholders:
            if p not in accepted_placeholder_patterns:
                raise AllelePatternError(f"{p} not a valid placeholder for allele_id_format pattern. Accepted placeholders : {accepted_placeholder_patterns}")
            
        ## check minimum required format
        if "start" not in pattern_placeholders:
            raise AllelePatternError(f"{p} not a valid placeholder for allele_id_format pattern. Must have at least 'chr' and 'start' placeholders.")
        
        if "alt" not in pattern_placeholders and "ref" in pattern_placeholders:
            log.warning("allele_id_format : 'ref' provided but 'alt' absent. This may be a mistake.")
        
        pattern = re.escape(allele_id_format)
    
        if not pattern:
            raise AllelePatternError("allele_id_format could not be parsed. Please check your format string and ensure it contains valid placeholders.")
    
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
            raise AllelePatternError("allele_id_format does not contain any valid placeholders (e.g., {chr}, {start}).")

        ## anchor the pattern to match the entire string
        pattern = f"^{pattern}$"
        
        return pattern

    @check_bed_paths
    @check_allele_id_format
    def get_interval_alleles( self,
                              allele_id_format : Union[str, None],
                              intervals : Union[List, None] = None,
                              intervals_bed : Union[str, None] = None,
                              allele_coords_bed : Union[str, None] = None,
                              delimiter : str = "\t"
                              ) -> List[Set[str]]:
        bisect = require_package("bisect", "bisect")
        
        if allele_id_format == None:
            if self._allele_id_format == None:
                raise IntervalError("Please provide an allele id format")
            else:
                allele_id_format = self._allele_id_format
        
        ## check intervals were provided
        if not intervals and not intervals_bed:
            raise IntervalError("Please provide either a list of intervals using 'intervals' or the path to a bed file with interval positions using 'intervals_bed'.")
                
        ## get intervals
        cleaned_intervals = self._get_intervals(intervals=intervals,
                                                allele_id_format=allele_id_format,
                                                intervals_bed=intervals_bed,
                                                allele_coords_bed=allele_coords_bed,
                                                delimiter=delimiter)
        
        ## get allele positions
        allele_positions = self.compute_allele_positions( allele_id_format,
                                                          allele_coords_bed )

        ## find alleles which fall within the given intervals
        interval_alleles = []
            
        ## handle interval input
        for [chr_key, start, end] in cleaned_intervals:
            if chr_key not in allele_positions:
                log.warning(f"Chromosome {chr_key} not found in allele positions. Skipping...")
                continue

            chr_pos_dict = allele_positions[chr_key]
            if not chr_pos_dict:
                log.warning(f"No alleles found for chromosome {chr_key}. Skipping...")
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
            log.warning(f"No alleles found within the given intervals. Please check your intervals and that allele IDs are in the expected form defined by allele_id_format (currently {allele_id_format}).")
            
        log.info(f"Found {len(interval_alleles)} alleles.")
        return list(set(interval_alleles))


    def _get_intervals( self,
                        allele_id_format : str,
                        intervals : Union[List, None],
                        intervals_bed : Union[str, None] = None,
                        allele_coords_bed : Union[str, None] = None,
                        delimiter : str = " "
                        ) -> List[List]:
        log.info(f"Generating a list of intervals.")
        
        ## get intervals passed from list
        cleaned_intervals_list = self._get_intervals_from_list(intervals)
        
        ## get intervals from bed file
        cleaned_intervals_bed = self._get_intervals_from_bed(intervals_bed=intervals_bed,
                                                             delimiter=delimiter)
                    
        ## concatenate intervals from list and bed
        intervals_list = self._make_intervals_list(intervals_from_list=cleaned_intervals_list,
                                                   intervals_from_bed=cleaned_intervals_bed)
        
        return intervals_list
    
    
    def _get_intervals_from_list( self,
                                  intervals : Union[List, None]
                                  ) -> List:
        if intervals != None:
            log.info(f"Found {len(intervals)} records in 'intervals'.")
            interval_list = self._check_intervals(intervals=intervals,
                                                  from_bed=False)
            log.info(f"Found {len(interval_list)} valid intervals in 'intervals'.")
        else:
            interval_list = []
            log.info(f"No intervals passed using 'intervals'.")
        
        return interval_list
        
    
    def _get_intervals_from_bed( self,
                                 intervals_bed : Union[str, None], 
                                 delimiter : str = " ",
                                 ) -> List:
        """ Read the coords BED file containing genomic intervals.
        Must be structures as:
        chr start end
        """
        cleaned_intervals = []
        
        log.debug(f"Intervals bed path : {intervals_bed}")
         
        if intervals_bed != None:
            with open(intervals_bed, 'r') as fin:
                lines = fin.readlines()
            
            log.info(f"Found {len(lines)} entries in {intervals_bed}")
            
            ## construct a generic list of bed entries
            intervals = [l.strip("\n").split() for l in lines]
            cleaned_intervals = self._check_intervals(intervals=intervals,
                                                      from_bed=True)
            log.info(f"Found {len(cleaned_intervals)} valid intervals in {intervals_bed}")
        
        return cleaned_intervals


    @staticmethod
    def _check_intervals( intervals : list,
                          from_bed : bool = False,
                          ) -> List[Tuple[str, int, int]]:
        ## clean intervals
        cleaned_intervals = []
            
        for i, interval in enumerate(intervals):
            if not isinstance(interval, list) and not isinstance(interval, tuple):
                if not from_bed:
                    raise IntervalError(f"intervals should be list of lists/tuples like [ [$interval_1], [$interval_2], ... ]")
                else:
                    raise IntervalError(f"bed file contains malformed lines. ")
            try:
                if len(interval) == 3:
                    chr_key = str(interval[0])
                    start = int(interval[1])
                    end = int(interval[2])
                elif len(interval) == 2:
                    chr_key = "None"
                    start = int(interval[0])
                    end = int(interval[1])
                elif len(interval) == 1:
                    chr_key = str(interval[0])
                    start = None
                    end = None
                else:
                    if not from_bed: 
                        raise IntervalError(f"Malformed record at {i}. Intervals must be either [chr_key, start, end], [start, end] or [chr_key] : not {interval}")
                    else:
                        raise IntervalError(f"Malformed line at {i}. Bed lines must be either $chr_key$, $chr_key start end$ or $start end$ : not {interval}")
                    
            except (ValueError, TypeError, IndexError) as e:
                raise IntervalError(f"Failed to parse interval {interval} : {e}")
                        
            cleaned_intervals.append([chr_key, start, end])
            
        return cleaned_intervals
    
    
    @staticmethod
    def _make_intervals_list( intervals_from_list : List,
                              intervals_from_bed : List
                              ) -> List:
        ## intervals only from list
        if len(intervals_from_list) > 0 and len(intervals_from_bed) == 0:
            return intervals_from_list
        
        ## intervals only from bed
        elif len(intervals_from_list) == 0 and len(intervals_from_bed) > 0:
            return intervals_from_bed
        
        ## no intervals for some reason
        if len(intervals_from_list) == 0 and len(intervals_from_bed) == 0:
            raise IntervalError("No valid intervals found.")
        
        ## intervals from both are assumed to be of the same structure due to passing through
        ## _check_intervals and so can be safely concatenated
        else:
            all_intervals = intervals_from_list + intervals_from_bed
            log.info(f"Total number of intervals from all sources : {len(all_intervals)}")
            return all_intervals
        

    @staticmethod
    def _read_allele_coords_bed( allele_coords_bed : str 
                                 ) -> Dict:
        """ Read the coords BED file containing position : allele_id mappings.
        Must be structures as:
        chr start end allele_id
        """
        coords = {}
        with open(allele_coords_bed, 'r') as fin:
            lines = fin.readlines()
            
        for i, line in enumerate(lines):
            parts = line.strip().split()
            if len(parts) != 4:
                log.warning(f"Skipping malformed line {i+1} in {allele_coords_bed}: '{line.strip()}'")
                continue
            
            chrom = str(parts[0])
            start = int(parts[1])
            end = int(parts[2])
            allele_id = str(parts[3])
            if allele_id in coords:
                log.warning(f"Duplicate allele_id '{allele_id}' found in {allele_coords_bed}. Overwriting previous entry.")
            coords[allele_id] = [chrom, start, end]
    
        return coords

    
    @check_guids_list
    def check_guids( self,
                     guids : list,
                     filter : bool = False
                     ) -> Union[List[str], None]:
        """ Check guids are present within the matrix and construct a list of present guids to proceed with
        if filter is True
        """
        absent_guids = set()
        present_guids = set()
        for id in guids:
            if id not in self.headers["guids"]:
                absent_guids.add(id)
            else:
                present_guids.add(id)

        if filter:
            return list(present_guids)

        if len(absent_guids) > 0:
            raise InvalidGUIDQueryError(f"guids {absent_guids} not found in allele matrix.")


    @check_alleles_list
    def check_alleles( self,
                       alleles : list,
                       filter : bool = False
                       ) -> Union[List[str], None]:
        """ Check alleles are present within the matrix and construct a list of present alleles to proceed with
        if filter is True
        """
        absent_alleles = []
        present_alleles = []
        for allele_id in alleles:
            if allele_id not in self.headers["alleles"]:
                absent_alleles.append(allele_id)
            else:
                present_alleles.append(allele_id)

        if filter:
            return present_alleles
        
        if len(absent_alleles) > 0:
            raise InvalidAlleleQueryError(f"alleles {absent_alleles} not found in allele matrix.")

    
    def get_allele_positions( self
                              ) -> Dict:
        """ Returns allele positions by inverting the allele_positions dictionary
        """
        if self._allele_position_map_cache is None:
            self._allele_position_map_cache = {
                allele: (chr_key, pos)
                for chr_key, pos_dict in self.get_allele_posmap().items()
                for pos, alleles in pos_dict.items()
                for allele in alleles
            }
        return self._allele_position_map_cache


    def _init_missing_masks(self, missing_masks: Union[Dict, None]) -> None:
        guid_list = self.headers.get("guids", []) if isinstance(self.headers, dict) else []
        per_guid = {guid: [] for guid in guid_list}

        if not missing_masks or not isinstance(missing_masks, dict):
            self._missing_masks = per_guid
            return

        raw_per_guid = missing_masks.get("column_masks") or missing_masks.get("missing_masks")
        if isinstance(raw_per_guid, dict):
            for guid, sites in raw_per_guid.items():
                if guid not in per_guid:
                    continue
                if isinstance(sites, (list, tuple, set)):
                    per_guid[guid] = list(sites)
                elif sites is None:
                    per_guid[guid] = []
                else:
                    per_guid[guid] = [sites]

        self._missing_masks = per_guid


    def has_missing_mask(self) -> bool:
        return any(self._missing_masks.values())


    def get_guid_missing_mask(self, guid: str) -> List:
        return self._missing_masks.get(guid, [])


    def get_missing_mask_rows(self) -> Union[List[List[int]], None]:
        """
        Return a row-aligned list of missing-column indices for each GUID.
        """
        if not self.has_missing_mask():
            return None
        if self._missing_mask_rows is not None:
            return self._missing_mask_rows

        rows: List[List[int]] = []
        for guid in self.headers.get("guids", []):
            cols = self._missing_masks.get(guid, [])
            rows.append(list(cols))

        self._missing_mask_rows = rows
        return self._missing_mask_rows


    def count_sites_for_alleles(self,
                                alleles: List[str],
                                allele_id_format: Union[str, None] = None) -> int:
        """
        Count unique genomic sites represented by the provided allele IDs.

        Args:
            alleles: Allele identifiers to include in the site count.
            allele_id_format: Optional override for decoding allele positions if they
                have not already been computed.

        Returns:
            Number of unique (chromosome, position) sites covered by the alleles.
        """
        if alleles is None:
            return 0

        if not isinstance(alleles, list):
            alleles = list(alleles)
        validate_type(alleles, list, "alleles")

        if alleles is self.headers.get("alleles") and self._site_count_all_cache is not None:
            return self._site_count_all_cache

        have_positions = bool(self._allele_positions_from_bed or self._allele_positions_from_ids)
        fmt_to_use = allele_id_format or self._allele_id_format

        if not have_positions:
            if fmt_to_use is None:
                raise ParameterError(
                    "Normalizing distances requires allele positions. Provide 'allele_positions' "
                    "or supply 'allele_id_format' so positions can be decoded."
                )
            self.ensure_id_positions(fmt_to_use)
        elif allele_id_format and fmt_to_use != self._allele_id_format:
            self.ensure_id_positions(allele_id_format, force=True)

        allele_position_map = self.get_allele_positions()
        unique_sites = set()
        missing_alleles = []

        for allele in alleles:
            coord = allele_position_map.get(allele)
            if coord is None:
                missing_alleles.append(allele)
                continue
            chr_key, pos = coord
            if chr_key is None or pos is None:
                continue
            unique_sites.add((str(chr_key), int(pos)))

        if missing_alleles:
            log.warning(
                "Could not determine positions for %d allele(s); they were excluded from the site count.",
                len(missing_alleles),
            )

        site_count = len(unique_sites)
        if alleles is self.headers.get("alleles"):
            self._site_count_all_cache = site_count
        return site_count


    def has_id_positions(self) -> bool:
        """Return True if allele positions have been decoded from allele IDs."""
        return bool(self._allele_positions_from_ids and self._allele_id_format)


    def get_cached_allele_id_format(self) -> Union[str, None]:
        """Return the cached allele ID format used for allele ID decoding."""
        return self._allele_id_format


    def ensure_id_positions(self,
                            allele_id_format: str,
                            *,
                            force: bool = False) -> None:
        """
        Ensure allele positions are decoded from allele IDs for the given format.

        Args:
            allele_id_format: The allele ID grammar to use for decoding.
            force: When True, recompute positions even if already cached for the same format.
        """
        validate_type(allele_id_format, str, "allele_id_format")

        recompute = force or (self._allele_positions_from_ids and self._allele_id_format != allele_id_format)
        if recompute or not self._allele_positions_from_ids or self._allele_id_format != allele_id_format:
            self.compute_allele_positions(
                allele_id_format=allele_id_format,
                allele_coords_bed=self._allele_coords_bed,
                recompute_positions=recompute
            )
        
        
    def get_allele_posmap( self
                           ) -> Dict:
        """ Get the allele_positions dictionary.
        """
        return self.allele_positions


    def encode_guid( self,
                     guid : str
                     ) -> int:
        """ Encode a GUID to its corresponding index in the headers dictionary.
        """
        validate_type(guid, str, "guid")
        return self._decoded_headers["guids"][guid]


    def decode_guid( self,
                     row_coord : int
                     ) -> str:
        """ Decode a row coordinate to its corresponding GUID in the headers dictionary.
        """
        validate_type(row_coord, Union[int, np.int8, np.int16, np.int32, np.int64], "row_coord")
        return self.headers["guids"][row_coord]


    def encode_allele( self,
                       allele : str
                       ) -> int:
        """ Encode an allele to its corresponding index in the headers dictionary.
        """
        validate_type(allele, str, "allele")
        return self._decoded_headers["alleles"][allele]


    def decode_allele( self,
                       col_coord : int
                       ) -> str:
        """ Decode a column coordinate to its corresponding allele in the headers dictionary.
        """
        validate_type(col_coord, Union[int, np.int8, np.int16, np.int32, np.int64, np.uint, np.uint8, np.uint16, np.uint32, np.uint64], "col_coord")
        return self.headers["alleles"][col_coord]
