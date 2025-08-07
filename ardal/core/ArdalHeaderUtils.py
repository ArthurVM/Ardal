import os
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Union, List, Tuple, Set, Dict

from ..utils.misc import require_package
from ..utils.decorators import check_allele_id_format, check_bed_paths, check_alleles_list, check_guids_list, validate_type
from ..utils.exceptions import AllelePatternError, IntervalError, InvalidGUIDQueryError, InvalidAlleleQueryError
from ..utils.logger import get_logger

log = get_logger()


# core/ArdalHeaderUtils.py
class ArdalHeaderUtils:

    def __init__( self,
                  headers : dict,
                  allele_coords_bed : Union[str, None] = None,
                  allele_id_format : Union[str, None] = None,
                  allele_positions : Union[Dict, None] = None ):

        self.headers = headers
        self._decoded_headers = self._decode_headers()
        
        ## allele positions were not provided
        if not allele_positions:
            self.allele_positions = defaultdict(lambda: defaultdict(list))
            
            ## position cache switches
            self._allele_positions_from_bed = False   ## positions provided from bed
            self._allele_positions_from_ids = False   ## positions decoded from allele_ids
            
            ## if allele_id_format is provided then just decode positions now
            if allele_id_format:
                self.allele_positions = self.compute_allele_positions(allele_id_format = allele_id_format,
                                                                      allele_coords_bed = allele_coords_bed)
                self._allele_coords_bed = allele_coords_bed
                self._allele_id_format = allele_id_format
        
        ## allele positions were provided in a dict
        else:
            ## TODO: some checks on this
            self.allele_positions = allele_positions
            self._allele_positions_from_bed = False
            self._allele_positions_from_ids = True


    def _decode_headers( self ) -> Dict:
        return  { "guids" : dict(zip(self.headers["guids"], range(len(self.headers["guids"])))),
                  "alleles" : dict(zip(self.headers["alleles"], range(len(self.headers["alleles"]))))
                }


    @check_allele_id_format
    @check_bed_paths
    def compute_allele_positions( self,
                                  allele_id_format: str = "{chr}.{start}.{ref}.{alt}",
                                  allele_coords_bed: Union[str, None] = None,
                                  recompute_positions: bool = False
                                  ) -> Dict:
        """
        Compute or retrieve allele positions from either BED file or allele ID grammar.

        BED file takes precedence: any allele defined in the BED will not be parsed from ID.
        Positions are cached separately for BED and ID-derived alleles.
        
        Returns:
            Dict[str, Dict[int, List[str]]]: Dictionary of chromosome keys to position-to-allele mappings.
        """
        log.debug(f"Request to compute allele positions with id_format={allele_id_format}, bed={allele_coords_bed}, recompute={recompute_positions}")

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
            coords = self._read_allele_coords_bed(allele_coords_bed)
            for allele in self.headers["alleles"]:
                if allele in coords:
                    chr_key, start, end = coords[allele]
                    allele_positions[str(chr_key)][int(start)].append(allele)
                    bed_defined_alleles.add(allele)
            self._allele_positions_from_bed = True
            self._allele_coords_bed = allele_coords_bed

        ## decode from allele ID format if requested or not already cached
        if recompute_positions or not self._allele_positions_from_ids:
            log.debug("Decoding allele positions from ID format.")
            for allele in self.headers["alleles"]:
                if allele in bed_defined_alleles:
                    continue  ## skip if already defined via BED
                try:
                    chr_key, start, end, ref, alt = self._decode_allele_position(
                        allele_id=allele,
                        allele_id_format=allele_id_format
                    )
                    if chr_key is not None and start is not None:
                        allele_positions[str(chr_key)][int(start)].append(allele)
                except ValueError as e:
                    log.debug(f"Skipping allele ID '{allele}': {e}")
            self._allele_positions_from_ids = True
            self._allele_id_format = allele_id_format

        self.allele_positions = allele_positions
        return self.allele_positions


    def _decode_allele_position( self,
                                 allele_id : str,
                                 allele_id_format : str = "{chr}.{start}.{ref}.{alt}"
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
        
        pattern = self._check_allele_format_grammar(allele_id_format=allele_id_format)
        
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
        if "start" not in pattern_placeholders or "alt" not in pattern_placeholders:
            raise AllelePatternError(f"{p} not a valid placeholder for allele_id_format pattern. Must have at least 'chr' and 'start' placeholders.")
        
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
                              intervals : Union[List, None] = None,
                              intervals_bed : Union[str, None] = None,
                              allele_coords_bed : Union[str, None] = None,
                              allele_id_format : str = "{chr}.{start}.{ref}.{alt}",
                              delimiter : str = " "
                              ) -> List[Set[str]]:
        bisect = require_package("bisect", "bisect")
        
        ## check intervals were provided
        if not intervals and not intervals_bed:
            raise IntervalError("Please provide either a list of intervals using 'intervals' or the path to a bed file with interval positions using 'intervals_bed'.")
                
        ## get intervals
        cleaned_intervals = self._get_intervals(intervals=intervals,
                                                intervals_bed=intervals_bed,
                                                allele_coords_bed=allele_coords_bed,
                                                allele_id_format=allele_id_format,
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
            
        log.debug(f"Found {len(interval_alleles)} alleles.")
        return list(set(interval_alleles))


    def _get_intervals( self,
                        intervals : Union[List, None],
                        intervals_bed : Union[str, None] = None,
                        allele_coords_bed : Union[str, None] = None,
                        allele_id_format : str = "{chr}.{start}.{ref}.{alt}",
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
            log.debug(f"Found {len(intervals)} records in 'intervals'.")
            interval_list = self._check_intervals(intervals=intervals,
                                                  from_bed=False)
            log.debug(f"Found {len(interval_list)} valid intervals in 'intervals'.")
        else:
            interval_list = []
            log.debug(f"No intervals passed using 'intervals'.")
        
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
        
        log.debug(f"Bed path : {intervals_bed}")
         
        if intervals_bed != None:
            with open(intervals_bed, 'r') as fin:
                lines = fin.readlines()
            
            log.debug(f"Found {len(lines)} entries in {intervals_bed}")
            
            ## construct a generic list of bed entries
            intervals = [l.strip("\n").split(delimiter) for l in lines]
            cleaned_intervals = self._check_intervals(intervals=intervals,
                                                      from_bed=True)
            log.debug(f"Found {len(cleaned_intervals)} valid intervals in {intervals_bed}")
        
        return cleaned_intervals


    @staticmethod
    def _check_intervals( intervals : list,
                          from_bed : bool = False,
                          ) -> List[Tuple[str, int, int]]:
        ## clean intervals
        cleaned_intervals = []
            
        for i, interval in enumerate(intervals):
            if not isinstance(interval, list):
                if not from_bed:
                    raise IntervalError(f"intervals should be a nested list like [ [$interval_1], [$interval_2], ... ]")
                else:
                    raise IntervalError(f"bed file contains malformed lines. ")
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
            log.debug(f"Total number of intervals from all sources : {len(all_intervals)}")
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
            
            chrom, start, end, allele_id = parts
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
        absent_guids = []
        present_guids = []
        for id in guids:
            if id not in self.headers["guids"]:
                absent_guids.append(id)
            else:
                present_guids.append(id)

        if filter:
            return present_guids

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

    
    def get_allele_positions( self ) -> Dict:
        """ Returns the allele_positions dictionary.
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
        validate_type(row_coord, int, "row_coord")
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
        validate_type(col_coord, int, "col_coord")
        return self.headers["alleles"][col_coord]