import os
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Union

from .utilities import *
from ..exceptions.exceptions import *


# core/ArdalHeaderUtils.py
class ArdalHeaderUtils:

    def __init__( self,
                  headers : dict,
                  coords_bed : Union[str, None] = None,
                  allele_id_format : Union[str, None] = None ):

        self.headers = headers
        self.pos_decoded_bool = False
        self._decoded_headers = self._decode_headers()
        self.allele_positions = {}
        
        ## if allele_id_format is provided then just decode positions now
        if allele_id_format:
            self.allele_positions = self.getAllelePositions(allele_id_format=allele_id_format, coords_bed=coords_bed)
            self._coords_bed = coords_bed
            self._allele_id_format = allele_id_format


    def _decode_headers( self ) -> dict:
        return  { "guids" : dict(zip(self.headers["guids"], range(len(self.headers["guids"])))),
                  "alleles" : dict(zip(self.headers["alleles"], range(len(self.headers["alleles"]))))
                }


    def getAllelePositions( self,
                            allele_id_format : str = "{chr}.{start}.{ref}.{alt}",
                            coords_bed : Union[str, None] = None,
                            recompute_positions : bool = False ) -> dict:
        # print(f"Decoding allele positions with allele_id_format : {allele_id_format}")
        ## check that the allele_positions havent already been computed
        ## if not then compute them
        if not self.pos_decoded_bool or recompute_positions:
            print("Computing allele positions.")
            if coords_bed:
                coords = self._readCoordsBED(coords_bed)
            else:
                coords = None

            ## get the position of each allele
            allele_positions = defaultdict(lambda: defaultdict(list))
            for allele in self.headers["alleles"]:
                try:
                    ## if the bed was provided then just get positions from that
                    if coords:
                        if allele not in coords:
                            continue
                        chr_key, start, end = coords[allele]
                        allele_positions[str(chr_key)][int(start)].append(allele)
                        
                    ## otherwise compute them using allele_IDs and the allele_id_format grammar
                    else:
                        chr_key, start, end, ref, alt = self._decodeAllelePosition(allele, allele_id_format)
                        if start is not None:
                            allele_positions[chr_key][start].append(allele)
                            
                except ValueError as e:
                    ardal_warn(f"Could not parse allele ID '{allele}': {e}")
                    continue
            
            ## assign for future calls
            self.allele_positions = allele_positions
        
        else:
            print("Found precomputed allele positions.")
        self.pos_decoded_bool = True
        return self.allele_positions
        

    def _decodeAllelePosition( self,
                              allele_id : str,
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
        
        pattern = self._checkAlleleFormatGrammar(allele_id_format=allele_id_format)
        
        if not allele_id:
            raise ValueError("allele_id cannot be empty.")

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
        

    def _checkAlleleFormatGrammar( self,
                                   allele_id_format : str ) -> str:
        """ Checks the allele id format grammar and constructs a pattern for regex-ing.
        """
        
        re = require_package("re", "re")
        
        accepted_placeholder_patterns = ["chr", "start", "end", "ref", "alt"]
        
        if not allele_id_format:
            raise ValueError("allele_id_format cannot be empty.")
        
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
        
        return pattern


    def getIntervalAlleles( self,
                            intervals : list,
                            allele_id_format : str = "{chr}.{start}.{ref}.{alt}",
                            coords_bed : Union[str, None] = None ) -> list:
        bisect = require_package("bisect", "bisect")
        
        ## check and clean intervals
        cleaned_intervals = self._checkIntervals(intervals)
        
        ## get allele positions
        allele_positions = self.getAllelePositions( allele_id_format,
                                                    coords_bed )

        ## find alleles which fall within the given intervals
        interval_alleles = []
            
        ## handle interval input
        for [chr_key, start, end] in cleaned_intervals:
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
            
        return list(set(interval_alleles))


    def _checkIntervals( self,
                         intervals : list ) -> list:
        ## clean intervals
        cleaned_intervals = []
        for interval in intervals:
            if not isinstance(interval, list):
                raise IntervalError(f"intervals parameter should be a nested list like [ [$interval_1], [$interval_2], ... ]")
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
            
            cleaned_intervals.append([chr_key, start, end])
            
        return cleaned_intervals
            
            
    def _readCoordsBED( self,
                       coords_bed : str ) -> dict:
        """ Read the coords BED file containing position : allele_id mappings.
        Must be structures as:
        chr start end allele_id
        """
        coords = {}
        with open(coords_bed, 'r') as fin:
            for i, line in enumerate(fin):
                parts = line.strip().split()
                if len(parts) != 4:
                    ardal_warn(f"Skipping malformed line {i+1} in {coords_bed}: '{line.strip()}'")
                    continue
                
                chrom, start, end, allele_id = parts
                if allele_id in coords:
                    ardal_warn(f"Duplicate allele_id '{allele_id}' found in {coords_bed}. Overwriting previous entry.")
                coords[allele_id] = [chrom, start, end]
        
        return coords


    def checkGUIDs( self,
                    guids : list,
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
            if id not in self.headers["guids"]:
                absent_guids.append(id)
            else:
                present_guids.append(id)

        if filter:
            return present_guids

        if len(absent_guids) > 0:
            raise InvalidGUIDQueryError(f"guids {absent_guids} not found in allele matrix.")


    def checkAlleles( self,
                      alleles : list,
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
            if allele_id not in self.headers["alleles"]:
                absent_alleles.append(allele_id)
            else:
                present_alleles.append(allele_id)

        if filter:
            return present_alleles
        
        if len(absent_alleles) > 0:
            raise InvalidAlleleQueryError(f"alleles {absent_alleles} not found in allele matrix.")



    def encodeGuid( self, guid : str ) -> int:
        """ Encode a GUID to its corresponding index in the headers dictionary.
        """
        return self._decoded_headers["guids"][guid]


    def decodeGuid( self, row_coord : int ) -> str:
        """ Decode a row coordinate to its corresponding GUID in the headers dictionary.
        """
        return self.headers["guids"][row_coord]


    def encodeAllele( self, allele : str ) -> int:
        """ Encode an allele to its corresponding index in the headers dictionary.
        """
        return self._decoded_headers["alleles"][allele]


    def decodeAllele( self, col_coord : int ) -> str:
        """ Decode a column coordinate to its corresponding allele in the headers dictionary.
        """
        return self.headers["alleles"][col_coord]