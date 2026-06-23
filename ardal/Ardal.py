""" Frontend for the Ardal allele matrix toolkit.
"""
import logging
import numpy as np
import pandas as pd
from typing import Union, Dict, Tuple

from .core.ArdalParser import ArdalParser
from .core.ArdalHeaderUtils import ArdalHeaderUtils
from .core.ArdalIO import ArdalIO
from .core.ArdalGet import ArdalGet
from .core.ArdalAllele import ArdalAllele
from .core.ArdalDistance import ArdalDistance
from .core.ArdalRecombination import ArdalRecombination
from .core.ArdalStats import ArdalStats
from .utils.decorators import check_density_threshold, check_Ardal_roaring_param, check_allele_id_format, check_bed_paths
from .utils.exceptions import MatrixParseError
from .utils.logger import get_logger

log = get_logger()
SILENT = logging.CRITICAL + 1
VERBOSITY_LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warn": logging.WARNING,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL,
    "silent": SILENT
}

## Ardal.py
class Ardal(object):
    """ Frontend class for managing Ardal allele matrices, providing interfaces for IO, querying,
    statistics, and transformation operations on binary allele matrices using Roaring or Bit-packed backends.
    """
    
    @check_density_threshold
    @check_Ardal_roaring_param
    @check_allele_id_format
    @check_bed_paths
    def __init__( self,
                  data_source : Union[str, pd.DataFrame, Tuple[str, str], Tuple[np.ndarray, Dict]],
                  roaring : Union[str, bool] = "auto",
                  density_threshold : float = 0.20,
                  allele_id_format : Union[str, None] = None,
                  allele_coords_bed : Union[str, None] = None,
                  verbosity : Union[int, str] = "error",
                  quiet_init : bool = False,
                  allele_positions : Union[Dict, None] = None ):
        """ Initializes an Ardal object from a given allele matrix source. Selects appropriate backend
        based on matrix density and user preference. Also sets up downstream API components.

        Args:
            data_source (str): Path to matrix input file or in-memory source.
            roaring (Union[str, bool]): 'auto', 'true', or 'false' to control backend selection.
            density_threshold (float): Threshold to control backend switch-over.
            allele_id_format (str or None): Optional allele ID format specifier.
            allele_coords_bed (str or None): Optional BED file to map allele positions.
            verbosity (int or str): Logging verbosity level.
            quiet_init (bool): If True, suppress initial matrix stats printout.
            allele_positions (dict or none): a dictionary of allele positions which can be passed to headerUtils.
        """
        from . import backends
        
        if isinstance(verbosity, str):
            ## set backend verbosity here to avoid double ardal import and circular issues
            backends.ardal.set_backend_verbosity(verbosity.lower())
            verbosity = VERBOSITY_LEVELS.get(verbosity.lower(), logging.WARNING)

        log.setLevel(verbosity)

        self._hybrid_matrix = None
        self._headers = None
        self._missing_sites = {"site_keys": {}, "per_guid": {}}
        is_packed_mem = False
        if isinstance(data_source, list):
            for item in data_source:
                if self._looks_like_bitpacked_matrix(item):
                    is_packed_mem = True
                    break

        parser = ArdalParser(data_source, is_packed_mem=is_packed_mem)

        ## do some parameter fiddling to enforce the generation of a roaring matrix
        self._roaring_setter(roaring, density_threshold)
                    
        if parser.matrix is None:
            ## raise an error if parsing fails to prevent unexpected behaviour down the line
            raise MatrixParseError(f"Failed to parse data from: {data_source}") 

        self._headers = parser.headers
        self._meta = parser.meta
        self._missing_masks = parser.missing_masks

        log.info("Initialising Ardal component classes.")
        self._headerUtils = ArdalHeaderUtils(headers = self._headers,
                                             meta = self._meta,
                                             allele_coords_bed = allele_coords_bed,
                                             allele_id_format = allele_id_format,
                                             allele_positions = allele_positions,
                                             missing_masks = self._missing_masks)

        if hasattr(self._missing_masks, "to_backend_payload"):
            mask_arg = self._missing_masks.to_backend_payload()
        else:
            mask_rows = self._headerUtils.get_missing_mask_rows()
            mask_arg = mask_rows if mask_rows else None

        log.info(f"""Initialising ardal::HybridMatrix backend with:
                                matrix = {parser.matrix.shape} {type(parser.matrix)}
                                is_bitpacked = {parser.is_bitpacked} {type(parser.is_bitpacked)}
                                n_cols_bits = {parser.meta["n_cols"]} {type(parser.meta["n_cols"])}
                                build_roaring = {self.build_roaring} {type(self.build_roaring)}
                                density_threshold = {self.density_threshold} {type(self.density_threshold)}""")
        self._hybrid_matrix = backends.ardal.HybridMatrix( parser.matrix,
                                                           is_bitpacked=parser.is_bitpacked,
                                                           n_cols_bits=parser.meta["n_cols"],
                                                           use_roaring_if_sparse=self.build_roaring,
                                                           density_threshold=self.density_threshold,
                                                           missing_mask=mask_arg )
        
        ## set this as a signal from the backend
        self.roaring = self._hybrid_matrix.roaringEnabled()
        
        self.io = ArdalIO(headerUtils = self._headerUtils,
                          hybrid_matrix = self._hybrid_matrix,
                          roaring_enabled = self.roaring)
        
        self.get = ArdalGet(headerUtils = self._headerUtils,
                            hybrid_matrix = self._hybrid_matrix,
                            roaring_enabled = self.roaring)
        
        self.allele = ArdalAllele(headerUtils = self._headerUtils,
                                  hybrid_matrix = self._hybrid_matrix,
                                  roaring_enabled = self.roaring)
        
        self.distance = ArdalDistance(headerUtils = self._headerUtils,
                                      hybrid_matrix = self._hybrid_matrix,
                                      roaring_enabled = self.roaring)
        
        self.stats = ArdalStats(headerUtils = self._headerUtils,
                                hybrid_matrix = self._hybrid_matrix,
                                roaring_enabled = self.roaring)

        self.recomb = ArdalRecombination(headerUtils = self._headerUtils,
                                         hybrid_matrix = self._hybrid_matrix,
                                         roaring_enabled = self.roaring)
        
        if not quiet_init and verbosity < SILENT:
            self.get.matrix_stats(print_table=True)
            
    
    def _roaring_setter( self,
                         roaring : str,
                         density_threshold : float
                         ) -> None:
        """ Interprets the `roaring` user parameter and determines whether to build
        a RoaringBitmap backend or a Bit-packed backend, along with appropriate
        density threshold.

        Args:
            roaring (str): 'auto', 'true', or 'false'.
            density_threshold (float): Threshold value between 0 and 1.

        Raises:
            ValueError: If roaring mode string is invalid.
        """
        if roaring == "true":
            ## do some parameter fiddling to enforce the generation of a roaring matrix
            self.build_roaring = True
            self.density_threshold = 1.0
        elif roaring == "auto":
            ## let Ardal decide
            ## TODO the whole backend dynamic profiling stuff
            self.build_roaring = True
            self.density_threshold = density_threshold
        elif roaring == "false":
            ## dont build roaring backend
            self.build_roaring = False
            self.density_threshold = density_threshold
        else:
            raise ValueError(f"Invalid roaring mode: {roaring}")
        
        log.debug(f"Roaring setter : roaring={roaring}; density_threshold={self.density_threshold}; build_roaring={self.build_roaring}")


    @staticmethod
    def _looks_like_bitpacked_matrix(obj) -> bool:
        """Heuristic to detect in-memory bitpacked matrices."""
        if not isinstance(obj, np.ndarray) or obj.ndim != 2:
            return False
        dt = obj.dtype
        return dt.kind == "u" and dt.itemsize == 8


    def set_verbosity( self,
                       level_param: Union[int, str]
                       ) -> None:
        """ Update logging level post hoc for Ardal and all downstream components.

        Args:
            level (int or str): Logging level name or constant (e.g., 'debug', logging.INFO).
        """
        if isinstance(level_param, str):
            level = VERBOSITY_LEVELS.get(level_param.lower(), logging.WARNING)
        log.info(f"Switching verbosity to {level_param}.")  
        log.setLevel(level)
