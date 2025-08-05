""" Frontend for the Ardal allele matrix toolkit.
"""
from typing import Union

from .core.ArdalParser import ArdalParser
from .core.ArdalHeaderUtils import ArdalHeaderUtils
from .core.ArdalIO import ArdalIO
from .core.ArdalGet import ArdalGet
from .core.ArdalAllele import ArdalAllele
from .core.ArdalCompare import ArdalCompare
from .core.ArdalStats import ArdalStats
from .core.utilities import *
from .exceptions.exceptions import *


## Ardal.py
class Ardal(object):
    """ Class for handling Ardal-neighbourhood objects. """


    def __init__( self,
                  data_source : str,
                  coords_bed : Union[str, None] = None,
                  allele_id_format : Union[str, None] = None,
                  roaring : Union[str, bool] = "auto",
                  density_threshold : float = 0.02,
                  file_format : Union[str, None] = None,
                  quiet : bool = False ):
        """ Ardal constructor
        """
        from . import _ardal

        super(Ardal, self).__init__()

        self._hybrid_matrix = None
        self._headers = None

        parser = ArdalParser(data_source, file_format)

        ## do some parameter fiddling to enforce the generation of a roaring matrix
        roaring_param = validate_Ardal_roaring_param(roaring)
        if roaring_param == "true":
            use_roaring_if_sparse = True
            density_threshold = 1.0
        elif roaring_param == "auto":
            use_roaring_if_sparse = True
        elif roaring_param == "false":
            use_roaring_if_sparse = False
            
        if parser.matrix is not None:
            self._hybrid_matrix = _ardal.HybridMatrix(parser.matrix,
                                                       use_roaring_if_sparse=use_roaring_if_sparse,
                                                       density_threshold=density_threshold)
            self._headers = parser.headers
        else:
            ## raise an error if parsing fails to prevent unexpected behaviour down the line
            raise MatrixParseError(f"Failed to parse data from: {data_source}") 
        
        self.roaring = self._hybrid_matrix.roaringEnabled()
        
        self.headerUtils = ArdalHeaderUtils(self._headers, coords_bed, allele_id_format)
        self.io = ArdalIO(self.headerUtils, self._hybrid_matrix, self.roaring)
        self.get = ArdalGet(self.headerUtils, self._hybrid_matrix, self.roaring)
        self.allele = ArdalAllele(self.headerUtils, self._hybrid_matrix, self.roaring)
        self.compare = ArdalCompare(self.headerUtils, self._hybrid_matrix, self.roaring)
        self.stats = ArdalStats(self.headerUtils, self._hybrid_matrix, self.roaring)

        if not quiet:
            self.get.matrixStats(print_stats=True)