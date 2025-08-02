from .core.ArdalParser import ArdalParser
from .core.ArdalIO import ArdalIO
from .core.ArdalGet import ArdalGet
from .core.ArdalAllele import ArdalAllele
from .core.ArdalCompare import ArdalCompare
from .core.ArdalStats import ArdalStats
from .core.utilities import *
from .exceptions.exceptions import *


class Ardal(object):
    """ Class for handling Ardal-neighbourhood objects. """


    def __init__( self,
                  data_source : str,
                  use_roaring_if_sparse : bool=True,
                  density_threshold : float=0.02,
                  force_roaring : bool=False,
                  __ref : bool=False, 
                  file_format : str=None,
                  quiet : bool = False ):
        """ Ardal constructor
        """
        from . import _ardal

        super(Ardal, self).__init__()

        self.__hybrid_matrix = None
        self.__headers = None
        self.__ref = __ref

        parser = ArdalParser(data_source, file_format)

        ## do some parameter fiddling to enforce the generation of a roaring matrix
        if force_roaring:
            use_roaring_if_sparse = True
            density_threshold = 1.0

        if parser.matrix is not None:
            self.__hybrid_matrix = _ardal.HybridMatrix(parser.matrix,
                                                       use_roaring_if_sparse=use_roaring_if_sparse,
                                                       density_threshold=density_threshold)
            self.__headers = parser.headers
        else:
            ## raise an error if parsing fails to prevent unexpected behaviour down the line
            raise MatrixParseError(f"Failed to parse data from: {data_source}") 
        
        self.roaring = self.__hybrid_matrix.roaringEnabled()

        self.io = ArdalIO(self.__headers, self.__hybrid_matrix, self.roaring)
        self.get = ArdalGet(self.__headers, self.__hybrid_matrix, self.roaring)
        self.allele = ArdalAllele(self.__headers, self.__hybrid_matrix, self.roaring)
        self.compare = ArdalCompare(self.__headers, self.__hybrid_matrix, self.roaring)
        self.stats = ArdalStats(self.__headers, self.__hybrid_matrix, self.roaring)

        if not quiet:
            self.get.matrixStats(print_stats=True)