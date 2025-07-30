# ardal/Exceptions/exceptions.py

class ArdalError(Exception):
    """Base class for all Ardal-specific exceptions."""
    pass


## I/O and file errors
class FileFormatError(ArdalError):
    """Raised when a file format is unsupported or malformed."""
    pass

class MissingFileError(ArdalError):
    """Raised when a required file is not found."""
    pass

class HeaderMismatchError(ArdalError):
    """Raised when header data does not align with matrix dimensions."""
    pass

class MatrixParseError(ArdalError):
    """Raised by Ardal base when matrix parsing fails"""
    pass

class MatrixWriteError(ArdalError):
    """Raised by Ardal base when matrix writing fails"""
    pass


## matrix construction / validation
class MatrixBuildError(ArdalError):
    """Raised when matrix construction fails."""
    pass

class SparseMatrixError(ArdalError):
    """Raised when sparse matrix conversion fails or thresholds are violated."""
    pass

class IncompatibleMatrixError(ArdalError):
    """Raised when attempting an operation between incompatible matrices."""
    pass

class RoaringError(ArdalError):
    """Raised when attempting to use a roaring backend when it was not initialised during matrix parsing."""
    pass


## parsing errors
class VCFParseError(ArdalError):
    """Raised when a VCF file cannot be parsed correctly."""
    pass

class AlleleEncodingError(ArdalError):
    """Raised when allele encoding fails."""
    pass


## query errors
class InvalidAlleleQueryError(ArdalError):
    """Raised when a query requests a non-existent or invalid allele."""
    pass

class InvalidGUIDQueryError(ArdalError):
    """Raised when a query requests a non-existent or invalid GUID."""
    pass

class IndexOutOfBoundsError(ArdalError):
    """Raised when row/column index exceeds matrix dimensions."""
    pass

class EmptySelectionError(ArdalError):
    """Raised when an operation requires a non-empty input but receives none."""
    pass


## computational/statistical errors
class DivergenceComputationError(ArdalError):
    """Raised when KL/JS divergence or IG stats fail to compute."""
    pass

class SIMDNotSupportedError(ArdalError):
    """Raised when SIMD operations are requested but not available/compiled."""
    pass

class ThreadingError(ArdalError):
    """Raised when thread allocation or behavior fails."""
    pass


## general errors
class InvalidTypeError(ArdalError):
    """Raised when data of the wrong type is passed to a function."""
    pass

class ParameterError(ArdalError):
    """Raised when invalid parameters are passed to a function."""
    pass