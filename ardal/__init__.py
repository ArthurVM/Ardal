import os as _os

## put a clamp on numpy and blas threading
for var in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    _os.environ.setdefault(var, "1")
_os.environ.setdefault("OMP_PROC_BIND", "spread")
_os.environ.setdefault("OMP_PLACES", "cores")

from ._version import __version__
from .Ardal import Ardal
from .core.ArdalParser import ArdalParser

__all__ = ["Ardal", "__version__"]
