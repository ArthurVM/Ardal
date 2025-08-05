""" utilities.py
This module provides utility functions for working with GUIDs and alleles in the Ardal framework.
"""
import importlib
from typing import Union
from collections import defaultdict

from ..exceptions.exceptions import *


def require_package(package_name, import_as = None, attr: Union[str, None] = None, error_message=None):
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
        

def validate_backend_argument( func ):
    functools = require_package("functools", "functools")
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        signature = require_package("inspect", attr="signature")
        
        ALLOWED_BACKENDS = {"auto", "bit", "roaring"}

        sig = signature(func)
        bound_args = sig.bind(*args, **kwargs)
        bound_args.apply_defaults()

        backend = bound_args.arguments.get("backend", None)

        if backend is None:
            raise ValueError("Missing required argument: 'backend'")

        if isinstance(backend, str):
            backend_lower = backend.lower()
            if backend_lower not in ALLOWED_BACKENDS:
                raise InvalidBackendError(f"Invalid backend '{backend}'. Must be one of {ALLOWED_BACKENDS}")
            bound_args.arguments["backend"] = backend_lower
        else:
            raise InvalidBackendError(f"'backend' must be a string, not {type(backend).__name__}")

        return func(*bound_args.args, **bound_args.kwargs)

    return wrapper


def validate_Ardal_roaring_param( roaring_param : Union[str, bool] ) -> str:
    """ Checks the 'roaring' parameter passed to Ardal during object creation is acceptable.
    """
    ALLOWED_PARAMETERS = {"auto", "true", "false"}
    
    if not isinstance(roaring_param, str) and not isinstance(roaring_param, bool):
        raise ParameterError(f"'roaring' argument must be of type 'string' or 'bool' and one of {ALLOWED_PARAMETERS}.")
    
    lower_param = str(roaring_param).lower()
    if lower_param not in ALLOWED_PARAMETERS:
        raise ParameterError(f"'roaring' argument must be one of {ALLOWED_PARAMETERS}.")
    
    return lower_param