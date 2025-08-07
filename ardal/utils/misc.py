""" utilities.py
This module provides generic utility functions for the Ardal framework.
"""
import importlib
from typing import Union


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