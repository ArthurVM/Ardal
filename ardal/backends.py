""" Module for importing the compiled C++ backends. 
This exists to avoid circular import issues. """

try:
    from . import ardal
except ImportError as e:
    raise ImportError(
        "Could not import the compiled Ardal C++ backends. "
        "Please ensure the project has been built correctly. "
        f"Original error: {e}"
    )