## ardal/utils/validators.py
from typing import Union, get_origin, get_args
import types as _types

from .exceptions import InvalidBackendError, ParameterError, AllelePatternError, IntervalError


def validate_backend_argument( backend ) -> str:
    ALLOWED_BACKENDS = {"auto", "bit", "roaring"}

    if backend is None:
        raise ValueError("Missing required argument: 'backend'")

    if isinstance(backend, str):
        backend_lower = backend.lower()
        if backend_lower not in ALLOWED_BACKENDS:
            raise InvalidBackendError(f"Invalid backend '{backend}'. Must be one of {ALLOWED_BACKENDS}.")
        return backend_lower
    else:
        raise InvalidBackendError(f"'backend' must be a string, not {type(backend).__name__}.")
    

def validate_density_threshold( density_threshold ) -> None:
    ## handles 0 or 1 ints
    if isinstance(density_threshold, int):
        density_threshold = float(density_threshold)
    
    validate_type(density_threshold, float, name="density_threshold")
    
    if density_threshold < 0.0 or density_threshold > 1.0:
        raise ValueError("Density threshold must be a float between 0-1.")
    
    
def validate_generic_threshold( threshold ) -> None:
    ## handles 0 or 1 ints
    if isinstance(threshold, int):
        threshold = float(threshold)
    
    validate_type(threshold, float, name="threshold")
    
    if threshold < 0.0 or threshold > 1.0:
        raise ValueError("threshold must be a float between 0-1.")
    
    
def validate_Ardal_roaring_param( roaring_param ) -> str:
    ALLOWED_PARAMETERS = {"auto", "true", "false"}
    
    ## function specific type handling
    if not isinstance(roaring_param, str) and not isinstance(roaring_param, bool):
        raise InvalidBackendError(f"'roaring' argument must be of type 'string' or 'bool' and one of {ALLOWED_PARAMETERS}.")
        
    roaring_param_lower = str(roaring_param).lower()
    if roaring_param_lower not in ALLOWED_PARAMETERS:
        raise InvalidBackendError(f"'roaring' argument must be one of {ALLOWED_PARAMETERS}.")
    
    return roaring_param_lower


def validate_allele_id_format( allele_id_format ) -> None:
    if allele_id_format != None:
        validate_type(allele_id_format, str, name="allele_id_format")
        
        ## check placeholder patterns are contained within braces
        c1 = allele_id_format.count("{")
        c2 = allele_id_format.count("}")
        
        if c1 == 0 or c2 == 0 or c1 != c2:
            raise AllelePatternError("'allele_id_format' argument malformed, please ensure placeholders are contained within curly braces (e.g.{chr}.{start}.{alt}).")


def validate_thread_count( threads ) -> None:
    validate_type(threads, int, name="threads")
    if threads < 1:
        raise ParameterError("threads must be at least 1.")
    

def validate_use_simd( use_simd ) -> None:
    validate_type(use_simd, bool, name="use_simd")


def validate_coords_bed( coords_bed ) -> None:
    validate_type(coords_bed, str, name="coords_bed")
    

def validate_guids( guids ) -> None:
    validate_type(guids, Union[list, None], "guids")
    
    if guids:
        for i in guids:
            validate_type(i, Union[str, None], f"guids[{i}]")
        

def validate_alleles( alleles ) -> None:
    validate_type(alleles, Union[list, None], "alleles")
    
    if alleles:
        for i in alleles:
            validate_type(i, Union[str, None], f"alleles[{i}]")


def validate_type( value,
                   expected_type,
                   name: str = "value"
                   ) -> None:
    """
    Validates that a given value is of the expected type.

    Args:
        value: The value to check.
        expected_type: The expected type or tuple of types.
        name (str): The name of the variable (for error messages).

    Raises:
        TypeError: If the value is not of the expected type.
    """
    origin = get_origin(expected_type)

    if origin in (Union, _types.UnionType):              # Union[X, Y] or X | Y
        types_ = tuple(get_args(expected_type))
        if not all(isinstance(t, type) for t in types_):
            raise TypeError("Union must contain only runtime-checkable classes.")
    elif isinstance(expected_type, type):                 # single class
        types_ = (expected_type,)
    else:
        raise TypeError("expected_type must be a class or a Union[...] / X|Y.")

    if not isinstance(value, types_):
        names = " | ".join(t.__name__ for t in types_)
        raise TypeError(f"{name} must be of type {names}, not {type(value).__name__}.")