import functools
from inspect import signature
from typing import Union

from .exceptions import InvalidBackendError, ParameterError, AllelePatternError
from .validators import *

def check_backend_argument( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()

        val = bound_args.arguments.get("backend", None)
        bound_args.arguments["backend"] = validate_backend_argument(val)
        
        return func(*bound_args.args, **bound_args.kwargs)

    return wrapper


def check_density_threshold( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        
        val = bound_args.arguments.get("density_threshold", None)        
        validate_density_threshold(val)
        
        return func(*args, **kwargs)
    
    return wrapper


def check_generic_threshold( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        
        val = bound_args.arguments.get("threshold", None)        
        validate_generic_threshold(val)
        
        return func(*args, **kwargs)
    
    return wrapper
    

def check_Ardal_roaring_param( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """ Checks the 'roaring' parameter passed to Ardal during object creation is acceptable.
        """        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        
        val = bound_args.arguments.get("roaring", None)
        bound_args.arguments["roaring"] = validate_Ardal_roaring_param(val)
        
        return func(*bound_args.args, **bound_args.kwargs)
    
    return wrapper


def check_allele_id_format( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """ Basic checks on the 'allele_id_format' parameter.
        """        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        
        val = bound_args.arguments.get("allele_id_format", None)
        validate_allele_id_format(val)
            
        return func(*bound_args.args, **bound_args.kwargs)
    
    return wrapper


def check_thread_count( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """ Basic checks on the 'threads' parameter.
        """        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        
        val = bound_args.arguments.get("threads", None)
        validate_thread_count(val)
            
        return func(*bound_args.args, **bound_args.kwargs)
    
    return wrapper


def check_use_simd( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """ Basic checks on the 'use_simd' parameter.
        """        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        
        val = bound_args.arguments.get("use_simd", None)
        validate_use_simd(val)
            
        return func(*bound_args.args, **bound_args.kwargs)
    
    return wrapper


def check_bed_paths( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """ Checks made on input beds.
        """        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        
        if "intervals_bed" in bound_args.arguments:
            val = bound_args.arguments.get("intervals_bed", None)
            validate_type(val, Union[str, None], "intervals_bed")
        
        if "allele_coords_bed" in bound_args.arguments:
            val = bound_args.arguments.get("allele_coords_bed", None)
            validate_type(val, Union[str, None], "allele_coords_bed")
            
        return func(*bound_args.args, **bound_args.kwargs)
    
    return wrapper


def check_guids_list( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """ Basic checks on the 'guids' parameter.
        """        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        
        val = bound_args.arguments.get("guids", None)
        validate_guids(val)
            
        return func(*bound_args.args, **bound_args.kwargs)
    
    return wrapper


def check_alleles_list( func ):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """ Basic checks on the 'alleles' parameter.
        """        
        bound_args = signature(func).bind(*args, **kwargs)
        bound_args.apply_defaults()
        
        val = bound_args.arguments.get("alleles", None)
        validate_alleles(val)
            
        return func(*bound_args.args, **bound_args.kwargs)
    
    return wrapper