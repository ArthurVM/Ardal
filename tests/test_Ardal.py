import os
import numpy as np
import pytest
from ardal import Ardal
from ardal.utils.exceptions import MalformedInputError, UnsupportedFormatError, LoadMatrixError, InvalidBackendError
    

## check initialisation with different inputs

## SUCCESSES
def test_ardal_init_success__mem(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    assert ardal._hybrid_matrix is not None
    assert ardal._headerUtils is not None
    
    
def test_ardal_init_success__npy(test_data_matrix_npy):
    ardal = Ardal(data_source=test_data_matrix_npy, quiet_init=True)
    assert ardal._hybrid_matrix is not None
    assert ardal._headerUtils is not None


def test_ardal_init_success__npz(test_data_matrix_npz):
    ardal = Ardal(data_source=test_data_matrix_npz, quiet_init=True)
    assert ardal._hybrid_matrix is not None
    assert ardal._headerUtils is not None
    

def test_ardal_init_success__csv(test_data_matrix_csv):
    ardal = Ardal(data_source=test_data_matrix_csv, quiet_init=True)
    assert ardal._hybrid_matrix is not None
    assert ardal._headerUtils is not None
    

def test_ardal_init_success__makes_contiguous():
    ## check Ardal makes a non-contigous array contiguous
    tmp_matrix = np.ones((3, 4))
    non_c_matrix = tmp_matrix.T
    headers = {"guids" : [f"GUID{i}" for i in range(non_c_matrix.shape[0])], \
               "alleles" : [f"A.{i*100}.T" for i in range(non_c_matrix.shape[1])]}
    ardal = Ardal(data_source=[non_c_matrix, headers], quiet_init=True)
    assert ardal._hybrid_matrix is not None
    assert ardal._headerUtils is not None


## FAILS
def test_ardal_init_fail__None_input():
    with pytest.raises(MalformedInputError, match="Input data structure cannot be None."):
        Ardal(data_source=None, quiet_init=True)
        
    
def test_ardal_init_fail__malformed_string():
    with pytest.raises(FileNotFoundError, match=f"File does not exist: not_a_file.npy"):
        Ardal(data_source="not_a_file.npy", quiet_init=True)
        
    with pytest.raises(UnsupportedFormatError, match="Unsupported file format: must be csv, parquet, npy, or npz"):
        Ardal(data_source="./data/wrong_extension.txt", quiet_init=True)
        

def test_ardal_init_fail__malformed_list(test_data_mem):
    with pytest.raises(MalformedInputError, match="Input list must contain two elements: matrix and headers."):
        Ardal(data_source=["too", "many", "elements"], quiet_init=True)

    with pytest.raises(FileNotFoundError, match="One or more file paths do not exist: ./data/not_a_file.npy, ./data/sim_headers.json"):
        Ardal(data_source=["./data/not_a_file.npy", "./data/sim_headers.json"], quiet_init=True)
        
    with pytest.raises(MalformedInputError, match="If list input, must contain either \\[headers::json, matrix::np.matrix\\] or two file paths."):
        Ardal(data_source=[1, 2], quiet_init=True)

    with pytest.raises(MalformedInputError, match="If list input, must contain either \\[headers::json, matrix::np.matrix\\] or two file paths."):
        Ardal(data_source=[[], test_data_mem[1]], quiet_init=True)

def test_ardal_init_fail__malformed_data(test_data_mem):
    tmp_matrix = np.ones((3, 4))
    
    with pytest.raises(LoadMatrixError, match="Matrix must be 2-dimensional."):
        Ardal(data_source=[np.zeros((2,3,4)), test_data_mem[1]], quiet_init=True)
        
    with pytest.raises(LoadMatrixError, match="Headers must contain 'guids' and 'alleles' keys."):
        headers = {"wrong" : [f"GUID{i}" for i in range(tmp_matrix.shape[0])], \
                   "keys" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet_init=True)
    
    with pytest.raises(LoadMatrixError, match="GUIDs must be a list of strings."):
        headers = {"guids" : [i for i in range(tmp_matrix.shape[0])], \
                   "alleles" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet_init=True)
        
    with pytest.raises(LoadMatrixError, match="Alleles must be a list of strings."):
        headers = {"guids" : [f"GUID{i}" for i in range(tmp_matrix.shape[0])], \
                   "alleles" : [i for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet_init=True)
    
    with pytest.raises(LoadMatrixError, match=f"Mismatch: {tmp_matrix.shape[0]} matrix rows vs {tmp_matrix.shape[0]+1} GUIDs."):
        headers = {"guids" : [f"GUID{i}" for i in range(tmp_matrix.shape[0]+1)], \
                   "alleles" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet_init=True)
    
    with pytest.raises(LoadMatrixError, match=f"Mismatch: {tmp_matrix.shape[1]} matrix columns vs {tmp_matrix.shape[1]+1} alleles."):
        headers = {"guids" : [f"GUID{i}" for i in range(tmp_matrix.shape[0])], \
                   "alleles" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1]+1)]}
        Ardal(data_source=[tmp_matrix, headers], quiet_init=True)
        
    with pytest.raises(LoadMatrixError, match=f"GUIDs must be unique."):
        headers = {"guids" : [f"NON-UNIQUE-GUID" for i in range(tmp_matrix.shape[0])], \
                   "alleles" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet_init=True)
        
    with pytest.raises(LoadMatrixError, match=f"Alleles must be unique."):
        headers = {"guids" : [f"GUID{i}" for i in range(tmp_matrix.shape[0])], \
                   "alleles" : ["NON-UNIQUE-ALLELE" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet_init=True)
    

## check roaring flags are handled correctly
def test_roaring_flag(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, roaring=True, quiet_init=True)
    assert ardal.roaring is True
    
    ardal = Ardal(data_source=test_data_mem, roaring=False, quiet_init=True)
    assert ardal.roaring is False
    
    ardal = Ardal(data_source=test_data_mem, roaring="auto", density_threshold=1.0, quiet_init=True)
    assert ardal.roaring is True
    
    ardal = Ardal(data_source=test_data_mem, roaring="auto", density_threshold=0.0, quiet_init=True)
    assert ardal.roaring is False
    
    with pytest.raises(InvalidBackendError, match=f"'roaring' argument must be of type 'string' or 'bool' and one of"):
        ardal = Ardal(data_source=test_data_mem, roaring=None, quiet_init=True)
        
    with pytest.raises(InvalidBackendError, match=f"'roaring' argument must be of type 'string' or 'bool' and one of"):
        ardal = Ardal(data_source=test_data_mem, roaring=1, quiet_init=True)
        
    with pytest.raises(InvalidBackendError, match=f"'roaring' argument must be of type 'string' or 'bool' and one of"):
        ardal = Ardal(data_source=test_data_mem, roaring=[], quiet_init=True)
        
    with pytest.raises(InvalidBackendError, match=f"'roaring' argument must be one of"):
        ardal = Ardal(data_source=test_data_mem, roaring="invalid", quiet_init=True)

    
    

## general small tests
def test_stats_methods(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    entropy = ardal.stats.alleleEntropy()
    assert isinstance(entropy, dict)
    

def test_distance_pairwise(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    df = ardal.distance.pairwise()
    assert hasattr(df, "shape")
    

def test_get_matrix_stats(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    stats = ardal.get.matrix_stats()
    assert isinstance(stats, dict)
    

# def test_with_coords_bed(matrix_npy, coords_bed_file):
#     ardal = Ardal(data_source=matrix_npy, coords_bed=coords_bed_file, quiet_init=True)
#     assert ardal._headerUtils is not None