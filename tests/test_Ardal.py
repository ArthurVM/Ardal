import os
import numpy as np
import pytest
from ardal import Ardal
from ardal.exceptions.exceptions import MalformedInputError, UnsupportedFormatError, LoadMatrixError
    

## check initialisation with different inputs

## SUCCESSES
def test_ardal_init_success__mem(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet=True)
    assert ardal._hybrid_matrix is not None
    assert ardal.headerUtils is not None
    
    
def test_ardal_init_success__npy(test_data_matrix_npy):
    ardal = Ardal(data_source=test_data_matrix_npy, quiet=True)
    assert ardal._hybrid_matrix is not None
    assert ardal.headerUtils is not None


def test_ardal_init_success__npz(test_data_matrix_npz):
    ardal = Ardal(data_source=test_data_matrix_npz, quiet=True)
    assert ardal._hybrid_matrix is not None
    assert ardal.headerUtils is not None
    

def test_ardal_init_success__csv(test_data_matrix_csv):
    ardal = Ardal(data_source=test_data_matrix_csv, quiet=True)
    assert ardal._hybrid_matrix is not None
    assert ardal.headerUtils is not None
    

def test_ardal_init_success__makes_contiguous():
    ## check Ardal makes a non-contigous array contiguous
    tmp_matrix = np.ones((3, 4))
    non_c_matrix = tmp_matrix.T
    headers = {"guids" : [f"GUID{i}" for i in range(non_c_matrix.shape[0])], \
               "alleles" : [f"A.{i*100}.T" for i in range(non_c_matrix.shape[1])]}
    ardal = Ardal(data_source=[non_c_matrix, headers], quiet=True)
    assert ardal._hybrid_matrix is not None
    assert ardal.headerUtils is not None


## FAILS
def test_ardal_init_fail__None_input():
    with pytest.raises(MalformedInputError, match="Input data structure cannot be None."):
        Ardal(data_source=None, quiet=True)
        
    
def test_ardal_init_fail__malformed_string():
    with pytest.raises(FileNotFoundError, match=f"File does not exist: not_a_file.npy"):
        Ardal(data_source="not_a_file.npy", quiet=True)
        
    with pytest.raises(UnsupportedFormatError, match="Unsupported file format: must be csv, parquet, npy, or npz"):
        Ardal(data_source="./data/wrong_extension.txt", quiet=True)
        

def test_ardal_init_fail__malformed_list(test_data_mem):
    with pytest.raises(MalformedInputError, match="Input list must contain two elements: matrix and headers."):
        Ardal(data_source=["too", "many", "elements"], quiet=True)

    with pytest.raises(FileNotFoundError, match="One or more file paths do not exist: ./data/not_a_file.npy, ./data/sim_headers.json"):
        Ardal(data_source=["./data/not_a_file.npy", "./data/sim_headers.json"], quiet=True)
        
    with pytest.raises(MalformedInputError, match="If list input, must contain either \\[headers, np.matrix\\] or two file paths."):
        Ardal(data_source=[1, 2], quiet=True)

    with pytest.raises(MalformedInputError, match="If list input, must contain either \\[headers, np.matrix\\] or two file paths."):
        Ardal(data_source=[[], test_data_mem[1]], quiet=True)

def test_ardal_init_fail__malformed_data(test_data_mem):
    tmp_matrix = np.ones((3, 4))
    
    with pytest.raises(LoadMatrixError, match="Matrix must be 2-dimensional."):
        Ardal(data_source=[np.zeros((2,3,4)), test_data_mem[1]], quiet=True)
        
    with pytest.raises(LoadMatrixError, match="Headers must contain 'guids' and 'alleles' keys."):
        headers = {"wrong" : [f"GUID{i}" for i in range(tmp_matrix.shape[0])], \
                   "keys" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet=True)
    
    with pytest.raises(LoadMatrixError, match="GUIDs must be a list of strings."):
        headers = {"guids" : [i for i in range(tmp_matrix.shape[0])], \
                   "alleles" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet=True)
        
    with pytest.raises(LoadMatrixError, match="Alleles must be a list of strings."):
        headers = {"guids" : [f"GUID{i}" for i in range(tmp_matrix.shape[0])], \
                   "alleles" : [i for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet=True)
    
    with pytest.raises(LoadMatrixError, match=f"Mismatch: {tmp_matrix.shape[0]} matrix rows vs {tmp_matrix.shape[0]+1} GUIDs."):
        headers = {"guids" : [f"GUID{i}" for i in range(tmp_matrix.shape[0]+1)], \
                   "alleles" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet=True)
    
    with pytest.raises(LoadMatrixError, match=f"Mismatch: {tmp_matrix.shape[1]} matrix columns vs {tmp_matrix.shape[1]+1} alleles."):
        headers = {"guids" : [f"GUID{i}" for i in range(tmp_matrix.shape[0])], \
                   "alleles" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1]+1)]}
        Ardal(data_source=[tmp_matrix, headers], quiet=True)
        
    with pytest.raises(LoadMatrixError, match=f"GUIDs must be unique."):
        headers = {"guids" : [f"NON-UNIQUE-GUID" for i in range(tmp_matrix.shape[0])], \
                   "alleles" : [f"A.{i*100}.T" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet=True)
        
    with pytest.raises(LoadMatrixError, match=f"Alleles must be unique."):
        headers = {"guids" : [f"GUID{i}" for i in range(tmp_matrix.shape[0])], \
                   "alleles" : ["NON-UNIQUE-ALLELE" for i in range(tmp_matrix.shape[1])]}
        Ardal(data_source=[tmp_matrix, headers], quiet=True)
    

## check roaring is intialised correctly
def test_roaring_flag(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, roaring=True, quiet=True)
    assert ardal.roaring is True
    
    ardal = Ardal(data_source=test_data_mem, roaring=False, quiet=True)
    assert ardal.roaring is False
    
    ardal = Ardal(data_source=test_data_mem, roaring="auto", density_threshold=1.0, quiet=True)
    assert ardal.roaring is True
    
    ardal = Ardal(data_source=test_data_mem, roaring="auto", density_threshold=0.0, quiet=True)
    assert ardal.roaring is False
    

## general small tests
def test_stats_methods(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet=True)
    entropy = ardal.stats.alleleEntropy()
    assert isinstance(entropy, dict)
    

def test_compare_pairwise(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet=True)
    df = ardal.compare.pairwise()
    assert hasattr(df, "shape")
    

def test_get_matrix_stats(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet=True)
    stats = ardal.get.matrixStats()
    assert isinstance(stats, dict)
    

# def test_with_coords_bed(matrix_npy, coords_bed_file):
#     ardal = Ardal(data_source=matrix_npy, coords_bed=coords_bed_file, quiet=True)
#     assert ardal.headerUtils is not None