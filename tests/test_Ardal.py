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
        
    with pytest.raises(UnsupportedFormatError, match="Unsupported file extension: .txt"):
        Ardal(data_source="./data/wrong_extension.txt", quiet_init=True)
        

def test_ardal_init_fail__malformed_list(test_data_mem):
    with pytest.raises(MalformedInputError, match="Input list must contain two elements: matrix and headers."):
        Ardal(data_source=["too", "many", "elements"], quiet_init=True)

    with pytest.raises(FileNotFoundError, match="One or more file paths do not exist: data/not_a_file.npy, data/sim_headers.json"):
        Ardal(data_source=["./data/not_a_file.npy", "./data/sim_headers.json"], quiet_init=True)
        
    with pytest.raises(MalformedInputError, match="If list input, it must be \\[headers::dict, matrix::np.ndarray\\] or \\[headers.json, matrix.bin\\] or dense pairs."):
        Ardal(data_source=[1, 2], quiet_init=True)

    with pytest.raises(MalformedInputError, match="If list input, it must be \\[headers::dict, matrix::np.ndarray\\] or \\[headers.json, matrix.bin\\] or dense pairs."):
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
    entropy = ardal.stats.allele_entropy()
    assert isinstance(entropy, dict)
    

def test_distance_pairwise(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    df = ardal.distance.pairwise()
    assert hasattr(df, "shape")


def test_distance_pairwise_jaccard(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, roaring=True, quiet_init=True)
    if not hasattr(ardal._hybrid_matrix, "jaccard") or not ardal.roaring:
        pytest.skip("Jaccard distance backend not available in this build.")
    dist = ardal.distance.pairwise(metric="jaccard", backend="roaring", threads=1)
    dense = test_data_mem[0].astype(np.uint8, copy=False)
    dense_bool = dense.astype(bool, copy=False)
    intersection = np.logical_and(
        dense_bool[:, None, :], dense_bool[None, :, :]
    ).sum(axis=2, dtype=np.int32)
    union = np.logical_or(
        dense_bool[:, None, :], dense_bool[None, :, :]
    ).sum(axis=2, dtype=np.int32)
    with np.errstate(divide="ignore", invalid="ignore"):
        similarity = np.divide(
            intersection,
            union,
            out=np.ones_like(intersection, dtype=np.float32),
            where=union > 0,
        )
    expected_full = similarity.astype(np.float32, copy=False)
    iu, ju = np.triu_indices(expected_full.shape[0], 1)
    expected_condensed = expected_full[iu, ju]
    np.testing.assert_allclose(dist, expected_condensed, rtol=1e-6, atol=1e-7)


def test_distance_pairwise_inner_product(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    dist = ardal.distance.pairwise(metric="inner_product", backend="bit", threads=1)
    dense = test_data_mem[0].astype(np.uint32, copy=False)
    gram = dense @ dense.T
    iu, ju = np.triu_indices(gram.shape[0], 1)
    expected_condensed = gram[iu, ju].astype(np.uint32, copy=False)
    np.testing.assert_array_equal(dist, expected_condensed)


def test_get_matrix_stats(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    stats = ardal.get.matrix_stats()
    assert isinstance(stats, dict)


def test_cosine_distance_matrix(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    if not hasattr(ardal._hybrid_matrix, "cosineDistance"):
        pytest.skip("cosineDistance backend not available in this build.")
    mat = ardal._hybrid_matrix.cosineDistance(use_simd=True, threads=1, backend="bit")
    assert mat.ndim == 1
    dense = test_data_mem[0].astype(np.float64, copy=False)
    dot = dense @ dense.T
    norms = np.linalg.norm(dense, axis=1)
    denom = norms[:, None] * norms[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = np.divide(dot, denom, out=np.zeros_like(dot, dtype=np.float64), where=denom > 0.0)
    zero_mask = norms == 0.0
    only_one_zero = np.logical_xor(zero_mask[:, None], zero_mask[None, :])
    cosine = np.where(only_one_zero, 0.0, cosine)
    expected = 1.0 - cosine
    expected = np.where(np.logical_and(zero_mask[:, None], zero_mask[None, :]), 0.0, expected)
    iu, ju = np.triu_indices(expected.shape[0], 1)
    expected_condensed = expected[iu, ju]
    np.testing.assert_allclose(mat, expected_condensed, rtol=1e-7, atol=1e-9)


def test_distance_pairwise_cosine(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    if not hasattr(ardal._hybrid_matrix, "cosineDistance"):
        pytest.skip("cosineDistance backend not available in this build.")
    dist = ardal.distance.pairwise(metric="cosine", threads=1)
    dense = test_data_mem[0].astype(np.float64, copy=False)
    dot = dense @ dense.T
    norms = np.linalg.norm(dense, axis=1)
    denom = norms[:, None] * norms[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = np.divide(dot, denom, out=np.zeros_like(dot, dtype=np.float64), where=denom > 0.0)
    zero_mask = norms == 0.0
    only_one_zero = np.logical_xor(zero_mask[:, None], zero_mask[None, :])
    cosine = np.where(only_one_zero, 0.0, cosine)
    expected = 1.0 - cosine
    expected = np.where(np.logical_and(zero_mask[:, None], zero_mask[None, :]), 0.0, expected)
    iu, ju = np.triu_indices(expected.shape[0], 1)
    expected_condensed = expected[iu, ju]
    np.testing.assert_allclose(dist, expected_condensed, rtol=1e-7, atol=1e-9)


def test_distance_pairwise_snv(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    if not hasattr(ardal._hybrid_matrix, "snvHamming"):
        pytest.skip("SNV distance backend not available in this build.")
    hamming = ardal.distance.pairwise(metric="hamming", threads=1)
    snv = ardal.distance.pairwise(
        metric="snv",
        threads=1,
        allele_id_format="{chr}.{start}.{ref}.{alt}",
    )
    np.testing.assert_array_equal(snv, hamming)


def test_knn_metric_selection(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    if not hasattr(ardal._hybrid_matrix, "knn_hamming"):
        pytest.skip("metric-aware knn backend not available in this build.")
    metrics = ["hamming", "inner_product", "jaccard", "cosine"]
    if hasattr(ardal._hybrid_matrix, "knnSnv"):
        metrics.append("snv")
    for metric in metrics:
        kwargs = {"metric": metric, "backend": "bit"}
        if metric == "snv":
            kwargs["allele_id_format"] = "{chr}.{start}.{ref}.{alt}"
        neighbours = ardal.distance.knn("GUID1", 3, **kwargs)
        assert len(neighbours) == 3


def test_knn_metric_value_types(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    if not hasattr(ardal._hybrid_matrix, "knn_cosine"):
        pytest.skip("cosine knn backend not available in this build.")

    cosine_neighbours = ardal.distance.knn("GUID1", 3, metric="cosine", backend="bit")
    assert cosine_neighbours
    assert all(isinstance(val, float) for val in cosine_neighbours.values())

    if hasattr(ardal._hybrid_matrix, "knnSnv"):
        snv_neighbours = ardal.distance.knn(
            "GUID1",
            3,
            metric="snv",
            backend="bit",
            allele_id_format="{chr}.{start}.{ref}.{alt}",
        )
        assert snv_neighbours
        assert all(isinstance(val, int) for val in snv_neighbours.values())


def test_neighbourhood_metric_selection(test_data_mem):
    ardal = Ardal(data_source=test_data_mem, quiet_init=True)
    ham = ardal.distance.neighbourhood("GUID1", 3, metric="hamming", backend="bit")
    assert isinstance(ham, dict) and len(ham) > 0
    ip = ardal.distance.neighbourhood("GUID1", 1, metric="inner_product", backend="bit")
    assert isinstance(ip, dict)
    if hasattr(ardal._hybrid_matrix, "snvHamming"):
        snv = ardal.distance.neighbourhood(
            "GUID1",
            3,
            metric="snv",
            backend="bit",
            allele_id_format="{chr}.{start}.{ref}.{alt}",
        )
        assert isinstance(snv, dict)

# def test_with_coords_bed(matrix_npy, coords_bed_file):
#     ardal = Ardal(data_source=matrix_npy, coords_bed=coords_bed_file, quiet_init=True)
#     assert ardal._headerUtils is not None
