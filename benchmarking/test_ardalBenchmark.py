import pytest
import _ardal as ardal
import numpy as np

def test_access(benchmark, allele_matrix):
    """
    Benchmarks the access method of the AlleleMatrix class.
    """
    coords = np.array([[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]], dtype=np.uint64)

    # Benchmark the access method
    result = benchmark(allele_matrix.access, coords)
    assert isinstance(result, np.ndarray)
    assert result.shape == (5,)

def test_accessGUID(benchmark, allele_matrix):
    """
    Benchmarks the accessGUID method of the AlleleMatrix class.
    """
    # Benchmark the accessGUID method
    result = benchmark(allele_matrix.accessGUID, 0)
    assert isinstance(result, set)

def test_getMatrix(benchmark, allele_matrix):
    """
    Benchmarks the getMatrix method of the AlleleMatrix class.
    """
    # Benchmark the getMatrix method
    result = benchmark(allele_matrix.getMatrix)
    assert isinstance(result, np.ndarray)

def test_getNumRows(benchmark, allele_matrix):
    """
    Benchmarks the getNumRows method of the AlleleMatrix class.
    """
    # Benchmark the getNumRows method
    result = benchmark(allele_matrix.getNumRows)
    assert isinstance(result, int)

def test_getNumCols(benchmark, allele_matrix):
    """
    Benchmarks the getNumCols method of the AlleleMatrix class.
    """
    # Benchmark the getNumCols method
    result = benchmark(allele_matrix.getNumCols)
    assert isinstance(result, int)

def test_getMass(benchmark, allele_matrix):
    """
    Benchmarks the getMass method of the AlleleMatrix class.
    """
    # Benchmark the getMass method
    result = benchmark(allele_matrix.getMass, 0)
    assert isinstance(result, int)

def test_hamming(benchmark, distance, distance_cache):
    """
    Benchmarks the hamming method of the Distance class.
    """
    ## clear cache
    distance_cache.clear()
    ## Benchmark the hamming method
    result = benchmark(distance.hamming, nocache=1)
    assert isinstance(result, np.ndarray)

def test_jaccard(benchmark, distance, distance_cache):
    """
    Benchmarks the jaccard method of the Distance class.
    """
    ## clear cache
    distance_cache.clear()
    ## Benchmark the jaccard method
    result = benchmark(distance.jaccard)
    assert isinstance(result, np.ndarray)

def test_neighbourhood(benchmark, neighbourhood, distance_cache):
    """
    Benchmarks the neighbourhood method of the Neighbourhood class.
    """
    ## clear cache
    distance_cache.clear()
    ## Benchmark the neighbourhood method
    result = benchmark(neighbourhood.neighbourhood, 0, 2, nocache=1)
    assert isinstance(result, np.ndarray)

def test_neighbourhoodSIMD(benchmark, neighbourhood, distance_cache):
    """
    Benchmarks the neighbourhoodSIMD method of the Neighbourhood class.
    """
    ## clear cache
    distance_cache.clear()
    ## Benchmark the neighbourhoodSIMD method
    result = benchmark(neighbourhood.neighbourhoodSIMD, 0, 2, nocache=1)
    assert isinstance(result, list)

def test_gatherSNPs(benchmark, allele_matrix):
    """
    Benchmarks the gatherSNPs method of the AlleleMatrix class.
    """
    guid_indices = np.array([0, 1, 2, 3, 4], dtype=np.int32)

    # Benchmark the gatherSNPs method
    result = benchmark(allele_matrix.gatherSNPs, guid_indices)
    assert isinstance(result, list)

def test_cache_get(benchmark, distance_cache):
    """
    Benchmarks the get method of the DistanceCache class.
    """
    distance_cache.put(0, 1, 5)
    # Benchmark the get method
    result = benchmark(distance_cache.get, 0, 1)
    assert isinstance(result, int)

def test_cache_put(benchmark, distance_cache):
    """
    Benchmarks the put method of the DistanceCache class.
    """
    # Benchmark the put method
    benchmark(distance_cache.put, 0, 1, 5)

def test_cache_clear(benchmark, distance_cache):
    """
    Benchmarks the clear method of the DistanceCache class.
    """
    distance_cache.put(0, 1, 5)
    # Benchmark the clear method
    benchmark(distance_cache.clear)
