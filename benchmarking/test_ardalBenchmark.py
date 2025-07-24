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

def test_getMass(benchmark, allele_matrix):
    """
    Benchmarks the getMass method of the AlleleMatrix class.
    """
    # Benchmark the getMass method
    result = benchmark(allele_matrix.getMass)
    assert isinstance(result, int)

def test_hamming(benchmark, allele_matrix):
    """
    Benchmarks the hamming method of the Distance class.
    """
    ## clear cache
    allele_matrix.flushCache()
    ## Benchmark the hamming method
    result = benchmark(allele_matrix.hamming)
    assert isinstance(result, np.ndarray)

def test_jaccard(benchmark, allele_matrix):
    """
    Benchmarks the jaccard method of the Distance class.
    """
    ## clear cache
    allele_matrix.flushCache()
    ## Benchmark the jaccard method
    result = benchmark(allele_matrix.jaccard)
    assert isinstance(result, np.ndarray)

def test_neighbourhood(benchmark, allele_matrix):
    """
    Benchmarks the neighbourhood method of the Neighbourhood class.
    """
    ## clear cache
    allele_matrix.flushCache()
    ## Benchmark the neighbourhood method
    result = benchmark(allele_matrix.neighbourhood, 0, 2)
    assert isinstance(result, np.ndarray)

def test_neighbourhoodSIMD(benchmark, allele_matrix):
    """
    Benchmarks the neighbourhoodSIMD method of the Neighbourhood class.
    """
    ## clear cache
    allele_matrix.flushCache()
    ## Benchmark the neighbourhoodSIMD method
    result = benchmark(allele_matrix.neighbourhoodSIMD, 0, 2)
    assert isinstance(result, list)

def test_gatherSNPs(benchmark, allele_matrix):
    """
    Benchmarks the gatherSNPs method of the AlleleMatrix class.
    """
    guid_indices = np.array([0, 1, 2, 3, 4], dtype=np.int32)

    ## benchmark the gatherSNPs method
    result = benchmark(allele_matrix.gatherSNPs, guid_indices)
    assert isinstance(result, list)

def test_cache_clear(benchmark, allele_matrix):
    """
    Benchmarks the clear method of the DistanceCache class.
    """
    ## benchmark the clear method
    benchmark(allele_matrix.flushCache())
