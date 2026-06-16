import numpy as np


def _access_coords(ardal_obj, coords):
    dense = ardal_obj.get.bit_matrix()
    return dense[coords[:, 0], coords[:, 1]]


def _guid_alleles(ardal_obj, guid):
    guid_map = ardal_obj.io.to_dict(backend="auto")
    return guid_map[guid]


def _total_mass(ardal_obj):
    return int(np.sum(ardal_obj.get.row_masses()))


def _gather_snps(ardal_obj, guids):
    guid_map = ardal_obj.io.to_dict(backend="auto")
    return [guid_map[guid] for guid in guids]


def test_access(benchmark, ardal_obj):
    """
    Benchmarks coordinate access via Ardal.get.bit_matrix.
    """
    dense = ardal_obj.get.bit_matrix()
    limit = min(5, dense.shape[0], dense.shape[1])
    coords = np.array([[i, i] for i in range(limit)], dtype=np.uint64)

    # Benchmark access by coordinates
    result = benchmark(_access_coords, ardal_obj, coords)
    assert isinstance(result, np.ndarray)
    assert result.shape == (limit,)

def test_accessGUID(benchmark, ardal_obj):
    """
    Benchmarks per-GUID allele lookup via Ardal.io.to_dict.
    """
    guid = ardal_obj.get.headers()["guids"][0]
    # Benchmark per-GUID allele lookup
    result = benchmark(_guid_alleles, ardal_obj, guid)
    assert isinstance(result, list)

def test_getMatrix(benchmark, ardal_obj):
    """
    Benchmarks the get.bit_matrix method of the Ardal class.
    """
    # Benchmark get.bit_matrix
    result = benchmark(ardal_obj.get.bit_matrix)
    assert isinstance(result, np.ndarray)

def test_getMass(benchmark, ardal_obj):
    """
    Benchmarks total mass derived from get.row_masses.
    """
    # Benchmark total mass
    result = benchmark(_total_mass, ardal_obj)
    assert isinstance(result, int)

def test_hamming(benchmark, ardal_obj):
    """
    Benchmarks the hamming metric via Ardal.distance.pairwise.
    """
    ## Benchmark the hamming metric
    result = benchmark(ardal_obj.distance.pairwise, metric="hamming")
    assert isinstance(result, np.ndarray)

def test_jaccard(benchmark, ardal_obj):
    """
    Benchmarks the jaccard metric via Ardal.distance.pairwise.
    """
    ## Benchmark the jaccard metric
    result = benchmark(ardal_obj.distance.pairwise, metric="jaccard")
    assert isinstance(result, np.ndarray)

def test_neighbourhood(benchmark, ardal_obj):
    """
    Benchmarks the neighbourhood method of the ArdalDistance class.
    """
    guid = ardal_obj.get.headers()["guids"][0]
    ## Benchmark the neighbourhood method (scalar)
    result = benchmark(ardal_obj.distance.neighbourhood, guid, 2, use_simd=False)
    assert isinstance(result, dict)

def test_neighbourhoodSIMD(benchmark, ardal_obj):
    """
    Benchmarks the neighbourhood method (SIMD) of the ArdalDistance class.
    """
    guid = ardal_obj.get.headers()["guids"][0]
    ## Benchmark the neighbourhood method (SIMD)
    result = benchmark(ardal_obj.distance.neighbourhood, guid, 2, use_simd=True)
    assert isinstance(result, dict)

def test_gatherSNPs(benchmark, ardal_obj):
    """
    Benchmarks gathering alleles per GUID via Ardal.io.to_dict.
    """
    guids = ardal_obj.get.headers()["guids"][:5]

    ## benchmark the gather alleles per guid
    result = benchmark(_gather_snps, ardal_obj, guids)
    assert isinstance(result, list)
