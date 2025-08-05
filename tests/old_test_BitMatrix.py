import pytest
import numpy as np

## the tests here are for the public methods of the c++ bitmatrix class

@pytest.fixture
def bit_matrix(ardal_object):
    """Extracts the BitMatrix object from the Ardal object."""
    ## accessing private member for testing purposes
    return ardal_object._Ardal__bit_matrix


def test_get_matrix(bit_matrix, test_data):
    """Tests the getMatrix method of the BitMatrix class."""
    matrix_data, _ = test_data
    result = bit_matrix.getMatrix()
    assert isinstance(result, np.ndarray)
    np.testing.assert_array_equal(result, matrix_data)


def test_get_row_masses(bit_matrix):
    """Tests the getRowMasses method of the BitMatrix class."""
    result = bit_matrix.getRowMasses()
    assert isinstance(result, np.ndarray)
    assert len(result) == 10  ## from conftest.py, 10 guids
    ## from conftest.py, guid1 has snp1, snp5, snp6, snp7, snp8, snp10 -> 6
    assert result[0] == 6
    ## guid10 has snp4, snp5, snp7, snp9 -> 4
    assert result[9] == 4


def test_get_column_masses(bit_matrix):
    """Tests the getColumnMasses method of the BitMatrix class."""
    result = bit_matrix.getColumnMasses()
    assert isinstance(result, np.ndarray)
    assert len(result) == 10 ## from conftest.py, 10 snps
    ## from conftest.py, snp1 is in guid1, guid2 -> 2
    assert result[0] == 2
    ## snp5 is in all 10
    assert result[4] == 10


def test_get_set_bit_indices(bit_matrix):
    """Tests the getSetBitIndices method of the BitMatrix class."""
    ## for guid1 (index 0), snps are 1, 5, 6, 7, 8, 10. indices are 0, 4, 5, 6, 7, 9
    result = bit_matrix.getSetBitIndices(0)
    assert isinstance(result, np.ndarray)
    np.testing.assert_array_equal(np.sort(result), [0, 4, 5, 6, 7, 9])

    with pytest.raises(IndexError):
        bit_matrix.getSetBitIndices(100)


def test_hamming(bit_matrix):
    """
    Tests the hamming method of the BitMatrix class.
    """
    result = bit_matrix.hamming()
    assert isinstance(result, np.ndarray)
    ## total pairs for 10 rows is 10 * 9 / 2 = 45
    assert len(result) == 45


def test_neighbourhood(bit_matrix):
    """
    Tests the neighbourhood method of the BitMatrix class.
    """
    ## guid1 vs guid2: hamming distance is 0
    ## guid1 vs guid3: hamming distance is 2
    ## guid1 vs guid6: hamming distance is 4
    
    ## test with simd
    result_simd = bit_matrix.neighbourhood(row_idx=0, epsilon=2, use_simd=True)
    assert isinstance(result_simd, list)
    result_simd_dict = dict(result_simd)
    assert result_simd_dict == {1: 0, 2: 2, 3: 2, 4: 2} ## guid2, guid3, guid4, guid5

    ## test without simd
    result_scalar = bit_matrix.neighbourhood(row_idx=0, epsilon=2, use_simd=False)
    assert isinstance(result_scalar, list)
    result_scalar_dict = dict(result_scalar)
    assert result_scalar_dict == {1: 0, 2: 2, 3: 2, 4: 2}

    ## test invalid row coord
    with pytest.raises(RuntimeError):
        bit_matrix.neighbourhood(100, 2)

    ## test invalid epsilon
    with pytest.raises(RuntimeError):
        bit_matrix.neighbourhood(0, -1)


def test_inner_product_neighbourhood(bit_matrix):
    """Tests the innerProductNeighbourhood method of the BitMatrix class."""
    ## guid1 has 6 set bits
    ## guid2 has 6 set bits. ip is 6
    ## guid8 has 5 set bits. ip is 4
    
    ## test with simd
    result_simd = bit_matrix.innerProductNeighbourhood(row_idx=0, ip_epsilon=5, use_simd=True)
    assert isinstance(result_simd, list)
    result_simd_dict = dict(result_simd)
    assert 1 in result_simd_dict ## guid2
    assert result_simd_dict[1] == 6

    ## test without simd
    result_scalar = bit_matrix.innerProductNeighbourhood(row_idx=0, ip_epsilon=5, use_simd=False)
    assert isinstance(result_scalar, list)
    result_scalar_dict = dict(result_scalar)
    assert 1 in result_scalar_dict ## guid2
    assert result_scalar_dict[1] == 6

    ## test invalid row coord
    with pytest.raises(RuntimeError):
        bit_matrix.innerProductNeighbourhood(100, 2)

    ## test invalid epsilon
    with pytest.raises(RuntimeError):
        bit_matrix.innerProductNeighbourhood(0, -1)
