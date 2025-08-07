import pytest
import numpy as np

from ardal.utils.exceptions import InvalidGUIDQueryError, ParameterError, RoaringError


def test_subset(get_component):
    """
    Tests the subset method of the ArdalGet class.
    """
    # Test subset by GUIDs
    guid_subset = ["GUID1", "GUID10"]
    matrix, headers = get_component.subset(guids=guid_subset, data_only=True)
    assert headers["guids"] == guid_subset
    assert matrix.shape == (2, 10)

    # Test subset by Alleles
    allele_subset = ["SNP1", "SNP10"]
    matrix, headers = get_component.subset(alleles=allele_subset, data_only=True)
    assert headers["alleles"] == allele_subset
    assert matrix.shape == (10, 2)

    # Test subset by both
    matrix, headers = get_component.subset(guids=guid_subset, alleles=allele_subset, data_only=True)
    assert headers["guids"] == guid_subset
    assert headers["alleles"] == allele_subset
    assert matrix.shape == (2, 2)
    # GUID1 has SNP1 and SNP10, GUID10 does not have SNP1 or SNP10
    expected_matrix = np.array([[1, 1], [0, 0]], dtype=np.uint8)
    np.testing.assert_array_equal(matrix, expected_matrix)

    # Test invalid GUID
    with pytest.raises(InvalidGUIDQueryError, match=f"guids "):  ## doesnt seem to like the end of the warning string
        get_component.subset(guids=["INVALID_GUID"], data_only=True)

    # Test empty lists
    with pytest.raises(ParameterError, match="guids and alleles cannot both be empty."):
        get_component.subset(guids=[], alleles=[], data_only=True)


def test_matrix_stats(get_component, capsys):
    """
    Tests the matrixStats method of the ArdalGet class.
    """
    stats = get_component.matrixStats(print_stats=False)
    assert isinstance(stats, dict)
    assert "Number of GUIDs" in stats
    assert "Matrix Density" in stats
    assert stats["Number of GUIDs"] == 10
    assert stats["Number of Alleles"] == 10

    # Test printing
    get_component.matrixStats(print_stats=True)
    captured = capsys.readouterr()
    assert "Ardal Matrix Statistics" in captured.out


def test_density(get_component):
    """
    Tests the density method of the ArdalGet class.
    """
    # 56 set bits out of 100 total elements in test_data
    assert get_component.density() == pytest.approx(0.56)


def test_bit_matrix(get_component, test_data_mem):
    """
    Tests the bitMatrix method of the ArdalGet class.
    """
    matrix, _ = test_data_mem
    np.testing.assert_array_equal(get_component.bitMatrix(), matrix)


def test_roaring_matrix(get_component, get_component_no_roaring):
    """
    Tests the roaringMatrix method of the ArdalGet class.
    """
    ## test decoded output
    roaring_dict = get_component.roaringMatrix(decode=True)
    assert isinstance(roaring_dict, dict)
    assert "GUID1" in roaring_dict
    assert set(roaring_dict["GUID1"]) == {"SNP1", "SNP5", "SNP6", "SNP7", "SNP8", "SNP10"}

    ## test raw output
    raw_roaring = get_component.roaringMatrix(decode=False)
    assert isinstance(raw_roaring, list)
    assert len(raw_roaring) == 10
    assert isinstance(raw_roaring[0], np.ndarray)
    
    ## test roaringMatrix getter when roaring is not available
    with pytest.raises(RoaringError, match="Ardal object was instantialised with 'roaring=False'. Cannot retrieve roaring matrix."):
        get_component_no_roaring.roaringMatrix(decode=False)

