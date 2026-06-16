import pytest
import numpy as np

from ardal import Ardal
from ardal.core.ArdalHeaderUtils import ArdalHeaderUtils
from ardal.core.ArdalGet import ArdalGet
from ardal.utils.make_meta import make_meta
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
    allele_subset = ["chr1.1.A.T", "chr1.10.A.T"]
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
    Tests the matrix_stats method of the ArdalGet class.
    """
    stats = get_component.matrix_stats(print_table=False)
    assert isinstance(stats, dict)
    assert "Number of GUIDs" in stats
    assert "Matrix Density" in stats
    assert stats["Number of GUIDs"] == 10
    assert stats["Number of Alleles"] == 10

    # Test printing
    get_component.matrix_stats(print_table=True)
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
    np.testing.assert_array_equal(get_component.bit_matrix(), matrix)


def test_roaring_matrix(get_component, get_component_no_roaring):
    """
    Tests the roaringMatrix method of the ArdalGet class.
    """
    ## test decoded output
    roaring_dict = get_component.roaring_matrix(decode=True)
    assert isinstance(roaring_dict, dict)
    assert "GUID1" in roaring_dict
    assert set(roaring_dict["GUID1"]) == {"chr1.1.A.T", "chr1.5.A.T", "chr1.6.A.T", "chr1.7.A.T", "chr1.8.A.T", "chr1.10.A.T"}

    ## test raw output
    raw_roaring = get_component.roaring_matrix(decode=False)
    assert isinstance(raw_roaring, list)
    assert len(raw_roaring) == 10
    assert isinstance(raw_roaring[0], np.ndarray)
    
    ## test roaringMatrix getter when roaring is not available
    with pytest.raises(RoaringError, match="Ardal object was instantialised with 'roaring=False'. Cannot retrieve roaring matrix."):
        get_component_no_roaring.roaring_matrix(decode=False)


def test_missing_sites_with_positions(hybrid_matrix, headers, meta):
    missing_masks = {"column_masks": {"GUID1": [0, 7], "GUID2": [0], "GUID10": []}}
    header_utils = ArdalHeaderUtils(
        headers=headers,
        meta=meta,
        allele_id_format="{chr}.{start}.{ref}.{alt}",
        missing_masks=missing_masks,
    )
    get_component = ArdalGet(header_utils, hybrid_matrix, hybrid_matrix.roaringEnabled())

    out = get_component.missing_sites()
    assert out["GUID1"] == ["chr1:1", "chr1:8"]
    assert out["GUID2"] == ["chr1:1"]
    assert out["GUID10"] == []


def test_missing_sites_without_positions(hybrid_matrix, headers, meta):
    missing_masks = {"column_masks": {"GUID1": [0, 7], "GUID2": [0], "GUID10": []}}
    header_utils = ArdalHeaderUtils(
        headers=headers,
        meta=meta,
        allele_id_format=None,
        missing_masks=missing_masks,
    )
    get_component = ArdalGet(header_utils, hybrid_matrix, hybrid_matrix.roaringEnabled())

    out = get_component.missing_sites()
    assert out["GUID1"] == ["chr1.1.A.T", "chr1.8.A.T"]
    assert out["GUID2"] == ["chr1.1.A.T"]
    assert out["GUID10"] == []


def test_subset_carries_missing_masks(hybrid_matrix, headers, meta):
    missing_masks = {"column_masks": {"GUID1": [0, 7], "GUID2": [0], "GUID10": [7]}}
    header_utils = ArdalHeaderUtils(
        headers=headers,
        meta=meta,
        allele_id_format="{chr}.{start}.{ref}.{alt}",
        missing_masks=missing_masks,
    )
    get_component = ArdalGet(header_utils, hybrid_matrix, hybrid_matrix.roaringEnabled())

    _, headers_meta = get_component.subset(
        guids=["GUID1", "GUID2", "GUID10"],
        alleles=["chr1.1.A.T", "chr1.10.A.T"],
        drop_zero_cols=False,
        data_only=True,
    )

    assert "column_masks" in headers_meta
    assert headers_meta["column_masks"]["GUID1"] == [0]
    assert headers_meta["column_masks"]["GUID2"] == [0]
    assert headers_meta["column_masks"]["GUID10"] == []


def test_subset_child_missing_sites_with_position_cache(test_data_mem):
    matrix, headers = test_data_mem
    missing = {"GUID1": [0], "GUID2": [1]}
    headers_meta = {
        "meta": make_meta(
            matrix,
            headers,
            generated_by="test",
            format_name="ardal.dense.v1",
            matrix_file=None,
        ),
        "headers": headers,
        "column_masks": missing,
    }

    ard = Ardal(
        data_source=[matrix, headers_meta],
        roaring=True,
        quiet_init=True,
        allele_id_format="{chr}.{start}.{ref}.{alt}",
    )
    child = ard.get.subset(
        guids=["GUID1", "GUID2"],
        alleles=["chr1.1.A.T", "chr1.2.A.T"],
        drop_zero_cols=False,
    )
    out = child.get.missing_sites()
    assert out["GUID1"] == ["chr1:1"]
    assert out["GUID2"] == ["chr1:2"]


def test_subset_passes_child_ardal_kwargs(test_data_mem):
    ard = Ardal(
        data_source=test_data_mem,
        roaring=True,
        quiet_init=True,
        allele_id_format="{chr}.{start}.{ref}.{alt}",
    )

    child = ard.get.subset(
        guids=["GUID1", "GUID2"],
        alleles=["chr1.1.A.T", "chr1.2.A.T"],
        drop_zero_cols=False,
        child_ardal_kwargs={
            "roaring": False,
            "density_threshold": 0.75,
            "allele_positions": {"custom": {123: ["chr1.1.A.T"]}},
        },
    )

    assert child.build_roaring is False
    assert child.density_threshold == 0.75
    assert child._headerUtils.allele_positions == {"custom": {123: ["chr1.1.A.T"]}}
