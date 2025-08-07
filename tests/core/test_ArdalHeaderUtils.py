import os
import pytest

from ardal.utils.exceptions import InvalidAlleleQueryError, InvalidGUIDQueryError


def test_decode_headers(headerUtils_component):
    decoded = headerUtils_component._decoded_headers
    assert decoded["guids"]["g1"] == 0
    assert decoded["alleles"]["chr1.200.G.C"] == 1


def test_encode_decode_guid(headerUtils_component):
    idx = headerUtils_component.encodeGuid("g2")
    assert idx == 1
    guid = headerUtils_component.decodeGuid(2)
    assert guid == "g3"

def test_encode_decode_allele(headerUtils_component):
    idx = headerUtils_component.encodeAllele("chr2.150.T.G")
    assert idx == 2
    allele = headerUtils_component.decodeAllele(0)
    assert allele == "chr1.100.A.T"

def test_check_guids_valid(headerUtils_component):
    headerUtils_component.checkGUIDs(["g1", "g3"])  # Should not raise

def test_check_guids_invalid(headerUtils_component):
    with pytest.raises(InvalidGUIDQueryError):
        headerUtils_component.checkGUIDs(["g1", "gX"])

def test_check_alleles_valid(headerUtils_component):
    headerUtils_component.checkAlleles(["chr1.100.A.T", "chr2.150.T.G"])  # Should not raise

def test_check_alleles_invalid(headerUtils_component):
    with pytest.raises(InvalidAlleleQueryError):
        headerUtils_component.checkAlleles(["chr1.100.A.T", "chrX.999.A.C"])

def test_decode_allele_position(headerUtils_component):
    chr_, start, end, ref, alt = headerUtils_component._decodeAllelePosition("chr1.100.A.T")
    assert chr_ == "chr1"
    assert start == 100
    assert ref == "A"
    assert alt == "T"

def test_get_allele_positions(headerUtils_component):
    pos = headerUtils_component.getAllelePositions(allele_id_format="{chr}.{start}.{ref}.{alt}")
    assert "chr1" in pos
    assert 100 in pos["chr1"]
    assert "chr1.100.A.T" in pos["chr1"][100]

# def test_read_coords_bed(tmp_path, headerUtils_component):
#     bed_content = "chr1 100 101 chr1.100.A.T\nchr2 150 151 chr2.150.T.G\n"
#     bed_file = tmp_path / "coords.bed"
#     bed_file.write_text(bed_content)
#     coords = headerUtils_component._readCoordsBED(str(bed_file))
#     assert coords["chr1.100.A.T"] == ["chr1", "100", "101"]
#     assert coords["chr2.150.T.G"] == ["chr2", "150", "151"]