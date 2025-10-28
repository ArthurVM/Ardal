import os
import pytest

from ardal.utils.exceptions import InvalidAlleleQueryError, InvalidGUIDQueryError


def test_decode_headers(headerUtils_component):
    decoded = headerUtils_component._decoded_headers
    assert decoded["guids"]["GUID1"] == 0
    assert decoded["alleles"]["chr1.1.A.T"] == 0


def test_encode_decode_guid(headerUtils_component):
    idx = headerUtils_component.encode_guid("GUID2")
    assert idx == 1
    guid = headerUtils_component.decode_guid(2)
    assert guid == "GUID3"


def test_encode_decode_allele(headerUtils_component):
    idx = headerUtils_component.encode_allele("chr1.3.A.T")
    assert idx == 2
    allele = headerUtils_component.decode_allele(0)
    assert allele == "chr1.1.A.T"


def test_check_guids_valid(headerUtils_component):
    headerUtils_component.check_guids(["GUID1", "GUID2"])  # Should not raise


def test_check_guids_invalid(headerUtils_component):
    with pytest.raises(InvalidGUIDQueryError):
        headerUtils_component.check_guids(["GUID1", "GUIDX"])


def test_check_alleles_valid(headerUtils_component):
    headerUtils_component.check_alleles(["chr1.1.A.T", "chr1.2.A.T"])  # Should not raise


def test_check_alleles_invalid(headerUtils_component):
    with pytest.raises(InvalidAlleleQueryError):
        headerUtils_component.check_alleles(["chr1.1.A.T", "chrX.999.A.C"])


def test_decode_allele_position(headerUtils_component):
    allele_id_format = "{chr}.{start}.{ref}.{alt}"
    pattern = headerUtils_component._check_allele_format_grammar(allele_id_format=allele_id_format)
    chr_, start, end, ref, alt = headerUtils_component._decode_allele_position(allele_id="chr1.100.A.T", allele_id_format=allele_id_format, pattern=pattern)
    assert chr_ == "chr1"
    assert start == 100
    assert ref == "A"
    assert alt == "T"


def test_get_allele_positions(headerUtils_component):
    headerUtils_component.compute_allele_positions(allele_id_format="{chr}.{start}.{ref}.{alt}")
    pos = headerUtils_component.get_allele_positions()
    print(pos)
    assert pos["chr1.1.A.T"] == ("chr1", 1)
    assert pos["chr1.2.A.T"] == ("chr1", 2)
    assert pos["chr1.3.A.T"] == ("chr1", 3)
    assert pos["chr1.4.A.T"] == ("chr1", 4)
    assert pos["chr1.5.A.T"] == ("chr1", 5)
    assert pos["chr1.6.A.T"] == ("chr1", 6)
    assert pos["chr1.7.A.T"] == ("chr1", 7)
    assert pos["chr1.8.A.T"] == ("chr1", 8)
    assert pos["chr1.9.A.T"] == ("chr1", 9)
    
    ## we changed this earlier in test_ArdalAllele:test_get_positions
    ## can be commented out if it causes issues due to run order
    assert pos["chr1.10.A.T"] == ("chr2", 300)


def test_read_coords_bed(tmp_path, headerUtils_component):
    bed_content = "chr1 100 101 chr1.100.A.T\nchr2 150 151 chr2.150.T.G\n"
    bed_file = tmp_path / "coords.bed"
    bed_file.write_text(bed_content)
    coords = headerUtils_component._read_allele_coords_bed(str(bed_file))
    assert coords["chr1.100.A.T"] == ["chr1", 100, 101]
    assert coords["chr2.150.T.G"] == ["chr2", 150, 151]