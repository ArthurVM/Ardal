import pytest
import numpy as np

from ardal.utils.exceptions import InvalidGUIDQueryError, ParameterError, RoaringError


def test_snp_count(allele_component):
    """
    Tests the snpCount method of the ArdalGet class.
    """
    counts = allele_component.count()
    assert isinstance(counts, dict)
    assert counts["GUID1"] == 6
    assert counts["GUID10"] == 4
    
    
def test_interval_alleles_simple(allele_component_simdata):
    # Suppose interval "chr1:50-150" should match "chr1.100.A.T"
    intervals = [("chr1", 50, 150)]
    alleles = allele_component.interval(intervals=intervals)
    assert "chr1.100.A.T" in alleles
    assert all(isinstance(a, str) for a in alleles)


def test_interval_alleles_multiple(allele_component_simdata):
    # Should match both chr1.100.A.T and chr1.200.G.C
    intervals = [("chr1", 50, 250)]
    alleles = allele_component.interval(intervals=intervals)
    assert set(["chr1.100.A.T", "chr1.200.G.C"]).issubset(set(alleles))


def test_interval_alleles_no_match(allele_component_simdata):
    # No alleles in this interval
    intervals = [("chr3", 1, 1000)]
    alleles = allele_component.interval(intervals=intervals)
    assert alleles == []


def test_interval_alleles_with_format(allele_component_simdata):
    # Use a custom allele_id_format if supported
    intervals = [("chr2", 100, 200)]
    alleles = allele_component.interval(intervals=intervals, allele_id_format="{chr}.{start}.{ref}.{alt}")
    assert "chr2.150.T.G" in alleles


# If you have a BED file, you can add a test like this:
def test_interval_alleles_with_bed(allele_component_simdata, tmp_path):
    bed_content = "chr2\t300\t301\tchr2.300.C.A\n"
    bed_file = tmp_path / "test.bed"
    bed_file.write_text(bed_content)
    intervals = [("chr2", 250, 350)]
    alleles = allele_component.interval(intervals=intervals, intervals_bed=str(bed_file))
    assert "chr2.300.C.A" in alleles
    