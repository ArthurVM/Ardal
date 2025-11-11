import pytest
import numpy as np

from ardal.utils.exceptions import InvalidGUIDQueryError, ParameterError, RoaringError


ALLELE_ID_FORMAT = "{chr}.{start}.{ref}.{alt}"

def test_unique(allele_component_simdata):
    """
    Tests the unique method of the ArdalGet class.
    """
    unique = allele_component_simdata.unique(["GUID100"])
    assert unique == {"GUID100" : {'chr1.100.A.T'}}
    

def test_unique_core(allele_component):
    """
    Tests the unique method of the ArdalGet class.
    """
    unique = allele_component.unique_core(["GUID1", "GUID2"])
    assert unique == ['chr1.1.A.T']
    
    unique = allele_component.unique_core(["GUID3", "GUID4", "GUID5"])
    assert unique == ['chr1.2.A.T']
    
    unique = allele_component.unique_core(["GUID6", "GUID7"])
    assert unique == ['chr1.3.A.T']
    
    unique = allele_component.unique_core(["GUID8", "GUID9", "GUID10"])
    assert unique == ['chr1.4.A.T']
    
    unique = allele_component.unique_core(["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7", "GUID8", "GUID9", "GUID10"])
    assert sorted(unique, key = lambda x: int(x.split(".")[1])) == ['chr1.5.A.T', 'chr1.7.A.T']
    
    unique = allele_component.unique_core(["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7"])
    assert unique == ['chr1.6.A.T']
    
    unique = allele_component.unique_core(["GUID1", "GUID2", "GUID3", "GUID4", "GUID5"])
    assert unique == ['chr1.8.A.T']
    
    unique = allele_component.unique_core(["GUID6", "GUID7", "GUID8", "GUID9", "GUID10"])
    assert unique == ['chr1.9.A.T']
    
    unique = allele_component.unique_core(["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7", "GUID8", "GUID9"])
    assert unique == ['chr1.10.A.T']


def test_core_genome(allele_component):
    """
    Tests the new core genome helper.
    """
    expected = {'chr1.1.A.T', 'chr1.5.A.T', 'chr1.6.A.T', 'chr1.7.A.T', 'chr1.8.A.T', 'chr1.10.A.T'}
    core = allele_component.core(["GUID1", "GUID2"])
    assert core == expected

    relaxed_core = allele_component.core(["GUID1", "GUID2", "GUID3"], missingness=0.34)
    assert 'chr1.1.A.T' in relaxed_core
    assert 'chr1.2.A.T' not in relaxed_core

    core_counts = allele_component.core(["GUID1", "GUID2"], return_counts=True)
    assert core_counts['chr1.1.A.T'] == 2
    assert core_counts['chr1.7.A.T'] == 2


def test_pan_genome(allele_component):
    """
    Tests the pan genome helper.
    """
    pan = allele_component.pan(["GUID1", "GUID3"])
    expected = {'chr1.1.A.T', 'chr1.2.A.T', 'chr1.5.A.T', 'chr1.6.A.T', 'chr1.7.A.T', 'chr1.8.A.T', 'chr1.10.A.T'}
    assert pan == expected

    pan_counts = allele_component.pan(["GUID6", "GUID7"], return_counts=True)
    assert pan_counts['chr1.3.A.T'] == 2
    assert pan_counts['chr1.9.A.T'] == 2


def test_genome_summary(allele_component):
    """
    Tests the combined genome summary helper.
    """
    summary = allele_component.genome_summary(["GUID1", "GUID2"])
    assert summary["core"] == {'chr1.1.A.T', 'chr1.5.A.T', 'chr1.6.A.T', 'chr1.7.A.T', 'chr1.8.A.T', 'chr1.10.A.T'}
    assert summary["pan"] == summary["core"]
    assert summary["unique"] == set(allele_component.unique_core(["GUID1", "GUID2"]))


def test_guids_with_alleles(allele_component):
    """
    Tests the guids_with_alleles method of the ArdalGet class.
    """
    guids = allele_component.guids_with_alleles(['chr1.1.A.T'])
    assert sorted(guids, key=lambda x: int(x.split("GUID")[1])) == ["GUID1", "GUID2"]
    
    guids = allele_component.guids_with_alleles(['chr1.2.A.T'])
    assert sorted(guids, key=lambda x: int(x.split("GUID")[1])) == ["GUID3", "GUID4", "GUID5"]
    
    guids = allele_component.guids_with_alleles(['chr1.3.A.T'])
    assert sorted(guids, key=lambda x: int(x.split("GUID")[1])) == ["GUID6", "GUID7"]
    
    guids = allele_component.guids_with_alleles(['chr1.4.A.T'])
    assert sorted(guids, key=lambda x: int(x.split("GUID")[1])) == ["GUID8", "GUID9", "GUID10"]
    
    guids = allele_component.guids_with_alleles(['chr1.5.A.T', 'chr1.7.A.T'])
    assert sorted(guids, key=lambda x: int(x.split("GUID")[1])) == ["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7", "GUID8", "GUID9", "GUID10"]
    
    guids = allele_component.guids_with_alleles(['chr1.6.A.T'])
    assert sorted(guids, key=lambda x: int(x.split("GUID")[1])) == ["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7"]
    
    guids = allele_component.guids_with_alleles(['chr1.8.A.T'])
    assert sorted(guids, key=lambda x: int(x.split("GUID")[1])) == ["GUID1", "GUID2", "GUID3", "GUID4", "GUID5"]
    
    guids = allele_component.guids_with_alleles(['chr1.9.A.T'])
    assert sorted(guids, key=lambda x: int(x.split("GUID")[1])) == ["GUID6", "GUID7", "GUID8", "GUID9", "GUID10"]
    
    guids = allele_component.guids_with_alleles(['chr1.10.A.T'])
    assert sorted(guids, key=lambda x: int(x.split("GUID")[1])) == ["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6", "GUID7", "GUID8", "GUID9"]


def test_match_names(allele_component):
    """
    Tests the match_names method of the ArdalGet class.
    """
    alleles = allele_component.match_names("chr1*")
    assert sorted(alleles, key = lambda x: int(x.split(".")[1])) == [f"chr1.{i}.A.T" for i in range(1, 11)]
    
    alleles = allele_component.match_names("*.A.*")
    assert sorted(alleles, key = lambda x: int(x.split(".")[1])) == [f"chr1.{i}.A.T" for i in range(1, 11)]
    
    alleles = allele_component.match_names("*5.*.T")
    assert sorted(alleles, key = lambda x: int(x.split(".")[1])) == ["chr1.5.A.T"]
    
    alleles = allele_component.match_names("*.1*.T")
    assert sorted(alleles, key = lambda x: int(x.split(".")[1])) == ["chr1.1.A.T", "chr1.10.A.T"]


def test_allele_count(allele_component):
    """
    Tests the count method of the ArdalGet class.
    """
    counts = allele_component.count()
    assert isinstance(counts, dict)
    assert counts["GUID1"] == 6
    assert counts["GUID10"] == 4
    
    
def test_interval_alleles_simple(allele_component_simdata):
    # Suppose interval "chr1:50-150" should match "chr1.100.A.T"
    intervals = [("chr1", 100, 101)]
    alleles = allele_component_simdata.interval(intervals=intervals, allele_id_format=ALLELE_ID_FORMAT)
    assert alleles == ["chr1.100.A.T"]
    assert all(isinstance(a, str) for a in alleles)


def test_interval_alleles_multiple(allele_component_simdata):
    # Should match both chr1.100.A.T and chr1.200.G.C
    intervals = [("chr1", 10, 20)]
    alleles = allele_component_simdata.interval(intervals=intervals, allele_id_format=ALLELE_ID_FORMAT)
    assert sorted(alleles, key = lambda x: int(x.split(".")[1])) == [f"chr1.{i}.A.T" for i in range(10,21)]
    

def test_interval_alleles_polyinterval(allele_component_simdata):
    # Should match both chr1.100.A.T and chr1.200.G.C
    intervals = [("chr1", 10, 20), ["chr1", 90, 100]]
    alleles = allele_component_simdata.interval(intervals=intervals, allele_id_format=ALLELE_ID_FORMAT)
    assert sorted(alleles, key = lambda x: int(x.split(".")[1])) == [f"chr1.{i}.A.T" for i in range(10,21)] + [f"chr1.{i}.A.T" for i in range(90,101)]
        

def test_interval_alleles_no_match(allele_component_simdata):
    # No alleles in this interval
    intervals = [("chr3", 1, 1000)]
    alleles = allele_component_simdata.interval(intervals=intervals, allele_id_format=ALLELE_ID_FORMAT)
    assert alleles == []


def test_interval_alleles_chr(allele_component_simdata):
    # Use a custom allele_id_format if supported
    intervals = [["chr1"]]
    alleles = allele_component_simdata.interval(intervals=intervals, allele_id_format=ALLELE_ID_FORMAT)
    assert sorted(alleles, key = lambda x: int(x.split(".")[1])) == allele_component_simdata._headerUtils.headers["alleles"]


# If you have a BED file, you can add a test like this:
def test_interval_alleles_with_coords_bed(allele_component_simdata, tmp_path):
    bed_content = "chr2\t300\t301\tchr1.100.A.T\n"
    bed_file = tmp_path / "test_allele_coords.bed"
    bed_file.write_text(bed_content)
    intervals = [("chr2", 250, 350)]
    alleles = allele_component_simdata.interval(intervals=intervals, allele_coords_bed=str(bed_file))
    assert alleles == ["chr1.100.A.T"]
    
    
# If you have a BED file, you can add a test like this:
def test_interval_alleles_with_intervals_bed(allele_component_simdata, tmp_path):
    bed_content = "chr1\t10\t20\nchr1\t50\t60\nchr1\t80\t90\n"
    bed_file = tmp_path / "test_intervals.bed"
    bed_file.write_text(bed_content)
    alleles = allele_component_simdata.interval(intervals=None, intervals_bed=str(bed_file))
    assert sorted(alleles, key = lambda x: int(x.split(".")[1])) == [f"chr1.{i}.A.T" for i in range(10,21)] + [f"chr1.{i}.A.T" for i in range(50,61)] + [f"chr1.{i}.A.T" for i in range(80,91)]


def test_interval_alleles_with_dual_interval_input(allele_component_simdata, tmp_path):
    bed_content = "chr1\t10\t20\nchr1\t80\t90\n"
    bed_file = tmp_path / "test_intervals.bed"
    bed_file.write_text(bed_content)
    intervals = [("chr1", 50, 60)]
    alleles = allele_component_simdata.interval(intervals=intervals, intervals_bed=str(bed_file))
    assert sorted(alleles, key = lambda x: int(x.split(".")[1])) == [f"chr1.{i}.A.T" for i in range(10,21)] + [f"chr1.{i}.A.T" for i in range(50,61)] + [f"chr1.{i}.A.T" for i in range(80,91)]


def test_get_positions(allele_component, tmp_path):
    """
    Tests the count method of the ArdalGet class.
    """
    ## ensure alleles have been positionally decoded
    intervals = [("chr1", 1, 10)]
    alleles = allele_component.interval(intervals=intervals, allele_id_format=ALLELE_ID_FORMAT)
    
    ## general test
    positions = allele_component.get_positions()
    assert positions == {'chr1.1.A.T': ('chr1', 1),
                         'chr1.2.A.T': ('chr1', 2),
                         'chr1.3.A.T': ('chr1', 3),
                         'chr1.4.A.T': ('chr1', 4),
                         'chr1.5.A.T': ('chr1', 5),
                         'chr1.6.A.T': ('chr1', 6),
                         'chr1.7.A.T': ('chr1', 7),
                         'chr1.8.A.T': ('chr1', 8),
                         'chr1.9.A.T': ('chr1', 9),
                         'chr1.10.A.T': ('chr1', 10)}
    
    ## adjust position of allele chr1.10.A.T to ensure overwrite happens as expected
    bed_content = "chr2\t300\t301\tchr1.10.A.T\n"
    bed_file = tmp_path / "test_allele_coords.bed"
    bed_file.write_text(bed_content)
    alleles = allele_component.interval(intervals=intervals, allele_coords_bed=str(bed_file))
    
    ## test overwrite
    positions = allele_component.get_positions()
    assert positions["chr1.10.A.T"] == ("chr2", 300)
