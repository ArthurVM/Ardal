import os
import sys
import pytest
import numpy as np
import json
from ardal import Ardal
from ardal.core.ArdalHeaderUtils import ArdalHeaderUtils
from ardal.core.ArdalIO import ArdalIO
from ardal.core.ArdalGet import ArdalGet
from ardal.core.ArdalStats import ArdalStats
from ardal.core.ArdalAllele import ArdalAllele
from ardal.core.ArdalCompare import ArdalCompare


@pytest.fixture(scope="session")
def test_data_mem():
    """
    Creates a simulated test dataset for testing the Ardal functionality.

    OUTPUT:
        tuple: A tuple containing the matrix (NumPy array) and headers (dictionary).
    """

    ## define GUIDs and SNPs
    guids = [f"GUID{i}" for i in range(1, 11)]  # GUID1 to GUID10
    alleles = [f"SNP{i}" for i in range(1, 11)]  # SNP1 to SNP10

    matrix = np.zeros((len(guids), len(alleles)), dtype=np.uint8)

    ## set SNPs
    matrix[0, 0] = 1  ## SNP1 in GUID1
    matrix[1, 0] = 1  ## SNP1 in GUID2
    matrix[2, 1] = 1  ## SNP2 in GUID3
    matrix[3, 1] = 1  ## SNP2 in GUID4
    matrix[4, 1] = 1  ## SNP2 in GUID5
    matrix[5, 2] = 1  ## SNP3 in GUID6
    matrix[6, 2] = 1  ## SNP3 in GUID7
    matrix[7, 3] = 1  ## SNP4 in GUID8
    matrix[8, 3] = 1  ## SNP4 in GUID9
    matrix[9, 3] = 1  ## SNP4 in GUID10

    ## SNP5 in all GUIDs
    matrix[:, 4] = 1

    ## SNP6 in GUID1-7
    matrix[0:7, 5] = 1

    ## SNP7 in all GUIDs
    matrix[:, 6] = 1

    ## SNP8 in GUID1-5
    matrix[0:5, 7] = 1

    ## SNP9 in GUID6-10
    matrix[5:10, 8] = 1

    ## SNP10 present in GUID1-9
    matrix[0:9, 9] = 1

    ## header json
    headers = {
        "guids": guids,
        "alleles": alleles
    }

    return [matrix, headers]


@pytest.fixture(scope="session")
def test_data_matrix_npy():
    return ["./data/sim_matrix.npy", "./data/sim_headers.json"]


@pytest.fixture(scope="session")
def test_data_matrix_npz():
    return ["./data/sim_matrix.npz", "./data/sim_headers.json"]


@pytest.fixture(scope="session")
def test_data_matrix_csv():
    return "./data/sim_matrix.csv"


@pytest.fixture(scope="session")
def ardal_object(test_data_mem):
    """
    Creates an Ardal object for testing.
    """
    ardal_object = Ardal(data_source=test_data_mem, roaring=True, quiet=True)
    return ardal_object


@pytest.fixture(scope="session")
def hybrid_matrix(test_data_mem, ardal_object):
    """Creates a HybridMatrix object for testing."""
    # Assuming default parameters from the Ardal constructor
    return ardal_object.get.hybridMatrix()


@pytest.fixture(scope="session")
def headerUtils_component(test_data_mem):
    """Creates a headerUtils object for testing."""
    matrix, headers = test_data_mem
    return ArdalHeaderUtils(headers)


@pytest.fixture(scope="session")
def io_component(headerUtils_component, hybrid_matrix):
    """Creates an ArdalIO component for unit testing."""
    return ArdalIO(headerUtils_component, hybrid_matrix, hybrid_matrix.roaringEnabled())


@pytest.fixture(scope="session")
def get_component(headerUtils_component, hybrid_matrix):
    """Creates an ArdalGet component for unit testing."""
    return ArdalGet(headerUtils_component, hybrid_matrix, hybrid_matrix.roaringEnabled())


@pytest.fixture(scope="session")
def stats_component(headerUtils_component, hybrid_matrix):
    """Creates an ArdalStats component for unit testing."""
    return ArdalStats(headerUtils_component, hybrid_matrix, hybrid_matrix.roaringEnabled())


@pytest.fixture(scope="session")
def allele_component(headerUtils_component, hybrid_matrix):
    """Creates an ArdalAllele component for unit testing."""
    return ArdalAllele(headerUtils_component, hybrid_matrix, hybrid_matrix.roaringEnabled())


@pytest.fixture(scope="session")
def compare_component(headerUtils_component, hybrid_matrix):
    """Creates an ArdalCompare component for unit testing."""
    return ArdalCompare(headerUtils_component, hybrid_matrix, hybrid_matrix.roaringEnabled())
