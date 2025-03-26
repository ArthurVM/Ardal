import sys
print(sys.path)
import pytest
import numpy as np
import json
from ardal import Ardal, ArdalParser


@pytest.fixture(scope="session")
def test_data():
    """
    Creates a simulated test dataset for testing the Ardal.unique functionality.

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
def ardal_object(test_data):
    """
    Creates an Ardal object for testing.
    """
    ardal_object = Ardal(data_source=test_data)
    return ardal_object
