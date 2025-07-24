# Ardal: A Package for Allele Matrix Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/ArthurVM/Ardal/actions/workflows/conda.yml/badge.svg)](https://github.com/ArthurVM/Ardal/actions/workflows/conda.yml/)

**Ardal** is a Python package designed for efficient analysis of allele matrices, particularly in the context of genomics and population genetics. It provides tools for calculating distances between samples, identifying core and accessory alleles, and performing other common operations on allele data. Ardal leverages a C++ backend for performance-critical operations, using `pybind11` to create seamless Python bindings.

## Features

*   **Allele Matrix Representation:** Efficiently stores and manipulates binary allele matrices.
*   **Distance Calculations:**
    *   Hamming distance: Calculates the number of differing alleles between samples.
    *   Jaccard distance: Measures the dissimilarity between samples based on shared alleles.
*   **Neighborhood Analysis:** Identifies samples within a specified Hamming distance of a target sample, with or without SIMD acceleration.
*   **Core and Accessory Allele Identification:** Determines the core (present in all samples) and accessory (present in some samples) alleles within a group of samples.
*   **Unique Allele Identification:** Identifies alleles that are unique to a specific set of samples.
*   **Data Input:** Supports various data formats, including:
    *   CSV
    *   Parquet
    *   NPY/JSON pairs
    *   In-memory data structures (NumPy arrays and dictionaries)
*   **Data Output:**
    *   NumPy arrays
    *   Pandas DataFrames
    *   Dictionaries
*   **Caching:** Implements a caching mechanism to speed up repeated Hamming distance calculations.
* **SIMD**: Uses SIMD instructions to accelerate neighbourhood calculations.

## Installation

Ardal requires a C++ compiler to build the C++ extension.

1.  **Using `pip`:**

    ```bash
    pip install ardal
    ```

    This will install the package and its dependencies.

2.  **From Source:**

    ```bash
    git clone https://github.com/ArthurVM/ardal.git
    cd ardal
    pip install .
    ```

    This will clone the repository, navigate to the project directory, and install the package.

3. **Using Conda**
    ```bash
    conda build .
    conda install /path/to/your/conda-bld/noarch/ardal-0.1.0-py_0.tar.bz2
    ```
    This will build and install the package using conda.

## Usage

Some example usage for `ardal`:

### Creating an `Ardal` Object
From data in memory
```
from ardal import Ardal

## create some dummy data
matrix_data = np.array([
    [1, 0, 1, 0, 1],  # GUID1
    [1, 0, 1, 0, 0],  # GUID2
    [1, 0, 1, 1, 0],  # GUID3
    [0, 1, 1, 0, 0],  # GUID4
    [0, 0, 0, 1, 0],  # GUID5
], dtype=np.uint8)
headers = {
    "guids" : ["GUID1", "GUID2", "GUID3", "GUID4", "GUID5"],
    "alleles" : ["SNP1", "SNP2", "SNP3", "SNP4", "SNP5"]
}

## create an Ardal object from dummy data
ard_obj = Ardal(data_source=[headers, matrix_data])
```
From data on disk
```
headers_path = "/path/to/headers.json"
matrix_path = "/path/to/matrix.npy"
ard_obj = Ardal(data_source=[headers_path, matrix_path])
```
### Compute a Distance Matrix
You can compute a distance matrix easily using either Hamming or Jaccard distance.
```
## calculate a Hamming distance matrix
hamming_matrix = obj.pairwise(metric="hamming")

## calculate a Jaccard distance matrix
jaccard_matrix = obj.pairwise(metric="jaccard")
```
### Compute the SNP Neighbourhood of a Sample
A SNP neighbourhood can be computed by finding all GUIDs which lie within n SNPs of a query GUID.
This can be achieved with or without SIMD, depending on CPU architecture.
```
## find the neighborhood of GUID1 within a Hamming distance of 2 (using SIMD)
neighborhood_simd = obj.neighbourhood("GUID1", n=2, simd=True)

## find the neighborhood of GUID1 within a Hamming distance of 2 (without SIMD)
neighborhood = obj.neighbourhood("GUID1", n=2, simd=False)
```
### Identify Core and Accessory Alleles
You can calculate the core and accessory alleles for a given set of GUIDs.
```
## find the core alleles for GUID1 and GUID2
core_alleles = obj.core(["GUID1", "GUID2"])

## find the accessory alleles for GUID1 and GUID2
accessory_alleles = obj.accessory(["GUID1", "GUID2"])
```
### Identifying Unique Alleles
You can find the alleles which are unique to a given set of GUIDs within the input data.
```
## find the SNPs found only in both GUID1 and GUID2
unique_snps = obj.unique(["GUID1", "GUID2"])
```
### Other Functions
```
## get the allele matrix as a NumPy array
matrix = obj.getMatrix()

## get the GUID and Allele headers
headers = obj.getHeaders()

## get the allele matrix as a Pandas DataFrame
df = obj.toDataFrame()

## get the SNP count for each GUID
snp_counts = obj.snpCount()

## get stats about the Ardal object
stats = obj.stats()

## flush the cache
obj.flushCache()
```

# Ardal API Documentation

The Ardal API documentation can be found [here](./docs/api.md)