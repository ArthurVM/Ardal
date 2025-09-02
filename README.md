# Ardal: A Package for Allele Matrix Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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
    [1, 0, 1, 1, 1],  # GUID1
    [1, 0, 1, 1, 0],  # GUID2
    [1, 1, 1, 1, 0],  # GUID3
    [0, 1, 1, 0, 0],  # GUID4
    [0, 0, 0, 1, 0],  # GUID5
    [1, 0, 0, 1, 0],  # GUID6
], dtype=np.uint8)
headers = {
    "guids" : ["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6"],
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
hamming_matrix = ard_obj.compare.pairwise(metric="hamming")

## calculate a Jaccard distance matrix
jaccard_matrix = ard_obj.compare.pairwise(metric="jaccard")
```
### Compute the SNP Neighbourhood of a Sample
A SNP neighbourhood can be computed by finding all GUIDs which lie within n SNPs of a query GUID.
This can be achieved with or without SIMD, depending on CPU architecture.
```
## find the neighborhood of GUID1 within a Hamming distance of 2 (using SIMD)
neighborhood_simd = ard_obj.compare.neighbourhood("GUID1", n=2, simd=True)

## find the neighborhood of GUID1 within a Hamming distance of 2 (without SIMD)
neighborhood = ard_obj.compare.neighbourhood("GUID1", n=2, simd=False)
```
### Identifying Unique Alleles
You can find the alleles which are unique to a given set of GUIDs within the input data.
```
## find the SNPs found only in both GUID1 and GUID2
unique_snps = ard_obj.allele.unique(["GUID1", "GUID2"])
```
### Stats
Ardal supports some pretty useful and powerful techniques for statistical investigation of the database.
For example, I might want to find a set of alleles which can be used to model a population - perhaps
ones which represent a lineage or strain. This can be achieved using information theoretic functions such 
as `.stats.snpInform()`.
```
## define our set of ingroup guids
ingroup_guids = ["GUID1", "GUID2", "GUID3"]

## calculate the Jensen-Shannon divergence of each SNP
snps = ard.stats.snpInform(ingroup_guids, metric="jensenshannon")

## select the top 1% of SNPs
t = np.percentile(list(snps.values()), 99.999)
top_snps = [snp for snp, score in snps.items() if score >= t]

## compute the allele co-occurrence to ensure co-occurring SNPs dont result in overfitting
cooc_snps = ard.alleleCooc(selected_snps, threshold=0.95, threads=10)
na_snps = [snp for d in list(cooc_snps.values()) for snp in d]
model_snps = [snp for snp in top_snps if snp not in na_snps]
```
### Utlity Functions
There are a number of general utility functions for getting data, input/output, subsetting, and database inspection.
```
## get the allele matrix stats as a JSON
stats_json = ard_obj.get.matrixStats()

## or print as a rich text table
ard_obj.get.matrixStats(print_stats=True)

## get the allele matrix as a NumPy array
bit_matrix = ard_obj.get.bitMatrix()

## get the roaring matrix
roaring_matrix = ard_obj.get.roaringMatrix()

## get the GUID and Allele headers
headers = ard_obj.get.getHeaders()

## get the allele matrix as a Pandas DataFrame
df = ard_obj.io.toDataFrame()

## or as a dictionary (which will be identical to the .get.roaringMatrix() function in practice)
allele_dict = ard_obj.io.toDict()

## write the allele matrix to disk as a npy JSON pair
ard_obj.io.write.snpCount(file_path="./, output_prefix="out_db")

## or compress the npy to npz
ard_obj.io.write.snpCount(file_path="./, output_prefix="out_db", npz=True)
```

# Ardal API Documentation

The Ardal API documentation can be found [here](./docs/api.md)
