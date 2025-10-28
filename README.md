# Ardal: A Package for Allele Matrix Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/ArthurVM/Ardal/actions/workflows/conda.yml/badge.svg)](https://github.com/ArthurVM/Ardal/actions/workflows/conda.yml/)

**Ardal** is a Python package designed for efficient analysis of allele matrices, particularly in the context of genomics and population genetics. It provides tools for calculating distances between samples, identifying core and accessory alleles, and performing other common operations on allele data. Ardal leverages a C++ backend for performance-critical operations, using `pybind11` to create seamless Python bindings.

## Features

*   **Allele Matrix Representation:** Efficiently stores and manipulates binary allele matrices.
*   **Distance Calculations:**
    *   Hamming distance: Counts the number of positions at which two allele vectors differ. It is best for revealing overall genetic divergence or mutation load between genomes.
    *   Jaccard index: Measures similarity as the ratio of shared alleles to the union of alleles across two vectors. It is best for revealing shared allele patterns while down-weighting invariant or absent sites.
    *   Inner Product: Quantifies the raw co-occurrence of alleles between two binary vectors. It is best for highlighting densely co-occurring alleles or lineage-defining variant constellations.
    *   Cosine distance: Measures angular dissimilarity by normalizing the inner product by the vector magnitudes. It is best for revealing directional similarity in variant composition independent of absolute allele counts.
*   **Neighborhood Analysis:** Identifies samples within a specified Hamming distance of a target sample, with or without SIMD acceleration.
*   **Core and Accessory Allele Identification:** Determines the core (present in all samples) and accessory (present in some samples) alleles within a group of samples.
*   **Unique Allele Identification:** Identifies alleles that are unique to a specific set of samples.
*   **Data Input:** Supports various data formats, including:
    *   Bitpacked binary matrices
    *   CSV
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
```python
import numpy as np
from ardal import Ardal

# create some dummy data
matrix_data = np.array([
    [1, 0, 1, 1, 1],  # GUID1
    [1, 0, 1, 1, 0],  # GUID2
    [1, 1, 1, 1, 0],  # GUID3
    [0, 1, 1, 0, 0],  # GUID4
    [0, 0, 0, 1, 0],  # GUID5
    [1, 0, 0, 1, 0],  # GUID6
], dtype=np.uint8)
headers = {
    "guids": ["GUID1", "GUID2", "GUID3", "GUID4", "GUID5", "GUID6"],
    "alleles": ["SNP1", "SNP2", "SNP3", "SNP4", "SNP5"],
}

# create an Ardal object from the in-memory data
ard_obj = Ardal(
    data_source=[matrix_data, headers],
    roaring="auto",            # automatically provision the roaring backend if sparse enough
    density_threshold=0.05,    # override the default heuristic if required
    quiet_init=True,           # suppress the initial stats printout
)
```

From data on disk
```python
from ardal import Ardal

headers_path = "/path/to/headers.json"
matrix_path = "/path/to/matrix.npy"

ard_obj = Ardal(data_source=[matrix_path, headers_path])
```
### Compute a Distance Matrix
You can compute a distance matrix easily using either Hamming or Jaccard distance.
```python
# calculate a condensed Hamming distance vector
hamming = ard_obj.distance.pairwise(metric="hamming")

# return a square pandas DataFrame using the roaring backend where available
jaccard_df = ard_obj.distance.pairwise(metric="jaccard", backend="roaring", return_square=True, as_dataframe=True)
```
### Compute the SNP Neighbourhood of a Sample
A SNP neighbourhood can be computed by finding all GUIDs which lie within n SNPs of a query GUID.
This can be achieved with or without SIMD, depending on CPU architecture.
```python
# find the neighbourhood of GUID1 within a Hamming distance of 2 (using SIMD)
neighbourhood_simd = ard_obj.distance.neighbourhood("GUID1", n=2, use_simd=True)

# find the neighbourhood of GUID1 within a Hamming distance of 2 (without SIMD)
neighbourhood_scalar = ard_obj.distance.neighbourhood("GUID1", n=2, use_simd=False)
```
### Identifying Unique Alleles
You can find the alleles which are unique to a given set of GUIDs within the input data.
```python
# find the SNPs found only in GUID1 and GUID2
unique_snps = ard_obj.allele.unique(["GUID1", "GUID2"])

# alleles present in all of GUID1/GUID2 yet absent elsewhere
unique_core = ard_obj.allele.unique_core(["GUID1", "GUID2"])
```
### Stats
Ardal supports some pretty useful and powerful techniques for statistical investigation of the database.
For example, I might want to find a set of alleles which can be used to model a population - perhaps
ones which represent a lineage or strain. This can be achieved using information theoretic functions such 
as `.stats.snpInform()`.
```python
import numpy as np

# define our set of ingroup guids
ingroup_guids = ["GUID1", "GUID2", "GUID3"]

# calculate the Jensen–Shannon divergence of each SNP
scores = ard_obj.stats.allele_inform(ingroup_guids, method="jensenshannon")

# select the top 1% of SNPs
threshold = np.percentile(list(scores.values()), 99)
top_snps = [snp for snp, score in scores.items() if score >= threshold]

# compute the allele co-occurrence to ensure co-occurring SNPs don't result in overfitting
cooc = ard_obj.stats.allele_cooc(top_snps, threshold=0.95, threads=4)
na_snps = {snp for neighbours in cooc.values() for snp in neighbours}
model_snps = [snp for snp in top_snps if snp not in na_snps]
```
### Utlity Functions
There are a number of general utility functions for getting data, input/output, subsetting, and database inspection.
```python
# get the allele matrix stats as a dictionary
stats = ard_obj.get.matrix_stats()

# or print as a rich text table
ard_obj.get.matrix_stats(print_table=True)

# get the allele matrix as a NumPy array
bit_matrix = ard_obj.get.bit_matrix()

# get the roaring matrix (decoded to allele IDs)
roaring_matrix = ard_obj.get.roaring_matrix()

# get the GUID and allele headers
headers = ard_obj.get.headers()

# get the allele matrix as a Pandas DataFrame
df = ard_obj.io.to_dataframe()

# or as a dictionary (optionally choose the backend)
allele_dict = ard_obj.io.to_dict(backend="auto")

# write the allele matrix to disk
ard_obj.io.write(output_prefix="out_db", out_directory="./exports", format="npz")
```

# Ardal API Documentation

The Ardal API documentation can be found [here](./docs/api.md)
