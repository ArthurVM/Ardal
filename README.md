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

This document provides detailed API documentation for the `ardal` package, including the `Ardal`, `ArdalParser`, and `AlleleMatrix` classes.

## Ardal Class

The `Ardal` class is the main interface for interacting with allele matrix data. It provides methods for calculating distances, identifying core and accessory alleles, and performing other common operations.

### `__init__(data_source, __ref=False, file_format=None)`

*   **Description:** Constructor for the `Ardal` class. Initializes an `Ardal` object from a data source.
*   **Parameters:**
    *   `data_source` (str or list): The path to the data file or a list containing the headers and matrix.
    *   `__ref` (bool, optional): A flag (currently unused). Defaults to `False`.
    *   `file_format` (str, optional): The format of the data file (e.g., "csv", "parquet", "npy", "json"). If `None`, the format is inferred from the file extension. Defaults to `None`.
*   **Raises:**
    *   `ValueError`: If parsing fails.

### `pairwise(guids=[], metric="hamming", chunk_size=100)`

*   **Description:** Calculates the pairwise distance matrix between the specified GUIDs.
*   **Parameters:**
    *   `guids` (list, optional): A list of GUIDs to compare. If empty, all GUIDs are compared. Defaults to `[]`.
    *   `metric` (str, optional): The distance metric to use ("hamming" or "jaccard"). Defaults to `"hamming"`.
    * `chunk_size` (int, optional): The chunk size to use when calculating the pairwise distance. Defaults to `100`.
*   **Returns:**
    *   `pd.DataFrame`: A Pandas DataFrame representing the distance matrix.
*   **Raises:**
    *   `ValueError`: If an invalid distance function is specified.

### `neighbourhood(guid, n, simd=True)`

*   **Description:** Gets the SNP neighborhood of a GUID.
*   **Parameters:**
    *   `guid` (str): The target GUID.
    *   `n` (int): The maximum Hamming distance.
    *   `simd` (bool, optional): Whether to use SIMD acceleration. Defaults to `True`.
*   **Returns:**
    *   `dict`: A dictionary where keys are neighboring GUIDs and values are their Hamming distances to the target GUID.
* **Raises:**
    * `ValueError`: If the guid is not found.

### `unique(guids)`

*   **Description:** Finds the set of SNPs unique to a given set of GUIDs. A SNP is considered unique if it is present in all of the specified GUIDs and absent in all other GUIDs.
*   **Parameters:**
    *   `guids` (list): A list of GUIDs.
*   **Returns:**
    *   `set`: A set of unique SNPs.
*   **Raises:**
    *   `ValueError`: If `guids` is not a list, if `guids` is empty, or if any GUID is not found.

### `core(guids, missingness=0.0, return_counts=False)`

*   **Description:** Takes a set of GUIDs and returns alleles common to this subset.
*   **Parameters:**
    *   `guids` (list): A list of GUIDs.
    *   `missingness` (float, optional): The maximum proportion of missingness allowed (between 0.0 and 1.0). Defaults to `0.0`.
    * `return_counts` (bool, optional): Whether to return a dictionary containing counts on the number of guids which exhibit this allele. Defaults to `False`.
*   **Returns:**
    *   `set`: A set of core alleles.
*   **Raises:**
    *   `ValueError`: If `guids` is not a list or set, if `guids` is empty, if any GUID is not found, or if `missingness` is not between 0 and 1.

### `accessory(guids, missingness=0.0, return_counts=False)`

*   **Description:** Takes a set of GUIDs and returns the accessory alleles (the symmetric set of the core alleles).
*   **Parameters:**
    *   `guids` (list): A list of GUIDs.
    *   `missingness` (float, optional): The maximum proportion of missingness allowed (between 0.0 and 1.0). Defaults to `0.0`.
    * `return_counts` (bool, optional): Whether to return a dictionary containing counts on the number of guids which exhibit this allele. Defaults to `False`.
*   **Returns:**
    *   `set`: A set of accessory alleles.
*   **Raises:**
    *   `ValueError`: If `guids` is not a list or set, if `guids` is empty, if any GUID is not found, or if `missingness` is not between 0 and 1.

### `allele(alleles)`

*   **Description:** Takes a set of alleles and returns all GUIDs containing all of those alleles.
*   **Parameters:**
    *   `alleles` (list): A list of alleles.
*   **Returns:**
    *   `set`: A set of GUIDs.
*   **Raises:**
    *   `ValueError`: If `alleles` is not a list or set, if `alleles` is empty, or if any allele is not found.

### `matchAlleleNames(expression)`

*   **Description:** Returns all allele names that match the given expression with wildcards.
*   **Parameters:**
    *   `expression` (str): A string containing the expression to match (e.g., "SNP*", "SNP1", "SNP[1-5]").
*   **Returns:**
    *   `set`: A set of matching allele names.
*   **Raises:**
    *   `ValueError`: If `expression` is not a string.

### `subsetbyGUID(guid_list)`

*   **Description:** Takes a list of GUIDs and subsets the allele matrix to include only these GUIDs, allowing for standard operations. Returns an Ardal object with the subsetted matrix.
*   **Parameters:**
    * `guid_list` (list): A list of GUIDs.
* **Returns:**
    * `None`: This function is not yet implemented.

### `toDict()`

*   **Description:** Returns a dictionary containing present allele IDs mapped to their GUID.
* **Returns:**
    * `None`: This function is not yet implemented.

### `stats()`

*   **Description:** Returns a dictionary containing information about the database and its size in memory.
*   **Returns:**
    *   `dict`: A dictionary containing statistics about the `Ardal` object, including:
        *   `n_guids`: The number of GUIDs.
        *   `n_alleles`: The number of alleles.
        *   `matrix_size`: The size of the allele matrix in memory (human-readable format).

### `getMatrix()`

*   **Description:** Returns the allele matrix.
*   **Returns:**
    *   `np.array`: The allele matrix as a NumPy array.

### `getHeaders()`

*   **Description:** Returns the headers.
*   **Returns:**
    *   `dict`: The headers as a dictionary.

### `snpCount()`

*   **Description:** Returns a dictionary of SNP counts for each GUID.
*   **Returns:**
    *   `dict`: A dictionary where keys are GUIDs and values are their SNP counts.

### `toDataFrame()`

*   **Description:** Returns the allele matrix as a Pandas DataFrame.
*   **Returns:**
    *   `pd.DataFrame`: The allele matrix as a Pandas DataFrame.

### `flushCache()`

*   **Description:** Flushes the distance cache.

## ArdalParser Class

The `ArdalParser` class is responsible for parsing allele matrix data from various file formats.

### `__init__(input_file_structure, file_format=None, prefix=None)`

*   **Description:** Constructor for the `ArdalParser` class.
*   **Parameters:**
    *   `input_file_structure` (str or list): The path to the data file or a list containing the headers and matrix.
    *   `file_format` (str, optional): The format of the data file (e.g., "csv", "parquet", "npy", "json"). If `None`, the format is inferred from the file extension. Defaults to `None`.
    *   `prefix` (str, optional): A prefix (currently unused). Defaults to `None`.
* **Raises:**
    * `ValueError`: If the file format is invalid, or if the input file structure is invalid.

## AlleleMatrix Class

The `AlleleMatrix` class is a C++ class (with Python bindings) that efficiently stores and manipulates the allele matrix data. It is primarily used internally by the `Ardal` class.

### `__init__(matrix)`

*   **Description:** Constructor for the `AlleleMatrix` class.
*   **Parameters:**
    *   `matrix` (`pybind11::array_t<uint8_t>`): A 2D NumPy array representing the allele matrix.
*   **Raises:**
    *   `std::runtime_error`: If the matrix is not 2D, if the matrix dimensions are too large, or if the matrix data type is not `uint8`.

### `access(coords)`

*   **Description:** Accesses elements in the matrix at specified coordinates.
*   **Parameters:**
    *   `coords` (`pybind11::array_t<size_t>`): A 2D NumPy array of coordinates (shape (k, 2)), where each row represents a (row, col) coordinate pair.
*   **Returns:**
    *   `pybind11::array_t<uint8_t>`: A 1D NumPy array containing the elements at the specified coordinates.
*   **Raises:**
    *   `std::runtime_error`: If the coordinates array is not 2D or does not have a shape of (k, 2), or if any coordinate is out of bounds.

### `accessGUID(guid_index)`

*   **Description:** Accesses the set of SNPs for a given GUID.
*   **Parameters:**
    *   `guid_index` (int): The index of the GUID (row) in the matrix.
*   **Returns:**
    *   `std::set<int>`: A set containing the indices of the SNPs present in the specified GUID.

### `getMatrix()`

*   **Description:** Returns the allele matrix.
*   **Returns:**
    *   `pybind11::array_t<uint8_t>`: A 2D NumPy array representing the allele matrix.

### `getMass()`

*   **Description:** Returns the mass of each row in the matrix.
*   **Returns:**
    *   `std::vector<int>`: A vector containing the mass of each row in the matrix.

### `hamming()`

*   **Description:** Calculates the Hamming distances between all pairs of rows.
*   **Returns:**
    *   `pybind11::array_t<int>`: A 1D NumPy array representing the condensed distance matrix containing the pairwise Hamming distances.

### `jaccard()`

*   **Description:** Calculates the Jaccard distances between all pairs of rows.
*   **Returns:**
    *   `pybind11::array_t<double>`: A 1D NumPy array representing the condensed distance matrix containing the pairwise Jaccard distances.

### `neighbourhood(row_coord, epsilon)`

*   **Description:** Finds the epsilon-neighborhood of a row using Hamming distance.
*   **Parameters:**
    *   `row_coord` (size_t): The index of the target row.
    *   `epsilon` (int): The maximum Hamming distance threshold.
*   **Returns:**
    *   `pybind11::array_t<int>`: A 1D NumPy array containing the indices of the rows that are within the epsilon-neighborhood of the target row.
*   **Raises:**
    *   `std::runtime_error`: If `row_coord` is out of range, or if `epsilon` is negative.

### `neighbourhoodSIMD(row_coord, epsilon)`

*   **Description:** Finds the epsilon-neighborhood of a row using Hamming distance (SIMD optimized).
*   **Parameters:**
    *   `row_coord` (size_t): The index of the target row.
    *   `epsilon` (int): The maximum Hamming distance threshold.
*   **Returns:**
    *   `pybind11::list`: A list of tuples, where each tuple contains the index of a neighboring row and its Hamming distance to the target row.
*   **Raises:**
    *   `std::runtime_error`: If `row_coord` is out of range, or if `epsilon` is negative.

### `gatherSNPs(guid_indices)`

*   **Description:** Gathers SNPs from multiple GUIDs.
*   **Parameters:**
    *   `guid_indices` (`const pybind11::array_t<int>`): A NumPy array containing the indices of the GUIDs.
*   **Returns:**
    *   `std::vector<int>`: A vector containing the indices of all SNPs present in the specified GUIDs.

### `getMass()`

*   **Description:** Gets the mass of each row in the matrix.
*   **Returns:**
    *   `std::vector<int>`: A vector containing the mass of each row in the matrix.

### `flushCache()`

*   **Description:** Clears (flushes) the Hamming distance cache.
