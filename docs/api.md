# Ardal API Documentation

The Ardal framework provides modular, high-performance operations for binary SNP matrices and associated metadata in pathogen genomics. The following API documentation outlines the purpose and usage of each core module.

---

## `Ardal`
**Purpose**: Main user-facing interface to the Ardal API. Constructs a wrapped `HybridMatrix` and exposes submodules for I/O, querying, statistics, and comparison.

### Constructor:
```python
Ardal(data_source: str,
      use_roaring_if_sparse: bool = True,
      density_threshold: float = 0.02,
      force_roaring: bool = False,
      __ref: bool = False,
      file_format: str = None)
```

### Attributes:
- `.io`: Instance of `ArdalIO` for data export/import
- `.get`: Instance of `ArdalGet` for matrix subsetting and introspection
- `.allele`: Instance of `ArdalAllele` for SNP querying
- `.compare`: Instance of `ArdalCompare` for pairwise distances and neighbourhoods
- `.stats`: Instance of `ArdalStats` for statistical analysis

**Notes**:
- Internally uses `HybridMatrix` from the `_ardal` C++ backend.
- Automatically detects and constructs Roaring or bitpacked representations.
- Performs an initial `matrixStats()` printout for user feedback.

---

## `ArdalAllele`
**Purpose**: Provides allele-centric queries on the binary SNP matrix.

### Methods:

#### `unique(guids: list, force_bit_backend: bool = False) -> dict`
Find SNPs unique to each GUID (present in one, absent in all others).
- **Input**: list of GUIDs
- **Returns**: `{guid: set of unique SNPs}`
- **Raises**: `ValueError` on invalid input

#### `guidsWithAlleles(alleles: list, force_bit_backend: bool = False) -> set`
Find GUIDs that contain *all* specified alleles.
- **Input**: list/set of allele IDs
- **Returns**: set of GUIDs

#### `matchNames(expression: str) -> set`
Wildcard pattern match on allele names.
- `*` is used for globbing

---

## `ArdalCompare`
**Purpose**: Distance and neighbourhood computation.

### Methods:

#### `pairwise(metric='hamming', use_simd=True, threads=1, force_bit_backend=False) -> pd.DataFrame`
Calculate all-vs-all distance matrix.
- **Metrics**: `'hamming'`, `'jaccard'`, `'inner_product'`
- **Returns**: `pandas.DataFrame`
- **Raises**: `MemoryError`, `ValueError`

#### `snvNeighbourhood(guid: str, n: int) -> dict`
Find GUIDs within `n` differing SNP positions (by parsed string position).
- **Warning**: Not production ready, assumes `{ref}{pos}{alt}` allele ID format

#### `neighbourhood(guid: str, n: int, use_simd=True, threads=1, force_bit_backend=False) -> dict`
Hamming-distance neighbourhood within `n` allelic differences.

---

## `ArdalGet`
**Purpose**: Matrix subsetting, conversion, and summarisation.

### Methods:

#### `subset(guid_list=[], allele_list=[]) -> list[np.ndarray, dict]`
Subset matrix by GUID and/or allele.
- **Returns**: `[new_matrix, new_headers]`
- **Raises**: `ParameterError`, `InvalidGUIDQueryError`, `InvalidAlleleQueryError`

#### `matrixStats(print_stats=False) -> dict`
Return metadata and memory stats for the matrix.

#### `density() -> float`
Return sparsity ratio of the matrix.

#### `bitMatrix() -> np.ndarray`
Return matrix as dense bit array.

#### `roaringMatrix(decode=True) -> dict`
Return roaring bitmap matrix.
- **Raises**: `RoaringError` if not initialised with roaring

#### `headers() -> dict`
Return headers (metadata) dictionary.

#### `snpCount() -> dict`
Return SNP count per GUID.

---

## `ArdalIO`
**Purpose**: Matrix I/O, conversion to formats.

### Methods:

#### `toDataFrame() -> pd.DataFrame`
Convert matrix to Pandas DataFrame.

#### `toDict(force_bit_backend=False) -> dict`
Return dictionary `{guid: list of alleles}`

#### `write(file_path, output_prefix, npz=False)`
Save matrix to disk (as `.npy` or `.npz` + headers).
- **Raises**: `MatrixWriteError`

#### `makeFastas(guids=[], ref=None, allele_id_format="{ref}.{chr}.{start}.{alt}")`
Experimental. Generate simulated FASTA per sample using SNPs.

---

## `ArdalStats`
**Purpose**: Statistical measures on SNP distribution.

### Methods:

#### `af(guids=[]) -> dict`
Compute allele frequency across provided GUIDs (or all).

#### `alleleEntropy() -> dict`
Shannon entropy for each allele.

#### `alleleCooc(allele_indices=[], threshold=0.95, threads=1) -> dict`
Return alleles co-occurring above threshold with specified alleles (or all).

#### `snpInform(guids, metric='kullbackleibler') -> dict`
Score SNPs by in-group vs out-group discriminative power.
- **Metrics**: `'kullbackleibler'`, `'jensenshannon'`, `'informationgain'`

#### `testSnpAssociations(...)`
Stub for future statistical tests (chi2, Fisher’s exact).

---

## `ArdalParser`
**Purpose**: Load matrix + headers from disk or memory.

### Constructor:
```python
ArdalParser(input_data_structure, file_format=None, prefix=None)
```
- Supports CSV, Parquet, NPY/JSON pair, or in-memory `[np.ndarray, dict]`

### Attributes:
- `matrix`: ndarray allele matrix
- `headers`: dict with keys `guids`, `alleles`

---

## Utility Functions (from `utilities.py`)

- `checkGUIDs(guids, headers, filter=False)`
- `checkAlleles(alleles, headers, filter=False)`
- `encodeGuid(guid, headers) -> int`
- `decodeGuid(index, headers) -> str`
- `encodeAllele(allele, headers) -> int`
- `decodeAllele(index, headers) -> str`
- `decodeAlleleID(allele_id, allele_id_format) -> tuple`

---

## Notes
- All methods assume headers contain correct `guids` and `alleles` keys.
- BitMatrix backend is used by default unless RoaringMatrix is enabled and appropriate.
- Roaring matrix must be explicitly enabled at object creation.

---

For advanced usage and performance tuning (e.g. SIMD, multi-threading, cache tuning), refer to the `BitMatrix` and `HybridMatrix` C++ bindings exposed in the backend.
