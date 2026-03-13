# Ardal API Documentation

The Ardal framework provides modular, high-performance operations for binary SNP matrices and their associated metadata. The API is organised into a lightweight Python frontend with pluggable component classes backed by C++ bindings exposed through `pybind11`. This document reflects the current (0.3.x) Python surface.

---

## `Ardal`
Main user-facing entry-point. Wraps the `_ardal.HybridMatrix` backend and exposes component classes for I/O, querying, statistics, and distance calculations.

### Constructor
```python
Ardal(
    data_source: Union[
        str,
        Tuple[str, str],
        Tuple[np.ndarray, Mapping[str, Sequence[str]]],
        list
    ],
    roaring: Union[str, bool] = "auto",
    density_threshold: float = 0.1,
    allele_id_format: Optional[str] = None,
    allele_coords_bed: Optional[str] = None,
    verbosity: Union[int, str] = "error",
    quiet_init: bool = False,
    allele_positions: Optional[Dict[str, Tuple[str, int]]] = None,
)
```
- `data_source`: single path, pair of paths, or in-memory `[matrix, headers]`/`[headers, matrix]`.
- `roaring`: `"auto"`, `"true"`, `"false"`, or boolean to control roaring backend provisioning.
- `density_threshold`: heuristic for auto backend selection when `roaring="auto"`.
- `allele_id_format` / `allele_coords_bed`: optional position decoding helpers passed to the header utilities.
- `verbosity`: logging level (`"debug"`, `"info"`, `"warn"`, `"error"`, `"critical"`, `"silent"` or logging constant).
- `quiet_init`: suppress the initial `matrix_stats` table.
- `allele_positions`: optional pre-computed allele position map.

### Attributes
- `.io`: `ArdalIO` instance for export/import helpers.
- `.get`: `ArdalGet` instance for subsetting, metadata, and matrix materialisation.
- `.allele`: `ArdalAllele` instance for allele-centric queries.
- `.distance`: `ArdalDistance` instance for pairwise distances, neighbourhoods, and k-NN.
- `.stats`: `ArdalStats` instance for statistical summaries.
- `.roaring`: boolean indicating whether a roaring backend was initialised.

### Methods
- `set_verbosity(level_param: Union[int, str]) -> None`: adjust logging level for the root logger and subordinate components.

---

## `ArdalAllele`
Allele-centric utilities that operate on the binary matrix and decoded headers.

### Methods
- `unique(guids: List[str], backend: str = "auto") -> Dict[str, Set[str]]`  
  Returns alleles present in each provided GUID and absent from all others. `backend` accepts `"auto"`, `"bit"`, or `"roaring"`.

- `unique_core(guids: List[str], backend: str = "auto") -> List[str]`  
  Alleles present in **all** provided GUIDs and absent elsewhere. Logs a warning if none are found.

- `guids_with_alleles(alleles: List[str], backend: str = "auto") -> List[str]`  
  GUIDs containing every allele in `alleles`.

- `match_names(expression: str) -> List[str]`  
  Wildcard (`*`) expression match against allele IDs.

- `interval(...) -> List[str]`  
  Fetch alleles that fall within user-supplied genomic intervals or BED coordinates. Supports optional `intervals`, `intervals_bed`, `allele_coords_bed`, and `allele_id_format`.

- `count() -> Dict[str, int]`  
  Per-GUID allele counts (row sums).

- `get_positions() -> Dict[str, Tuple[str, int]]`  
  Mapping of allele ID → `(chromosome, position)` derived from header metadata or supplied lookup.

---

## `ArdalDistance`
Pairwise distance matrices, epsilon-neighbourhoods, and k-nearest neighbour queries.

### Methods
- `pairwise(...) -> Union[np.ndarray, pd.DataFrame]`  
  ```python
  pairwise(
      guids: Optional[List[str]] = None,
      alleles: Optional[List[str]] = None,
      intervals: Optional[List] = None,
      intervals_bed: Optional[str] = None,
      allele_coords_bed: Optional[str] = None,
      metric: str = "hamming",
      use_simd: bool = True,
      threads: int = 1,
      backend: str = "auto",
      allele_id_format: str = "{chr}.{start}.{ref}.{alt}",
      return_square: bool = False,
      as_dataframe: bool = False,
  )
  ```
  Produces a condensed 1D distance vector by default, or an `n × n` array/DataFrame when `return_square=True`. Supported metrics: `"hamming"`, `"jaccard"`, `"inner_product"`, `"cosine"`. Interval filters allow region-restricted computations. Raises `ParameterError` or `MemoryError` for invalid requests.

- `snv_neighbourhood(guid: str, n: int) -> Dict[str, int]`  
  Experimental positional neighbourhood based on decoded allele IDs. Assumes coordinate-encoded IDs and is not production ready.

- `neighbourhood(guid: str, n: int, use_simd: bool = True, threads: int = 1, backend: str = "auto", metric: str = "hamming") -> Dict[str, int]`  
  Hamming or inner-product neighbourhood search up to radius `n`. Returns `{guid: distance}`.

- `knn(guid: str, k: int, use_simd: bool = True, threads: int = 1, backend: str = "auto", metric: str = "hamming") -> List[Tuple[str, float]]`  
  k-nearest neighbours for the requested metric (`"hamming"`, `"inner_product"`, `"jaccard"`, `"cosine"`). Returns a list sorted by similarity/distance with per-neighbour scores.

---

## `ArdalGet`
Matrix subsetting, metadata inspection, and conversion helpers.

### Methods
- `subset(...) -> Union["Ardal", List[Any]]`  
  ```python
  subset(
      guids: List[str] = [],
      alleles: List[str] = [],
      data_only: bool = False,
      threads: int = 1,
      child_verbosity: str = "silent",
      child_quiet_init: bool = True,
  )
  ```
  Subset by GUID and/or allele. Returns a new `Ardal` instance unless `data_only=True`, where it returns `[matrix, headers]`. Accepts packed backends and can materialise dense rows when only GUIDs are filtered.

- `matrix_stats(print_table: bool = False) -> Dict[str, Any]`  
  Summary of counts, density, and byte sizes. Optionally prints a formatted table.

- `density() -> float`  
  Matrix density (fraction of set bits).

- `bit_matrix() -> np.ndarray`  
  Dense uint8 representation of the bit matrix.

- `packed_matrix() -> np.ndarray`  
  `<u8` bit-packed matrix (words-per-row × GUIDs).

- `roaring_matrix(decode: bool = True) -> Union[List[np.ndarray], Dict[str, List[str]]]`  
  Roaring bitmap view; optionally decoded to allele IDs. Raises `RoaringError` if the roaring backend was not initialised.

- `headers() -> Dict[str, List[str]]`  
  Access the GUID and allele headers.

---

## `ArdalIO`
I/O and serialisation helpers.

### Methods
- `to_dataframe() -> pd.DataFrame`  
  Dense DataFrame with GUID rows and allele columns.

- `to_dict(backend: str = "auto") -> Dict[str, List[str]]`  
  Dictionary mapping GUID → allele list. Supports `"auto"`, `"bit"`, and `"roaring"` backend selection.

- `write(output_prefix: str, out_directory: str = "./", format: str = "bin") -> None`  
  Persist the matrix plus headers. `format` supports `"npy"`, `"npz"`, and `"bin"` (bit-packed). Raises `MatrixWriteError` and `ParameterError` for invalid destinations.

- `make_fastas(...)`  
  Placeholder for FASTA generation; currently raises `NotImplementedError`.

---

## `ArdalStats`
Statistical summaries and divergence metrics.

### Methods
- `allele_frequency(guids: List[str] = []) -> Dict[str, float]`  
  Allele frequencies within the supplied GUID subset (defaults to all GUIDs).

- `allele_entropy() -> Dict[str, float]`  
  Shannon entropy for every allele.

- `allele_missingness(guids: List[str] = [], alleles: List[str] = [], normalise: bool = False, window_size: Optional[int] = None, window_step: Optional[int] = None) -> np.ndarray`  
  Per-allele missing counts across selected GUIDs. If `normalise=True`, returns proportions (count / n_guids). If `window_size` is provided, returns window-averaged missingness across the allele order.

- `allele_cooc(alleles: List[str] = [], threshold: float = 0.95, threads: int = 1) -> Dict[str, List[str]]`  
  Allele co-occurrence above `threshold`. Subset to specific alleles or compute across the entire matrix.

- `allele_inform(guids: List[str], method: str = "kullbackleibler") -> Dict[str, float]`  
  Divergence-style scores contrasting in-group vs out-group (`"kullbackleibler"`, `"jensenshannon"`, `"informationgain"`).

- `test_associations(guids: List[str], tests: Optional[List[str]] = None) -> Dict[str, Dict[str, Dict[str, float]]]`  
  Wrapper for chi-squared and Fisher exact tests on 2×2 allele contingency tables. Returns nested dictionaries keyed by allele → test → statistic/p-value.

---

## `ArdalParser`
Input loader responsible for normalising dense or bit-packed sources.

### Constructor
```python
ArdalParser(
    input_data_structure: Union[list, str],
    verify_hash: bool = False,
    is_packed_mem: bool = False,
)
```

### Attributes
- `matrix`: Dense `np.ndarray` or `<u8` memory map ready for the backend.
- `headers`: `{"guids": [...], "alleles": [...]}` dictionary.
- `meta`: Bit-pack metadata when applicable (row/column counts, strides).
- `file_format`: Detected source format (`"csv"`, `"parquet"`, `"npy"`, `"npz"`, `"bitpack"`, etc.).
- `is_bitpacked`: Boolean flag indicating whether the matrix is already packed.

Supported inputs include CSV, Parquet, NPY/JSON, NPZ/JSON, bit-packed `{headers.json, matrix.bin}` pairs, and in-memory `[matrix, headers]` combinations.

---

## Backend Selection Helpers
Most high-level methods accept a `backend` argument with values:
- `"auto"` (default): Heuristic selection using matrix density and build-time roaring availability.
- `"bit"`: Force the bit-packed backend.
- `"roaring"`: Force the roaring bitmap backend (raises if unavailable).

SIMD acceleration can be toggled via `use_simd` flags on distance-related calls; thread parallelism is exposed through `threads` parameters where supported.

---

For additional performance tuning, consult the C++ backend headers (`BitMatrix`, `RoaringMatrix`, `HybridMatrix`) and the utilities under `ardal/utils`.
