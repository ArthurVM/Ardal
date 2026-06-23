# Ardal TODOs

## Matrix formats

- Keep `ardal.bin.v2` as the default active-analysis format: uncompressed, sectioned, and memmappable.
- Treat `.zst` as archive/transfer format unless section-aware compression is implemented.
- Design a future section-compressed format where one file can contain:
  - compressed allele matrix section
  - uncompressed missing offsets section
  - uncompressed missing ranges section
- Make section metadata explicit about `compression`, `offset`, `nbytes`, `uncompressed_nbytes`, `dtype`, and `shape`.

## Missingness

- Keep lazy missing interval loading as the default.
- Add an explicit `lazy_missing=False` path to inflate all missingness at matrix load time.
- Keep missingness in binary sections to avoid metadata bloom.
- Continue optimising missing interval use in the C++ backends, especially direct range-to-roaring construction.

## Builder API

- Keep `misc/*_v2.py` as reference scripts until the module API is fully settled.
- Maintain strict builder module responsibilities:
  - `schemas.py`: sidecar schemas
  - `vcf.py`: VCF parsing and filtering
  - `annotation.py`: FASTA/GFF annotation
  - `sample_json.py`: sample JSON handling
  - `missingness.py`: missing interval handling
  - `matrix_plan.py`: indexing and site statistics
  - `writers.py`: matrix and metadata writing
  - `api.py`: orchestration
- Keep refusing `genic_snps` / `nonsyns` builds when JSON inputs lack the required annotations.
- Build any future CLI on top of the Python API, not inside it.

## Metadata

- Keep `allele_model` in sidecars.
- Add provenance fields for reference IDs/checksums and VCF filter settings.
- Keep large data out of metadata; sidecars should describe binary sections, not contain them.

## Distance correctness

- Keep SNV-specific methods restricted to nucleotide allele models.
- Preserve multi-ALT SNV handling.
- Maintain regression tests for pairwise SNV, KNN SNV, neighbourhood SNV, multi-ALT sites, and protein-matrix refusal.

## Performance

- Prefer direct `.bin v2` generation over dense `.npy` / `.npz` staging.
- Keep `.npy` / `.npz` only as redundant compatibility outputs for now.
- Benchmark load time, peak RSS, first-query latency, KNN/neighbourhood runtime, and missingness-heavy datasets.
- Investigate vectorised or C++ direct bitpacking in the builder path.

## Checkpoints

- Keep focused tests for:
  - `.bin v2` parser loading
  - lazy binary missing sections
  - builder API
  - annotated projection refusal
  - SNV correctness
  - amino-acid matrix refusal for SNV methods
- Every matrix schema emitted by `ardal.builder` must have a parser round-trip test in the same change.
