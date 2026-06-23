import json

import pytest

from ardal.builder import json_directory_to_matrices, vcf_to_sample_json
from ardal.core.ArdalParser import ArdalParser


def assert_builder_matrix_loads(matrix_path):
    parser = ArdalParser([str(matrix_path), str(matrix_path) + ".meta"])
    assert parser.is_bitpacked
    assert parser.headers["guids"]
    assert parser.headers["alleles"]
    assert parser.matrix.shape[0] == len(parser.headers["guids"])


def assert_builder_dense_matrix_loads(matrix_path):
    parser = ArdalParser([str(matrix_path), str(matrix_path) + ".meta"])
    assert not parser.is_bitpacked
    assert parser.headers["guids"]
    assert parser.headers["alleles"]
    assert parser.matrix.shape[0] == len(parser.headers["guids"])


def test_builder_duty_modules_expose_separate_responsibilities():
    from ardal.builder.annotation import VariantAnnotator
    from ardal.builder.matrix_plan import index_samples_and_alleles
    from ardal.builder.missingness import missing_blocks_from_bed
    from ardal.builder.sample_json import sample_payload
    from ardal.builder.vcf import sample_alleles_from_vcf
    from ardal.builder.writers import write_matrix_meta

    assert VariantAnnotator.__name__ == "VariantAnnotator"
    assert index_samples_and_alleles.__name__ == "index_samples_and_alleles"
    assert missing_blocks_from_bed.__name__ == "missing_blocks_from_bed"
    assert sample_payload.__name__ == "sample_payload"
    assert sample_alleles_from_vcf.__name__ == "sample_alleles_from_vcf"
    assert write_matrix_meta.__name__ == "write_matrix_meta"


def test_vcf_to_sample_json_writes_schema_v2_payload(tmp_path):
    vcf_path = tmp_path / "sample.vcf"
    depth_path = tmp_path / "sample.bed"
    outdir = tmp_path / "json"

    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr1\t2\t.\tA\tC\t60\tPASS\tDP=20\tGT\t1\n"
    )
    depth_path.write_text("chr1\t1\t9\nchr1\t2\t20\nchr1\t3\t0\n")

    vcf_to_sample_json(
        vcf_in=vcf_path,
        depth_in=depth_path,
        outdir=outdir,
        overwrite=True,
    )

    payload = json.loads((outdir / "s1.json").read_text())
    assert payload["schema_version"] == "ardal.sample_variants.v2"
    assert payload["alleles"] == [["chr1", 2, 2, "A", "C"]]
    assert payload["missing"] == [["chr1", 1, 1], ["chr1", 3, 3]]


def test_json_directory_to_matrices_writes_bin_v2_sidecar(tmp_path):
    json_dir = tmp_path / "json"
    json_dir.mkdir()
    output_prefix = tmp_path / "matrix"

    samples = {
        "s1": {
            "schema_version": "ardal.sample_variants.v2",
            "sample_id": "s1",
            "alleles": [["chr1", 1, 1, "A", "C"], ["chr1", 2, 2, "G", "T"]],
            "genic": [],
            "nonsynonymous": [],
            "missing": [["chr1", 5, 7]],
        },
        "s2": {
            "schema_version": "ardal.sample_variants.v2",
            "sample_id": "s2",
            "alleles": [["chr1", 1, 1, "A", "G"]],
            "genic": [],
            "nonsynonymous": [],
            "missing": [["chr1", 2, 2]],
        },
    }
    for sample_id, payload in samples.items():
        (json_dir / f"{sample_id}.json").write_text(json.dumps(payload))

    json_directory_to_matrices(
        json_dir,
        str(output_prefix),
        emit_agnostic=False,
        data_types=("snps",),
        json_indent=2,
    )

    meta = json.loads((tmp_path / "matrix.snps.ref.bin.meta").read_text())
    assert meta["meta"]["format"] == "ardal.bin.v2"
    assert meta["meta"]["allele_model"]["alphabet"] == "nucleotide"
    assert meta["missing"]["encoding"] == "binary_sections"
    assert "sections" in meta["meta"]
    assert_builder_matrix_loads(tmp_path / "matrix.snps.ref.bin")


def test_builder_dense_and_dense_path_bin_outputs_are_parser_supported(tmp_path):
    json_dir = tmp_path / "json"
    json_dir.mkdir()
    output_prefix = tmp_path / "matrix"

    samples = {
        "s1": {
            "schema_version": "ardal.sample_variants.v2",
            "sample_id": "s1",
            "alleles": [["chr1", 1, 1, "A", "C"]],
            "genic": [],
            "nonsynonymous": [],
            "missing": [],
        },
        "s2": {
            "schema_version": "ardal.sample_variants.v2",
            "sample_id": "s2",
            "alleles": [["chr1", 2, 2, "G", "T"]],
            "genic": [],
            "nonsynonymous": [],
            "missing": [],
        },
    }
    for sample_id, payload in samples.items():
        (json_dir / f"{sample_id}.json").write_text(json.dumps(payload))

    json_directory_to_matrices(
        json_dir,
        str(output_prefix),
        emit_npy=True,
        emit_npz=True,
        emit_agnostic=False,
        data_types=("snps",),
        json_indent=2,
    )

    assert_builder_dense_matrix_loads(tmp_path / "matrix.snps.ref.npy")
    assert_builder_dense_matrix_loads(tmp_path / "matrix.snps.ref.npz")
    assert_builder_matrix_loads(tmp_path / "matrix.snps.ref.bin")


def test_json_directory_to_matrices_refuses_unannotated_genic_requests(tmp_path):
    json_dir = tmp_path / "json"
    json_dir.mkdir()
    output_prefix = tmp_path / "matrix"

    payload = {
        "schema_version": "ardal.sample_variants.v2",
        "sample_id": "s1",
        "alleles": [["chr1", 1, 1, "A", "C"]],
        "genic": [],
        "nonsynonymous": [],
        "missing": [],
    }
    (json_dir / "s1.json").write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="Requested genic_snps"):
        json_directory_to_matrices(
            json_dir,
            str(output_prefix),
            emit_agnostic=False,
            data_types=("genic_snps",),
            json_indent=2,
        )

    assert not list(tmp_path.glob("matrix.*"))


def test_json_directory_to_matrices_refuses_unannotated_nonsyn_requests(tmp_path):
    json_dir = tmp_path / "json"
    json_dir.mkdir()
    output_prefix = tmp_path / "matrix"

    payload = {
        "schema_version": "ardal.sample_variants.v2",
        "sample_id": "s1",
        "alleles": [["chr1", 1, 1, "A", "C"]],
        "genic": [0],
        "nonsynonymous": [],
        "missing": [],
    }
    (json_dir / "s1.json").write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="Requested nonsyns"):
        json_directory_to_matrices(
            json_dir,
            str(output_prefix),
            emit_agnostic=False,
            data_types=("nonsyns",),
            json_indent=2,
        )

    assert not list(tmp_path.glob("matrix.*"))


def test_json_directory_to_matrices_builds_when_requested_annotations_exist(tmp_path):
    json_dir = tmp_path / "json"
    json_dir.mkdir()
    output_prefix = tmp_path / "matrix"

    payload = {
        "schema_version": "ardal.sample_variants.v2",
        "sample_id": "s1",
        "alleles": [["chr1", 1, 1, "A", "C"]],
        "genic": [0],
        "nonsynonymous": [[0, "geneA", 1, "K", "N"]],
        "missing": [],
    }
    (json_dir / "s1.json").write_text(json.dumps(payload))

    json_directory_to_matrices(
        json_dir,
        str(output_prefix),
        emit_agnostic=False,
        data_types=("genic_snps", "nonsyns"),
        json_indent=2,
    )

    assert (tmp_path / "matrix.genic_snps.ref.bin.meta").exists()
    assert (tmp_path / "matrix.nonsyns.ref.bin.meta").exists()
    assert_builder_matrix_loads(tmp_path / "matrix.genic_snps.ref.bin")
    assert_builder_matrix_loads(tmp_path / "matrix.nonsyns.ref.bin")
