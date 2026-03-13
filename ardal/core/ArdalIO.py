""" ArdalIO.py
This module provides functionality for reading and writing allele matrices in the Ardal framework.
"""
import os
import pathlib
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Union
from importlib.metadata import version, PackageNotFoundError

from ..utils.misc import require_package
from ..utils.make_meta import make_meta
from ..utils.decorators import check_backend_argument
from ..utils.exceptions import MatrixWriteError, ParameterError
from ..utils.logger import get_logger

log = get_logger()


## core/ArdalIO.py
class ArdalIO:
    """ Class for reading and writing allele matrices in the Ardal framework.
    """

    def __init__( self,
                  headerUtils,
                  hybrid_matrix,
                  roaring_enabled : bool ):

        self._headerUtils = headerUtils
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled

    def _build_column_masks_payload(self) -> dict:
        """Return canonical per-GUID missing-column masks."""
        payload = {}
        for guid in self._headerUtils.headers.get("guids", []):
            payload[guid] = list(self._headerUtils.get_guid_missing_mask(guid))
        return payload


    def _normalise_meta_for_write(self, meta: dict, matrix_file: str) -> dict:
        """Emit metadata with stable key names expected by Ardal write outputs."""
        return {
            "format": meta.get("format"),
            "dtype": meta.get("dtype"),
            "endianness": meta.get("endianness"),
            "row_major": meta.get("row_major"),
            "n_rows": meta.get("n_rows"),
            "n_cols": meta.get("n_cols"),
            "matrix_file": matrix_file,
            "data_nbytes": meta.get("data_nbytes"),
            "data_sha256": meta.get("data_sha256"),
            "words_per_row": meta.get("words_per_row"),
            "bits_per_word": meta.get("bits_per_word"),
            "row_stride_bytes": meta.get("row_stride_bytes"),
            "generated_by": meta.get("generated_by"),
        }


    def _build_headers_meta_payload(self, meta: dict, matrix_file: str) -> dict:
        return {
            "meta": self._normalise_meta_for_write(meta, matrix_file),
            "headers": self._headerUtils.headers,
            "column_masks": self._build_column_masks_payload(),
        }


    def to_dataframe( self ) -> pd.DataFrame:
        """ Return the allele matrix as a Pandas DataFrame.
        """
        return pd.DataFrame(self._matrix.getBitMatrix(), index=self._headerUtils.headers["guids"], columns=self._headerUtils.headers["alleles"])


    @check_backend_argument
    def to_dict( self,
                 backend : str = "auto"
                 ) -> dict:
        """ Return a dictionary containing present allele IDs mapped to their guid.
        """
        allele_dict = defaultdict(list)
        for guid_idx, guid_name in enumerate(self._headerUtils.headers["guids"]):
            allele_indices = self._matrix.getSetBitIndices(guid_idx, backend=backend)
            for allele_idx in allele_indices:
                allele_dict[guid_name].append(self._headerUtils.decode_allele(allele_idx))
        return dict(allele_dict)
    

    def write( self,
               output_prefix : str,
               out_directory : str = "./",
               format : str = "bin",
               ) -> None:
        """ 
        Write the allele matrix to disk.
        Writes as a numpy/JSON pair.
        The format flag defines what 
        """
        """ Write the allele matrix to disk as a matrix/headers pair in specified format.
        Format must be one of:
        `npy` : outputs a dense uncompressed .npy matrix and .json headers.
        `npz` : outputs a dense compressed .npz matrix and .json headers.
        `bin` : outputs a bitpacked uncompressed .bin matrix with a .json containing headers and packing metadata.

        Args:
            file_path (str): _description_
            headers (str): _description_
            format (str) : _description_

        Returns:
            None
        """
        json = require_package("json", "json")
        
        ALLOWED_FORMATS = ["npy", "npz", "bin"]

        if format not in ALLOWED_FORMATS:
            raise ParameterError(f"Invalid matrix output format : {format}. Must be one of {ALLOWED_FORMATS}.")
        
        if not os.path.exists(out_directory):
            raise ParameterError(f"Path '{out_directory}' does not exist.")
        
        ## generate file paths
        matrix_out_path = os.path.join(out_directory, output_prefix + "." + format)
        headers_out_path = os.path.join(matrix_out_path + ".meta")
        
        if os.path.exists(headers_out_path):
            raise MatrixWriteError(f"File '{headers_out_path}' already exists.")
        if os.path.exists(matrix_out_path):
            raise MatrixWriteError(f"File '{matrix_out_path}' already exists.")

        matrix_to_save = self._matrix.getBitMatrix()
        if format == "npy":
            log.info("Writing dense matrix as .npy")
            np.save(matrix_out_path, matrix_to_save)
            meta = make_meta(matrix_to_save,
                             self._headerUtils.headers,
                             generated_by="ardal::io::write",
                             format_name="ardal.dense.v1",
                             matrix_file=matrix_out_path)
            headers_meta = self._build_headers_meta_payload(meta, matrix_out_path)
        elif format == "npz":
            log.info("Writing dense matrix as .npz")
            np.savez_compressed(matrix_out_path, matrix=matrix_to_save)
            meta = make_meta(matrix_to_save,
                             self._headerUtils.headers,
                             generated_by="ardal::io::write",
                             format_name="ardal.dense.v1",
                             matrix_file=matrix_out_path)
            headers_meta = self._build_headers_meta_payload(meta, matrix_out_path)
        elif format == "bin":
            log.info("Writing packed matrix as .bin")
            headers_meta = self._write_packed(matrix_out_path)
            
        log.info(f"Wrote allele matrix to disk : {matrix_out_path}")

        with open(headers_out_path, 'w') as fout:
            json.dump(headers_meta, fout, indent=4)
            
        log.info(f"Wrote headers to disk : {headers_out_path}")
        
        
    def _write_packed( self,
                       matrix_out_path : str ):
        arr = self._matrix.getPackedMatrix()
        
        ## ensure little endian
        if arr.dtype.byteorder not in ('<', '='):
            arr = arr.byteswap().newbyteorder('<')
        
        ## write bin
        bin_path = pathlib.Path(matrix_out_path)
        arr.tofile(bin_path)  ## raw little-endian uint64 words

        meta = make_meta(arr,
                         self._headerUtils.headers,
                         generated_by="ardal::io::write",
                         format_name="ardal.bitpack.v1",
                         matrix_file=matrix_out_path)
        headers_meta = self._build_headers_meta_payload(meta, matrix_out_path)
        
        return headers_meta
        
    
    @staticmethod
    def _generated_by_string() -> str:
        return "ardal::io::write"
    

    def make_fastas( self,
                     guids : list = [],
                     ref : Union[str, None] = None,
                     allele_id_format : Union[str, None] = None,
                     out_directory : str = "./",
                     snp_only : bool = True
                     ) -> dict:
        """ Takes a set of guids and a reference fasta and constructs a simulated fasta from each dataset using alleles stored
        in the allele matrix. The position of each allele is determined by the allele_id_format argument, whereby the keywords:
        ref, alt, chr, start, end, are used to outline the allele id naming convention.

        E.g. for the allele `A.chr1.100.101.T` would be decoded using the allele_id_format string `{ref}.{chr}.{start}.{end}.{alt}`.

        Returns:
            dict: Mapping of guid to output FASTA file paths.
        """
        # SeqIO = require_package("Bio", attr="SeqIO")
        # Seq = require_package("Bio.Seq", import_as="Bio.Seq", attr="Seq")
        # SeqRecord = require_package("Bio.SeqRecord", import_as="Bio.SeqRecord", attr="SeqRecord")
        
        from Bio import SeqIO, Seq, SeqRecord

        if ref is None:
            raise ParameterError("Reference FASTA path is required.")
        if not os.path.exists(ref):
            raise ParameterError(f"Reference FASTA not found: {ref}")

        if not out_directory:
            raise ParameterError("Output directory is required.")
        if not os.path.exists(out_directory):
            raise ParameterError(f"Output directory not found: {out_directory}")
        if not os.path.isdir(out_directory):
            raise ParameterError(f"Output directory is not a directory: {out_directory}")

        if guids:
            self._headerUtils.check_guids(guids)
        else:
            guids = list(self._headerUtils.headers.get("guids", []))

        if not allele_id_format:
            allele_id_format = self._headerUtils.get_cached_allele_id_format()
        if not allele_id_format:
            raise ParameterError("allele_id_format is required to decode allele positions.")

        pattern = self._headerUtils._check_allele_format_grammar(allele_id_format=allele_id_format)

        try:
            ref_records = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
        except Exception as exc:
            raise MatrixWriteError(f"Failed to parse reference FASTA: {exc}") from exc

        if not ref_records:
            raise MatrixWriteError("Reference FASTA appears to be empty.")

        ref_seqs: dict[str, str] = {
            name: str(record.seq).upper()
            for name, record in ref_records.items()
        }

        try:
            ardal_version = version("ardal")
        except PackageNotFoundError:
            ardal_version = None
        generated_by = "ardal::io::make_fastas"

        ## resolve allele positions and decode alleles once
        self._headerUtils.ensure_id_positions(allele_id_format)
        allele_position_map = self._headerUtils.get_allele_positions()
        alleles = self._headerUtils.headers.get("alleles", [])

        ref_default_chr = None
        if len(ref_seqs) == 1:
            ref_default_chr = next(iter(ref_seqs.keys()))

        ## precompute allele edits (chr, pos0, span_len, alt)
        allele_info = [None] * len(alleles)
        for idx, allele_id in enumerate(alleles):
            coord = allele_position_map.get(allele_id)
            chr_key = coord[0] if coord else None
            start = coord[1] if coord else None

            try:
                dec_chr, dec_start, dec_end, ref_base, alt_base = self._headerUtils._decode_allele_position(
                    allele_id=allele_id,
                    pattern=pattern,
                    allele_id_format=allele_id_format
                )
            except ValueError:
                dec_chr = None
                dec_start = None
                dec_end = None
                ref_base = None
                alt_base = None

            if chr_key is None:
                chr_key = dec_chr
            if start is None:
                start = dec_start

            if start is None:
                log.warning(f"Skipping allele '{allele_id}': missing position.")
                continue

            if chr_key is None:
                if ref_default_chr is None:
                    log.warning(
                        f"Skipping allele '{allele_id}': missing chromosome and reference has multiple sequences."
                    )
                    continue
                chr_key = ref_default_chr

            chr_key = str(chr_key)
            if chr_key not in ref_seqs:
                log.warning(
                    f"Skipping allele '{allele_id}': chromosome '{chr_key}' not found in reference."
                )
                continue

            alt_base = (alt_base or "").upper()
            ref_base = (ref_base or "").upper()
            if not alt_base:
                log.warning(f"Skipping allele '{allele_id}': missing alt base.")
                continue

            if snp_only:
                if len(alt_base) != 1 or (ref_base and len(ref_base) != 1):
                    continue
                if dec_end is not None and dec_start is not None and dec_end != dec_start:
                    continue
                span_len = 1
            else:
                if dec_end is not None:
                    length_start = dec_start if dec_start is not None else start
                    if length_start is not None and dec_end >= length_start:
                        span_len = dec_end - length_start + 1
                    else:
                        span_len = 0
                elif ref_base:
                    span_len = len(ref_base)
                else:
                    span_len = len(alt_base)
                if span_len <= 0 or len(alt_base) != span_len:
                    log.warning(
                        f"Skipping allele '{allele_id}': unsupported length change for non-SNP allele."
                    )
                    continue

            ref_seq = ref_seqs[chr_key]
            start_int = int(start)
            pos_candidates = []
            for pos in (start_int - 1, start_int):
                if pos not in pos_candidates:
                    pos_candidates.append(pos)
            pos_candidates = [
                pos for pos in pos_candidates
                if 0 <= pos <= len(ref_seq) - span_len
            ]
            if not pos_candidates:
                log.warning(
                    f"Skipping allele '{allele_id}': position out of reference bounds."
                )
                continue

            if ref_base:
                match_pos = None
                for pos in pos_candidates:
                    if ref_seq[pos:pos + span_len] == ref_base:
                        match_pos = pos
                        break
                if match_pos is None:
                    log.warning(
                        f"Skipping allele '{allele_id}': reference base mismatch at {chr_key}:{start}."
                    )
                    continue
                pos0 = match_pos
            else:
                pos0 = pos_candidates[0]

            allele_info[idx] = (chr_key, pos0, span_len, alt_base)

        ## convert reference sequences to mutable lists once
        ref_seq_lists = {name: list(seq) for name, seq in ref_seqs.items()}

        ## guard against overwriting existing outputs
        out_paths = {guid: os.path.join(out_directory, f"{guid}.fasta") for guid in guids}
        existing = [path for path in out_paths.values() if os.path.exists(path)]
        if existing:
            raise MatrixWriteError(f"Output FASTA already exists: {existing[0]}")

        for guid in guids:
            guid_idx = self._headerUtils.encode_guid(guid)
            allele_indices = self._matrix.getSetBitIndices(guid_idx, backend="auto")

            ## apply allele edits to a per-guid copy
            guid_seq_lists = {name: seq[:] for name, seq in ref_seq_lists.items()}
            applied_sites = {}

            for allele_idx in allele_indices:
                info = allele_info[int(allele_idx)]
                if info is None:
                    continue
                chr_key, pos0, span_len, alt_base = info
                site_key = (chr_key, pos0, span_len)
                existing = applied_sites.get(site_key)
                if existing is not None and existing != alt_base:
                    log.warning(
                        f"GUID '{guid}' has multiple alleles at {chr_key}:{pos0 + 1}. "
                        f"Keeping '{existing}', skipping '{alt_base}'."
                    )
                    continue
                applied_sites[site_key] = alt_base

                seq_list = guid_seq_lists[chr_key]
                if span_len == 1:
                    seq_list[pos0] = alt_base
                else:
                    seq_list[pos0:pos0 + span_len] = list(alt_base)

            ## write modified sequences to FASTA
            records = []
            for name, seq_list in guid_seq_lists.items():
                base_record = ref_records[name]
                description = base_record.description or name
                meta_parts = []
                if f"guid={guid}" not in description:
                    meta_parts.append(f"guid={guid}")
                if "generated_by=" not in description:
                    meta_parts.append(f"generated_by={generated_by}")
                if ardal_version and "ardal_version=" not in description:
                    meta_parts.append(f"ardal_version={ardal_version}")
                if meta_parts:
                    description = f"| ref={ref} | length={len(seq_list)} | {' '.join(meta_parts)} snp_only={snp_only}"
                record = SeqRecord.SeqRecord(
                    Seq.Seq("".join(seq_list)),
                    id=base_record.id,
                    name=base_record.name,
                    description=description
                )
                records.append(record)
            out_path = out_paths[guid]
            SeqIO.write(records, out_path, "fasta")

        return out_paths


    def make_alignment( self,
                        output_prefix : str,
                        guids : list = [],
                        ref : Union[str, None] = None,
                        allele_id_format : Union[str, None] = None,
                        out_directory : str = "./",
                        snp_only : bool = True,
                        missing_char : str = "N"
                        ) -> str:
        """Write a multi-FASTA alignment of ordered alleles.

        When a reference FASTA is provided, absent alleles are represented by the
        reference base at each locus. When no reference is provided, the reference
        base encoded in the allele IDs (if present) is used; otherwise `missing_char`
        is emitted for absent alleles.
        """
        from Bio import SeqIO, Seq, SeqRecord

        if not output_prefix:
            raise ParameterError("Output prefix is required.")
        if not out_directory:
            raise ParameterError("Output directory is required.")
        if not os.path.exists(out_directory):
            raise ParameterError(f"Output directory not found: {out_directory}")
        if not os.path.isdir(out_directory):
            raise ParameterError(f"Output directory is not a directory: {out_directory}")

        if not isinstance(missing_char, str) or len(missing_char) != 1:
            raise ParameterError("missing_char must be a single character.")

        if guids:
            self._headerUtils.check_guids(guids)
        else:
            guids = list(self._headerUtils.headers.get("guids", []))

        if not guids:
            raise ParameterError("No GUIDs available for alignment.")

        if allele_id_format is None:
            allele_id_format = self._headerUtils.get_cached_allele_id_format()
        if not allele_id_format:
            raise ParameterError("allele_id_format is required to build an alignment.")

        pattern = self._headerUtils._check_allele_format_grammar(allele_id_format=allele_id_format)
        self._headerUtils.ensure_id_positions(allele_id_format)
        allele_position_map = self._headerUtils.get_allele_positions()
        alleles = self._headerUtils.headers.get("alleles", [])

        ref_records = None
        ref_seqs = None
        ref_default_chr = None
        if ref is not None:
            if not os.path.exists(ref):
                raise ParameterError(f"Reference FASTA not found: {ref}")
            try:
                ref_records = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
            except Exception as exc:
                raise MatrixWriteError(f"Failed to parse reference FASTA: {exc}") from exc
            if not ref_records:
                raise MatrixWriteError("Reference FASTA appears to be empty.")
            ref_seqs = {
                name: str(record.seq).upper()
                for name, record in ref_records.items()
            }
            if len(ref_seqs) == 1:
                ref_default_chr = next(iter(ref_seqs.keys()))

        ## build ordered allele list with metadata
        allele_records = []
        for idx, allele_id in enumerate(alleles):
            coord = allele_position_map.get(allele_id)
            chr_key = coord[0] if coord else None
            start = coord[1] if coord else None

            try:
                dec_chr, dec_start, dec_end, ref_base, alt_base = self._headerUtils._decode_allele_position(
                    allele_id=allele_id,
                    pattern=pattern,
                    allele_id_format=allele_id_format
                )
            except ValueError:
                dec_chr = None
                dec_start = None
                dec_end = None
                ref_base = None
                alt_base = None

            if chr_key is None:
                chr_key = dec_chr
            if start is None:
                start = dec_start

            if start is None:
                log.warning(f"Skipping allele '{allele_id}': missing position.")
                continue

            if chr_key is None:
                if ref_default_chr is None:
                    log.warning(
                        f"Skipping allele '{allele_id}': missing chromosome and reference has multiple sequences."
                    )
                    continue
                chr_key = ref_default_chr

            chr_key = str(chr_key)

            alt_base = (alt_base or "").upper()
            ref_base = (ref_base or "").upper()
            if not alt_base:
                log.warning(f"Skipping allele '{allele_id}': missing alt base.")
                continue

            if len(alt_base) != 1:
                if snp_only:
                    continue
                log.warning(f"Skipping allele '{allele_id}': non-SNP allele not supported in alignment.")
                continue
            if ref_base and len(ref_base) != 1:
                if snp_only:
                    continue
                log.warning(f"Skipping allele '{allele_id}': non-SNP allele not supported in alignment.")
                continue
            if snp_only and dec_end is not None and dec_start is not None and dec_end != dec_start:
                continue

            ref_char = None
            if ref_seqs is not None:
                if chr_key not in ref_seqs:
                    log.warning(
                        f"Skipping allele '{allele_id}': chromosome '{chr_key}' not found in reference."
                    )
                    continue
                ref_seq = ref_seqs[chr_key]
                start_int = int(start)
                pos_candidates = []
                for pos in (start_int - 1, start_int):
                    if pos not in pos_candidates:
                        pos_candidates.append(pos)
                pos_candidates = [pos for pos in pos_candidates if 0 <= pos < len(ref_seq)]
                if not pos_candidates:
                    log.warning(
                        f"Skipping allele '{allele_id}': position out of reference bounds."
                    )
                    continue

                if ref_base:
                    match_pos = None
                    for pos in pos_candidates:
                        if ref_seq[pos] == ref_base:
                            match_pos = pos
                            break
                    if match_pos is None:
                        log.warning(
                            f"Skipping allele '{allele_id}': reference base mismatch at {chr_key}:{start}."
                        )
                        continue
                    pos0 = match_pos
                else:
                    pos0 = pos_candidates[0]
                ref_char = ref_seq[pos0]
            else:
                if ref_base and len(ref_base) == 1:
                    ref_char = ref_base

            allele_records.append({
                "idx": idx,
                "chr": chr_key,
                "pos": int(start),
                "id": allele_id,
                "alt": alt_base,
                "ref": ref_char,
            })

        if not allele_records:
            raise ParameterError("No alleles available for alignment.")

        allele_records.sort(key=lambda rec: (rec["chr"], rec["pos"], rec["id"]))
        ordered_indices = [rec["idx"] for rec in allele_records]
        ordered_alt = [rec["alt"] for rec in allele_records]
        ordered_ref = [rec["ref"] for rec in allele_records]

        out_path = os.path.join(out_directory, f"{output_prefix}.fasta")
        if os.path.exists(out_path):
            raise MatrixWriteError(f"File '{out_path}' already exists.")

        missing_rows = self._headerUtils.get_missing_mask_rows()
        missing_sets = [set(row) for row in missing_rows] if missing_rows else None

        try:
            ardal_version = version("ardal")
        except PackageNotFoundError:
            ardal_version = None

        records = []
        for guid in guids:
            guid_idx = self._headerUtils.encode_guid(guid)
            allele_indices = self._matrix.getSetBitIndices(guid_idx, backend="auto")
            present = set(int(idx) for idx in allele_indices)
            missing_set = missing_sets[guid_idx] if missing_sets else None

            seq_chars = []
            for allele_idx, alt_base, ref_char in zip(ordered_indices, ordered_alt, ordered_ref):
                if missing_set is not None and allele_idx in missing_set:
                    seq_chars.append(missing_char)
                elif allele_idx in present:
                    seq_chars.append(alt_base)
                else:
                    seq_chars.append(ref_char if ref_char else missing_char)

            desc_parts = [f"sites={len(seq_chars)}", f"snp_only={snp_only}"]
            if ref:
                desc_parts.append(f"ref={ref}")
            if ardal_version:
                desc_parts.append(f"ardal_version={ardal_version}")
            description = " ".join(desc_parts)

            record = SeqRecord.SeqRecord(
                Seq.Seq("".join(seq_chars)),
                id=guid,
                name=guid,
                description=description
            )
            records.append(record)

        SeqIO.write(records, out_path, "fasta")
        return out_path


    def to_plink( self,
                  out_prefix: str,
                  out_directory: str = "./",
                  allele_id_format: Union[str, None] = None,
                  snp_only: bool = True,
                  chr_to_int: bool = False,
                  chunk_size: Union[int, None] = 10000
                  ) -> None:
        """Write PLINK1 binary files (.bed/.bim/.fam) for the current matrix.

        Notes:
            - Allele presence is encoded as heterozygous (01).
            - Absence is encoded as homozygous reference (10).
            - Missing masks are encoded as missing (11).
            - chr_to_int optionally recodes chromosome labels to PLINK numeric codes.
            - When chr_to_int is True, a sidecar <prefix>.chrmap file is written.
        """
        if not out_prefix:
            raise ParameterError("Output prefix is required.")
        if not out_directory:
            raise ParameterError("Output directory is required.")
        if not os.path.exists(out_directory):
            raise ParameterError(f"Output directory not found: {out_directory}")
        if not os.path.isdir(out_directory):
            raise ParameterError(f"Output directory is not a directory: {out_directory}")

        bed_path = os.path.join(out_directory, f"{out_prefix}.bed")
        bim_path = os.path.join(out_directory, f"{out_prefix}.bim")
        fam_path = os.path.join(out_directory, f"{out_prefix}.fam")

        for path in (bed_path, bim_path, fam_path):
            if os.path.exists(path):
                raise MatrixWriteError(f"File '{path}' already exists.")

        if allele_id_format is None:
            allele_id_format = self._headerUtils.get_cached_allele_id_format()
        if not allele_id_format:
            raise ParameterError("allele_id_format is required to generate PLINK outputs.")

        ## decode variant metadata from allele IDs
        pattern = self._headerUtils._check_allele_format_grammar(allele_id_format=allele_id_format)
        self._headerUtils.ensure_id_positions(allele_id_format)
        allele_position_map = self._headerUtils.get_allele_positions()

        ## optional chromosome recoding for PLINK
        def _standard_chr_code(raw_value: str) -> Union[int, None]:
            raw = str(raw_value)
            norm = raw[3:] if raw.lower().startswith("chr") else raw
            norm_upper = norm.upper()
            chr_map = {
                "X": 23,
                "Y": 24,
                "XY": 25,
                "M": 26,
                "MT": 26,
            }
            if norm_upper in chr_map:
                return chr_map[norm_upper]
            try:
                return int(norm)
            except ValueError:
                return None

        contig_map = {}
        used_codes = set()
        max_code = 0

        alleles = self._headerUtils.headers.get("alleles", [])
        variants = []
        for idx, allele_id in enumerate(alleles):
            try:
                chr_key, start, end, ref_base, alt_base = self._headerUtils._decode_allele_position(
                    allele_id=allele_id,
                    pattern=pattern,
                    allele_id_format=allele_id_format
                )
            except ValueError:
                coord = allele_position_map.get(allele_id)
                if coord:
                    chr_key, start = coord
                else:
                    chr_key = None
                    start = None
                end = None
                ref_base = None
                alt_base = None

            if chr_key is None or start is None:
                log.warning(f"Skipping allele '{allele_id}': missing position.")
                continue

            ref_base = (ref_base or "").upper()
            alt_base = (alt_base or "").upper()

            if snp_only:
                if not alt_base or len(alt_base) != 1:
                    continue
                if ref_base and len(ref_base) != 1:
                    continue
                if end is not None and end != start:
                    continue

            raw_chr = str(chr_key)
            if chr_to_int:
                std_code = _standard_chr_code(raw_chr)
                if std_code is not None:
                    used_codes.add(std_code)
                    if std_code > max_code:
                        max_code = std_code

            variants.append({
                "idx": idx,
                "id": allele_id,
                "chr_raw": raw_chr,
                "pos": int(start),
                "ref": ref_base or "0",
                "alt": alt_base or "0"
            })

        if not variants:
            raise ParameterError("No variants available for PLINK export.")

        if chr_to_int:
            ## assign numeric codes for non-standard contigs
            next_code = max(26, max_code) + 1
            for variant in variants:
                raw = variant["chr_raw"]
                std_code = _standard_chr_code(raw)
                if std_code is not None:
                    variant["chr"] = str(std_code)
                    continue
                if raw in contig_map:
                    variant["chr"] = str(contig_map[raw])
                    continue
                while next_code in used_codes:
                    next_code += 1
                contig_map[raw] = next_code
                used_codes.add(next_code)
                variant["chr"] = str(next_code)
                next_code += 1
        else:
            for variant in variants:
                variant["chr"] = variant["chr_raw"]

        guids = list(self._headerUtils.headers.get("guids", []))
        if not guids:
            raise ParameterError("No GUIDs available for PLINK export.")

        ## write BIM and FAM using stable header order
        with open(bim_path, "w") as bim:
            for variant in variants:
                bim.write(
                    f"{variant['chr']}\t{variant['id']}\t0\t{variant['pos']}\t"
                    f"{variant['alt']}\t{variant['ref']}\n"
                )

        with open(fam_path, "w") as fam:
            for guid in guids:
                fam.write(f"{guid}\t{guid}\t0\t0\t0\t-9\n")

        if chr_to_int and contig_map:
            ## write contig mapping for non-standard chromosome labels
            map_path = os.path.join(out_directory, f"{out_prefix}.chrmap")
            if os.path.exists(map_path):
                raise MatrixWriteError(f"File '{map_path}' already exists.")
            with open(map_path, "w") as cmap:
                for raw, code in contig_map.items():
                    cmap.write(f"{raw}\t{code}\n")

        if chunk_size is None or chunk_size <= 0:
            chunk_size = len(variants)

        missing_rows = self._headerUtils.get_missing_mask_rows()
        missing_sets = [set(row) for row in missing_rows] if missing_rows else None

        guid_indices = list(range(len(guids)))
        n_samples = len(guids)

        hom_ref_code = np.uint8(2)
        het_code = np.uint8(1)
        missing_code = np.uint8(3)

        ## write SNP major BED stream
        with open(bed_path, "wb") as bed:
            bed.write(bytes([0x6C, 0x1B, 0x01]))

            for chunk_start in range(0, len(variants), chunk_size):
                chunk_vars = variants[chunk_start:chunk_start + chunk_size]
                allele_indices = [v["idx"] for v in chunk_vars]
                if not allele_indices:
                    continue

                chunk_packed = self._matrix.getSubsetPackedMatrix(guid_indices, allele_indices, 1)
                if chunk_packed.dtype.byteorder not in ('<', '='):
                    chunk_packed = chunk_packed.byteswap().newbyteorder('<')
                chunk_packed = np.asarray(chunk_packed, dtype=np.uint64, copy=False)
                if chunk_packed.ndim == 1:
                    chunk_packed = chunk_packed.reshape(1, chunk_packed.shape[0])

                for local_idx, variant in enumerate(chunk_vars):
                    ## extract presence per sample for this variant
                    word = local_idx // 64
                    bit = local_idx % 64
                    presence = (chunk_packed[:, word] >> np.uint64(bit)) & np.uint64(1)
                    genotypes = np.full(n_samples, hom_ref_code, dtype=np.uint8)
                    if presence.any():
                        genotypes[presence.astype(bool)] = het_code

                    if missing_sets:
                        ## apply missing mask (guid x allele index)
                        missing_mask = np.fromiter(
                            (variant["idx"] in missing_set for missing_set in missing_sets),
                            dtype=bool,
                            count=n_samples
                        )
                        if missing_mask.any():
                            genotypes[missing_mask] = missing_code

                    ## pack 4 samples per byte (PLINK bed spec)
                    n_bytes = (n_samples + 3) // 4
                    out_bytes = np.zeros(n_bytes, dtype=np.uint8)
                    for byte_idx in range(n_bytes):
                        base = byte_idx * 4
                        g0 = int(genotypes[base]) if base < n_samples else 0
                        g1 = int(genotypes[base + 1]) if base + 1 < n_samples else 0
                        g2 = int(genotypes[base + 2]) if base + 2 < n_samples else 0
                        g3 = int(genotypes[base + 3]) if base + 3 < n_samples else 0
                        out_bytes[byte_idx] = g0 | (g1 << 2) | (g2 << 4) | (g3 << 6)

                    bed.write(out_bytes.tobytes())
