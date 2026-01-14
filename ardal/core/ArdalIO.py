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
            headers_meta = {"meta": meta, "headers": self._headerUtils.headers}
        elif format == "npz":
            log.info("Writing dense matrix as .npz")
            np.savez_compressed(matrix_out_path, matrix=matrix_to_save)
            meta = make_meta(matrix_to_save,
                             self._headerUtils.headers,
                             generated_by="ardal::io::write",
                             format_name="ardal.dense.v1",
                             matrix_file=matrix_out_path)
            headers_meta = {"meta": meta, "headers": self._headerUtils.headers}
        elif format == "bin":
            log.info("Writing packed matrix as .bin")
            headers_meta = self._write_packed(matrix_out_path)
        
        ## update with missing sites data
        headers_meta.update(self._headerUtils._missing_masks)
            
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
        headers_meta = { "meta" : meta, "headers" : self._headerUtils.headers }
        
        return headers_meta
        
    
    @staticmethod
    def _generated_by_string() -> str:
        return "ardal::io::write"
    

    def make_fastas( self,
                     guids : list = [],
                     ref : Union[str, None] = None,
                     allele_id_format : str = "{ref}.{chr}.{start}.{alt}"
                     ) -> None:
        """ Takes a set of guids and a reference fasta and constructs a simulated fasta from each dataset using alleles stored
        in the allele matrix. The position of each allele is determined by the allele_id_format argument, whereby the keywords:
        ref, alt, chr, start, end, are used to outline the allele id naming convention.

        E.g. for the allele `A.chr1.100.101.T` would be decoded using the allele_id_format string `{ref}.{chr}.{start}.{end}.{alt}`.
        """
        SeqIO = require_package("Bio", attr="SeqIO")
        Seq = require_package("Bio.Seq", import_as="Bio.Seq", attr="Seq")
        SeqRecord = require_package("Bio.SeqRecord", import_as="Bio.SeqRecord", attr="SeqRecord")

        if ref is None:
            raise ParameterError("Reference FASTA path is required.")
        if not os.path.exists(ref):
            raise ParameterError(f"Reference FASTA not found: {ref}")

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

        ## decode alleles once
        alleles = self._headerUtils.headers.get("alleles", [])
        allele_info: list[Union[tuple[str, int, Union[int, None], str, str], None]] = [None] * len(alleles)
        for idx, allele_id in enumerate(alleles):
            try:
                chr_key, start, end, ref_base, alt_base = self._headerUtils._decode_allele_position(
                    allele_id=allele_id,
                    pattern=pattern,
                    allele_id_format=allele_id_format
                )
            except ValueError:
                continue
            if chr_key is None or start is None:
                continue
            chr_key = str(chr_key)
            if chr_key not in ref_seqs:
                continue
            ref_base = "" if ref_base is None else str(ref_base)
            alt_base = "" if alt_base is None else str(alt_base)
            allele_info[idx] = (chr_key, int(start), end, ref_base, alt_base)

        ## infer per contig coordinate offsets
        ## switch from 0based to 1based
        contig_offset: dict[str, int] = {}
        for contig, seq in ref_seqs.items():
            match_zero = 0
            match_one = 0
            checked = 0
            for info in allele_info:
                if not info or info[0] != contig:
                    continue
                _, start, end, ref_base, _ = info
                if not ref_base or start is None:
                    continue
                if end is None:
                    end = start + len(ref_base)
                try:
                    end_val = int(end)
                except (TypeError, ValueError):
                    continue
                if start < 0:
                    continue
                if start + len(ref_base) > len(seq) and start - 1 + len(ref_base) > len(seq):
                    continue
                checked += 1
                if start + len(ref_base) <= len(seq):
                    if seq[start:start + len(ref_base)] == ref_base.upper():
                        match_zero += 1
                if start - 1 >= 0 and start - 1 + len(ref_base) <= len(seq):
                    if seq[start - 1:start - 1 + len(ref_base)] == ref_base.upper():
                        match_one += 1
            if checked == 0:
                contig_offset[contig] = 0
                log.warning(f"No reference matches found for contig '{contig}'. Assuming 0-based coordinates.")
            else:
                contig_offset[contig] = -1 if match_one >= match_zero else 0

        ## generate fastas
        for guid in guids:
            guid_idx = self._headerUtils.encode_guid(guid)
            allele_indices = self._matrix.getSetBitIndices(guid_idx, backend="auto")

            contig_variants: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
            for allele_idx in allele_indices:
                info = allele_info[int(allele_idx)]
                if info is None:
                    continue
                contig, start, end, ref_base, alt_base = info
                if contig not in ref_seqs:
                    continue
                offset = contig_offset.get(contig, 0)
                start_idx = start + offset
                if start_idx < 0:
                    continue
                if end is None:
                    if ref_base:
                        end_val = start + len(ref_base)
                    elif alt_base:
                        end_val = start + max(1, len(alt_base))
                    else:
                        end_val = start + 1
                else:
                    end_val = int(end)
                end_idx = end_val + offset
                if end_idx < start_idx:
                    continue
                contig_variants[contig].append((start_idx, end_idx, ref_base, alt_base))

            records = []
            for contig, seq in ref_seqs.items():
                variants = contig_variants.get(contig, [])
                if not variants:
                    out_seq = seq
                else:
                    variants.sort(key=lambda v: v[0])
                    cursor = 0
                    out_chunks: list[str] = []
                    for start_idx, end_idx, ref_base, alt_base in variants:
                        if start_idx < cursor:
                            log.warning(f"Overlapping variants for {guid} on {contig}; skipping at {start_idx}.")
                            continue
                        if start_idx > len(seq):
                            continue
                        end_idx = min(end_idx, len(seq))
                        if ref_base:
                            ref_slice = seq[start_idx:end_idx]
                            if ref_slice != ref_base.upper():
                                log.warning(f"Reference mismatch for {guid} on {contig} at {start_idx}; skipping.")
                                continue
                        out_chunks.append(seq[cursor:start_idx])
                        if alt_base:
                            out_chunks.append(alt_base.upper())
                        cursor = end_idx
                    out_chunks.append(seq[cursor:])
                    out_seq = "".join(out_chunks)

                record_id = f"{guid}|{contig}"
                records.append(SeqRecord(Seq(out_seq), id=record_id, description=""))

            out_path = f"{guid}.fasta"
            try:
                SeqIO.write(records, out_path, "fasta")
            except Exception as exc:
                raise MatrixWriteError(f"Failed to write FASTA for {guid}: {exc}") from exc

            log.info(f"Wrote FASTA for GUID '{guid}' to {out_path}")

        return None
