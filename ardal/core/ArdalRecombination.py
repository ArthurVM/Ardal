""" ArdalRecombination.py
This module provides functionality for investigation signals of recombination in the Ardal framework.
"""
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from collections import defaultdict
from typing import Union, Tuple, List, Dict

from ..utils.decorators import check_thread_count, check_guids_list
from ..utils.exceptions import ParameterError
from ..utils.logger import get_logger

log = get_logger()

_RECOMB_SELF = None
_RECOMB_MODEL = None


def _init_recomb_worker(model: Dict) -> None:
    global _RECOMB_MODEL
    _RECOMB_MODEL = model


def _recomb_worker(guid: str,
                   min_tract_sites: int,
                   return_loci: bool,
                   threads: int) -> Tuple[str, Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]]]:
    recomb = _RECOMB_SELF
    model = _RECOMB_MODEL
    if recomb is None or model is None:
        raise RuntimeError("Recombination worker is not initialised.")

    obs, missing_mask = recomb._extract_guid_observations(
        guid=guid,
        allele_indices=model["allele_idx_arr"],
        threads=threads,
    )

    out: Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]] = {"windows": []}
    if return_loci:
        out["loci"] = []

    for chr_key, (start, end) in model["chr_slices"].items():
        if start == end:
            continue
        obs_slice = obs[start:end]
        missing_slice = missing_mask[start:end]
        emit = np.where(obs_slice, model["logp"][:, start:end], model["log1p"][:, start:end])
        if missing_slice.any():
            emit[:, missing_slice] = 0.0

        states = recomb._viterbi(emit, model["log_trans"], model["log_start"])

        if return_loci:
            for offset, state in enumerate(states):
                out["loci"].append({
                    "chr": chr_key,
                    "pos": model["ordered_pos"][start + offset],
                    "allele": model["ordered_alleles"][start + offset],
                    "population": model["pop_names"][int(state)],
                })

        if states.size == 0:
            continue

        segments: List[List[int]] = []
        seg_start = 0
        current = int(states[0])
        for t in range(1, states.size):
            if states[t] != current:
                segments.append([seg_start, t - 1, current])
                seg_start = t
                current = int(states[t])
        segments.append([seg_start, states.size - 1, current])

        segments = recomb._merge_short_segments(segments, min_tract_sites)
        segments = recomb._collapse_adjacent_segments(segments)

        for seg_start, seg_end, state in segments:
            abs_start = start + seg_start
            abs_end = start + seg_end
            window_scores = emit[:, seg_start:seg_end + 1].sum(axis=1)
            window_log_scores = window_scores + model["log_start"]
            max_score = float(np.max(window_log_scores))
            exp_scores = np.exp(window_log_scores - max_score)
            exp_sum = float(exp_scores.sum())
            if exp_sum == 0.0:
                posterior = [1.0 / len(model["pop_names"])] * len(model["pop_names"])
            else:
                posterior = (exp_scores / exp_sum).tolist()
            loglik = {
                model["pop_names"][idx]: float(score)
                for idx, score in enumerate(window_scores)
            }
            post = {
                model["pop_names"][idx]: float(score)
                for idx, score in enumerate(posterior)
            }
            out["windows"].append({
                "chr": chr_key,
                "start": model["ordered_pos"][abs_start],
                "end": model["ordered_pos"][abs_end],
                "n_sites": abs_end - abs_start + 1,
                "population": model["pop_names"][int(state)],
                "loglik": loglik,
                "posterior": post,
            })

    return guid, out

## core/ArdalRecombination.py
class ArdalRecombination:
    """ Class for investigation signals of recombination in the Ardal framework.
    """
    
    def __init__( self,
                  headerUtils,
                  hybrid_matrix,
                  roaring_enabled : bool ):

        self._headerUtils = headerUtils
        self._matrix = hybrid_matrix
        self.roaring = roaring_enabled
    
    
    def _build_ordered_snv_loci(self,
                                allele_id_format: str) -> Tuple[List[int], List[int], List[str], Dict[str, Tuple[int, int]]]:
        allele_to_locus, _, _ = self._headerUtils.get_snv_lookup(allele_id_format=allele_id_format)
        invalid = np.uint32(np.iinfo(np.uint32).max)
        valid_indices = np.flatnonzero(allele_to_locus != invalid)

        if valid_indices.size == 0:
            raise ParameterError("No biallelic SNV alleles were resolved for tract detection.")

        self._headerUtils.ensure_id_positions(allele_id_format)
        allele_positions = self._headerUtils.get_allele_positions()
        alleles = self._headerUtils.headers["alleles"]

        chr_groups: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
        missing_pos = 0
        for idx in valid_indices:
            allele_id = alleles[int(idx)]
            coord = allele_positions.get(allele_id)
            if coord is None:
                missing_pos += 1
                continue
            chr_key, pos = coord
            chr_groups[str(chr_key)].append((int(pos), int(idx), allele_id))

        if missing_pos:
            log.warning(f"Skipping {missing_pos} alleles with missing positions.")

        ordered_indices: List[int] = []
        ordered_pos: List[int] = []
        ordered_alleles: List[str] = []
        chr_slices: Dict[str, Tuple[int, int]] = {}

        cursor = 0
        for chr_key in sorted(chr_groups.keys()):
            items = sorted(chr_groups[chr_key], key=lambda row: row[0])
            start = cursor
            for pos, idx, allele_id in items:
                ordered_indices.append(idx)
                ordered_pos.append(pos)
                ordered_alleles.append(allele_id)
                cursor += 1
            chr_slices[chr_key] = (start, cursor)

        if not ordered_indices:
            raise ParameterError("No SNV loci with positions were available for tract detection.")

        return ordered_indices, ordered_pos, ordered_alleles, chr_slices


    def _extract_guid_observations(self,
                                   guid: str,
                                   allele_indices: np.ndarray,
                                   threads: int) -> Tuple[np.ndarray, np.ndarray]:
        guid_idx = self._headerUtils.encode_guid(guid)
        packed = self._matrix.getSubsetPackedMatrix_rows([guid_idx], threads)
        if packed.dtype.byteorder not in ('<', '='):
            packed = packed.byteswap().newbyteorder('<')
        packed = np.asarray(packed, dtype=np.uint64, copy=False)
        if packed.ndim == 1:
            packed = packed.reshape(1, packed.shape[0])

        word_idx = (allele_indices >> np.uint64(6)).astype(np.intp, copy=False)
        bit_offsets = (allele_indices & np.uint64(63)).astype(np.uint64, copy=False)
        selected = np.take(packed, word_idx, axis=1)
        shifted = np.right_shift(selected, bit_offsets[np.newaxis, :])
        obs = (shifted & np.uint64(1)).astype(np.uint8, copy=False)[0]

        if self._headerUtils.has_missing_mask():
            missing_cols = self._headerUtils.get_guid_missing_mask(guid)
            if missing_cols:
                missing_mask = np.isin(allele_indices, np.asarray(missing_cols, dtype=np.uint64))
            else:
                missing_mask = np.zeros(obs.shape[0], dtype=bool)
        else:
            missing_mask = np.zeros(obs.shape[0], dtype=bool)

        return obs, missing_mask


    @staticmethod
    def _viterbi(emit: np.ndarray,
                 log_trans: np.ndarray,
                 log_start: np.ndarray) -> np.ndarray:
        n_states, n_sites = emit.shape
        if n_sites == 0:
            return np.zeros(0, dtype=np.int16)

        dp = np.empty((n_states, n_sites), dtype=np.float64)
        ptr = np.empty((n_states, n_sites), dtype=np.int16)

        dp[:, 0] = log_start + emit[:, 0]
        ptr[:, 0] = -1

        for t in range(1, n_sites):
            scores = dp[:, t - 1][:, None] + log_trans
            best_prev = np.argmax(scores, axis=0)
            dp[:, t] = scores[best_prev, np.arange(n_states)] + emit[:, t]
            ptr[:, t] = best_prev

        states = np.empty(n_sites, dtype=np.int16)
        states[-1] = int(np.argmax(dp[:, -1]))
        for t in range(n_sites - 1, 0, -1):
            states[t - 1] = ptr[states[t], t]

        return states


    @staticmethod
    def _merge_short_segments(segments: List[List[int]],
                              min_sites: int) -> List[List[int]]:
        if min_sites <= 1:
            return segments

        i = 0
        while i < len(segments):
            start_idx, end_idx, _ = segments[i]
            length = end_idx - start_idx + 1
            if length >= min_sites or len(segments) == 1:
                i += 1
                continue

            if i == 0:
                merge_idx = 1
            elif i == len(segments) - 1:
                merge_idx = i - 1
            else:
                prev_len = segments[i - 1][1] - segments[i - 1][0] + 1
                next_len = segments[i + 1][1] - segments[i + 1][0] + 1
                merge_idx = i - 1 if prev_len >= next_len else i + 1

            if merge_idx < i:
                segments[merge_idx][1] = end_idx
            else:
                segments[merge_idx][0] = start_idx

            del segments[i]
            i = max(merge_idx - 1, 0)

        return segments


    @staticmethod
    def _collapse_adjacent_segments(segments: List[List[int]]) -> List[List[int]]:
        if not segments:
            return segments
        collapsed = [segments[0]]
        for start_idx, end_idx, state in segments[1:]:
            last = collapsed[-1]
            if last[2] == state and last[1] + 1 >= start_idx:
                last[1] = max(last[1], end_idx)
            else:
                collapsed.append([start_idx, end_idx, state])
        return collapsed


    @check_thread_count
    @check_guids_list
    def detect_tracts( self, 
                       anchor_pops : Dict, 
                       allele_id_format : Union[str, None] = None, 
                       guids : Union[List, None] = None, 
                       expected_tract_len : int = 200, 
                       min_tract_sites : int = 5, 
                       pseudocount : float = 0.5, 
                       return_loci : bool = False,
                       threads : int = 1 ) -> Dict[str, Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]]]:
        """Detect recombination tracts by assigning ancestry along ordered SNV loci.

        Args:
            guids: Query GUIDs to segment.
            anchor_pops: Mapping of population name -> list of GUIDs.
            allele_id_format: Format string used to decode SNV loci.
            expected_tract_len: Expected tract length in sites (controls switch penalty).
            min_tract_sites: Minimum tract length in sites for reporting.
            pseudocount: Pseudocount for allele frequency smoothing.
            return_loci: When True, return per-locus assignments alongside tracts.
            threads: Thread count for matrix extraction.

        Returns:
            Dict keyed by GUID, each containing "windows" entries (and optional "loci")
            with per-window population assignments, per-population log-likelihoods,
            and relative posteriors computed via softmax.
        """
        if not guids:
            guids = self._headerUtils.headers["guids"]
        self._headerUtils.check_guids(guids)
        if not anchor_pops or not isinstance(anchor_pops, dict):
            raise ParameterError("anchor_pops must be a non-empty dict of {population: [guids]}.")
        if expected_tract_len < 1:
            raise ParameterError("expected_tract_len must be at least 1.")
        if min_tract_sites < 1:
            raise ParameterError("min_tract_sites must be at least 1.")
        if pseudocount < 0:
            raise ParameterError("pseudocount must be non-negative.")
        
        ## check if we need to inherit allele_id_format from headerUtils
        if allele_id_format is None:
            allele_id_format = self._headerUtils._allele_id_format
        ## check it was inherited, else raise an error
        if allele_id_format is None:
            raise ParameterError("allele_id_format is None and cannot inherit from _headerUtils. Please provide an allele_id_format.")

        ordered_indices, ordered_pos, ordered_alleles, chr_slices = self._build_ordered_snv_loci(
            allele_id_format=allele_id_format
        )

        allele_idx_arr = np.asarray(ordered_indices, dtype=np.uint64)
        n_loci = allele_idx_arr.size
        pop_names = list(anchor_pops.keys())
        if not pop_names:
            raise ParameterError("anchor_pops must define at least one population.")

        log.info(f"Preparing recombination model with {len(pop_names)} populations and {n_loci} SNV loci.")

        eps = 1e-6
        logp = np.empty((len(pop_names), n_loci), dtype=np.float64)
        log1p = np.empty((len(pop_names), n_loci), dtype=np.float64)

        for i, pop in enumerate(pop_names):
            pop_guids = anchor_pops.get(pop, [])
            if not pop_guids:
                raise ParameterError(f"anchor_pops[{pop}] cannot be empty.")
            self._headerUtils.check_guids(pop_guids)
            pop_coords = [self._headerUtils.encode_guid(guid) for guid in pop_guids]
            freqs = np.asarray(self._matrix.colFrequency(pop_coords), dtype=np.float64)
            freqs = freqs[allele_idx_arr]

            n = float(len(pop_guids))
            p = (freqs * n + pseudocount) / (n + 2.0 * pseudocount)
            p = np.clip(p, eps, 1.0 - eps)
            logp[i, :] = np.log(p)
            log1p[i, :] = np.log1p(-p)

        if len(pop_names) == 1:
            log_trans = np.zeros((1, 1), dtype=np.float64)
            log_start = np.zeros(1, dtype=np.float64)
        else:
            p_switch = 1.0 / float(expected_tract_len)
            p_switch = min(max(p_switch, eps), 1.0 - eps)
            stay = 1.0 - p_switch
            switch = p_switch / float(len(pop_names) - 1)
            log_trans = np.full((len(pop_names), len(pop_names)), np.log(switch), dtype=np.float64)
            np.fill_diagonal(log_trans, np.log(stay))
            log_start = np.full(len(pop_names), -np.log(len(pop_names)), dtype=np.float64)

        results: Dict[str, Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]]] = {}

        for guid in guids:
            results[guid] = {"windows": []}
            if return_loci:
                results[guid]["loci"] = []
            obs, missing_mask = self._extract_guid_observations(
                guid=guid,
                allele_indices=allele_idx_arr,
                threads=threads,
            )

            for chr_key, (start, end) in chr_slices.items():
                if start == end:
                    continue
                obs_slice = obs[start:end]
                missing_slice = missing_mask[start:end]
                emit = np.where(obs_slice, logp[:, start:end], log1p[:, start:end])
                if missing_slice.any():
                    emit[:, missing_slice] = 0.0

                states = self._viterbi(emit, log_trans, log_start)

                if return_loci:
                    for offset, state in enumerate(states):
                        results[guid]["loci"].append({
                            "chr": chr_key,
                            "pos": ordered_pos[start + offset],
                            "allele": ordered_alleles[start + offset],
                            "population": pop_names[int(state)],
                        })

                if states.size == 0:
                    continue

                segments: List[List[int]] = []
                seg_start = 0
                current = int(states[0])
                for t in range(1, states.size):
                    if states[t] != current:
                        segments.append([seg_start, t - 1, current])
                        seg_start = t
                        current = int(states[t])
                segments.append([seg_start, states.size - 1, current])

                segments = self._merge_short_segments(segments, min_tract_sites)
                segments = self._collapse_adjacent_segments(segments)

                for seg_start, seg_end, state in segments:
                    abs_start = start + seg_start
                    abs_end = start + seg_end
                    window_scores = emit[:, seg_start:seg_end + 1].sum(axis=1)
                    window_log_scores = window_scores + log_start
                    max_score = float(np.max(window_log_scores))
                    exp_scores = np.exp(window_log_scores - max_score)
                    exp_sum = float(exp_scores.sum())
                    if exp_sum == 0.0:
                        posterior = [1.0 / len(pop_names)] * len(pop_names)
                    else:
                        posterior = (exp_scores / exp_sum).tolist()
                    loglik = {
                        pop_names[idx]: float(score)
                        for idx, score in enumerate(window_scores)
                    }
                    post = {
                        pop_names[idx]: float(score)
                        for idx, score in enumerate(posterior)
                    }
                    results[guid]["windows"].append({
                        "chr": chr_key,
                        "start": ordered_pos[abs_start],
                        "end": ordered_pos[abs_end],
                        "n_sites": abs_end - abs_start + 1,
                        "population": pop_names[int(state)],
                        "loglik": loglik,
                        "posterior": post,
                    })

        return results


    @check_thread_count
    @check_guids_list
    def detect_tracts_parallel( self,
                                anchor_pops: Dict,
                                allele_id_format: Union[str, None] = None,
                                guids: Union[List, None] = None,
                                expected_tract_len: int = 200,
                                min_tract_sites: int = 5,
                                pseudocount: float = 0.5,
                                return_loci: bool = False,
                                threads: int = 1,
                                max_workers: Union[int, None] = None,
                                backend: str = "process"
                                ) -> Dict[str, Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]]]:
        """Parallel recombination tract detection across GUIDs.

        Uses a process pool when fork is available, otherwise falls back to threads.
        """
        if not guids:
            guids = self._headerUtils.headers["guids"]
        self._headerUtils.check_guids(guids)
        if not anchor_pops or not isinstance(anchor_pops, dict):
            raise ParameterError("anchor_pops must be a non-empty dict of {population: [guids]}.")
        if expected_tract_len < 1:
            raise ParameterError("expected_tract_len must be at least 1.")
        if min_tract_sites < 1:
            raise ParameterError("min_tract_sites must be at least 1.")
        if pseudocount < 0:
            raise ParameterError("pseudocount must be non-negative.")

        if allele_id_format is None:
            allele_id_format = self._headerUtils._allele_id_format
        if allele_id_format is None:
            raise ParameterError("allele_id_format is None and cannot inherit from _headerUtils. Please provide an allele_id_format.")

        backend_lower = backend.lower()
        if backend_lower not in {"process", "thread"}:
            raise ParameterError("backend must be either 'process' or 'thread'.")

        ordered_indices, ordered_pos, ordered_alleles, chr_slices = self._build_ordered_snv_loci(
            allele_id_format=allele_id_format
        )

        allele_idx_arr = np.asarray(ordered_indices, dtype=np.uint64)
        n_loci = allele_idx_arr.size
        pop_names = list(anchor_pops.keys())
        if not pop_names:
            raise ParameterError("anchor_pops must define at least one population.")

        log.info(f"Preparing recombination model with {len(pop_names)} populations and {n_loci} SNV loci.")

        eps = 1e-6
        logp = np.empty((len(pop_names), n_loci), dtype=np.float64)
        log1p = np.empty((len(pop_names), n_loci), dtype=np.float64)

        for i, pop in enumerate(pop_names):
            pop_guids = anchor_pops.get(pop, [])
            if not pop_guids:
                raise ParameterError(f"anchor_pops[{pop}] cannot be empty.")
            self._headerUtils.check_guids(pop_guids)
            pop_coords = [self._headerUtils.encode_guid(guid) for guid in pop_guids]
            freqs = np.asarray(self._matrix.colFrequency(pop_coords), dtype=np.float64)
            freqs = freqs[allele_idx_arr]

            n = float(len(pop_guids))
            p = (freqs * n + pseudocount) / (n + 2.0 * pseudocount)
            p = np.clip(p, eps, 1.0 - eps)
            logp[i, :] = np.log(p)
            log1p[i, :] = np.log1p(-p)

        if len(pop_names) == 1:
            log_trans = np.zeros((1, 1), dtype=np.float64)
            log_start = np.zeros(1, dtype=np.float64)
        else:
            p_switch = 1.0 / float(expected_tract_len)
            p_switch = min(max(p_switch, eps), 1.0 - eps)
            stay = 1.0 - p_switch
            switch = p_switch / float(len(pop_names) - 1)
            log_trans = np.full((len(pop_names), len(pop_names)), np.log(switch), dtype=np.float64)
            np.fill_diagonal(log_trans, np.log(stay))
            log_start = np.full(len(pop_names), -np.log(len(pop_names)), dtype=np.float64)

        model = {
            "allele_idx_arr": allele_idx_arr,
            "ordered_pos": ordered_pos,
            "ordered_alleles": ordered_alleles,
            "chr_slices": chr_slices,
            "logp": logp,
            "log1p": log1p,
            "log_trans": log_trans,
            "log_start": log_start,
            "pop_names": pop_names,
        }

        global _RECOMB_SELF, _RECOMB_MODEL
        _RECOMB_SELF = self
        _RECOMB_MODEL = model

        results: Dict[str, Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]]] = {}

        if backend_lower == "process":
            if "fork" not in mp.get_all_start_methods():
                log.warning("Process backend requires fork; falling back to threads.")
                backend_lower = "thread"

        if backend_lower == "process":
            ctx = mp.get_context("fork")
            with ProcessPoolExecutor(mp_context=ctx, max_workers=max_workers) as executor:
                futures = [
                    executor.submit(_recomb_worker, guid, min_tract_sites, return_loci, threads)
                    for guid in guids
                ]
                for future in as_completed(futures):
                    guid, out = future.result()
                    results[guid] = out
        else:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(_recomb_worker, guid, min_tract_sites, return_loci, threads)
                    for guid in guids
                ]
                for future in as_completed(futures):
                    guid, out = future.result()
                    results[guid] = out

        return results
