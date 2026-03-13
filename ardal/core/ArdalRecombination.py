""" ArdalRecombination.py
This module provides functionality for investigation signals of recombination in the Ardal framework.
"""
import logging
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
        raise RuntimeError("Ancestry worker is not initialised.")

    ## extract observations for this guid
    obs, missing_mask = recomb._extract_guid_observations(
        guid=guid,
        allele_indices=model["allele_idx_arr"],
        threads=threads,
    )

    out: Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]] = {"windows": []}
    if return_loci:
        out["loci"] = []

    ## process each chromosome independently so tract calls cannot cross boundaries
    for chr_key, (start, end) in model["chr_slices"].items():
        if start == end:
            continue
        obs_slice = obs[start:end]
        missing_slice = missing_mask[start:end]
        ## emission log likelihoods per population
        emit = np.where(obs_slice, model["logp"][:, start:end], model["log1p"][:, start:end])
        if missing_slice.any():
            emit[:, missing_slice] = 0.0

        ## run viterbi to assign populations
        states = recomb._viterbi(emit, model["log_trans"], model["log_start"])
        state_post = recomb._forward_backward_posteriors(emit, model["log_trans"], model["log_start"])

        if return_loci:
            ## expand per locus assignments
            for offset, state in enumerate(states):
                out["loci"].append({
                    "chr": chr_key,
                    "pos": model["ordered_pos"][start + offset],
                    "allele": model["ordered_alleles"][start + offset],
                    "population": model["pop_names"][int(state)],
                })

        if states.size == 0:
            continue

        ## segment into runs of identical states
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

        ## score each segment and summarise chromosome-local state marginals
        for seg_start, seg_end, state in segments:
            abs_start = start + seg_start
            abs_end = start + seg_end
            window_scores = emit[:, seg_start:seg_end + 1].sum(axis=1)
            posterior = state_post[:, seg_start:seg_end + 1].mean(axis=1).tolist()
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
    """ Class for investigation signals of recombination and ancestry tracts in the Ardal framework.
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
        ## collect valid snv loci and preserve allele indices
        allele_to_locus, _, _ = self._headerUtils.get_snv_lookup(allele_id_format=allele_id_format)
        invalid = np.uint32(np.iinfo(np.uint32).max)
        valid_indices = np.flatnonzero(allele_to_locus != invalid)

        if valid_indices.size == 0:
            raise ParameterError("No biallelic SNV alleles were resolved for tract detection.")

        self._headerUtils.ensure_id_positions(allele_id_format)
        allele_positions = self._headerUtils.get_allele_positions()
        alleles = self._headerUtils.headers["alleles"]

        ## group alleles by chromosome
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

        ## build ordered loci and per chromosome slices
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
        ## load packed row and unpack alleles
        guid_idx = self._headerUtils.encode_guid(guid)
        packed = self._matrix.getSubsetPackedMatrix_rows([guid_idx], threads, True)
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

        ## apply missing mask if present
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
        ## dynamic programming for most likely states
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
    def _logsumexp(arr: np.ndarray,
                   axis: int) -> np.ndarray:
        max_arr = np.max(arr, axis=axis, keepdims=True)
        stable = np.exp(arr - max_arr)
        summed = np.sum(stable, axis=axis, keepdims=True)
        return np.squeeze(np.log(summed) + max_arr, axis=axis)


    @staticmethod
    def _forward_backward_posteriors(emit: np.ndarray,
                                     log_trans: np.ndarray,
                                     log_start: np.ndarray) -> np.ndarray:
        """Posterior state probabilities per locus from a chromosome-local HMM."""
        n_states, n_sites = emit.shape
        if n_sites == 0:
            return np.zeros((n_states, 0), dtype=np.float64)

        alpha = np.empty((n_states, n_sites), dtype=np.float64)
        beta = np.empty((n_states, n_sites), dtype=np.float64)

        alpha[:, 0] = log_start + emit[:, 0]
        for t in range(1, n_sites):
            scores = alpha[:, t - 1][:, None] + log_trans
            alpha[:, t] = emit[:, t] + ArdalRecombination._logsumexp(scores, axis=0)

        beta[:, -1] = 0.0
        for t in range(n_sites - 2, -1, -1):
            scores = log_trans + emit[:, t + 1][None, :] + beta[:, t + 1][None, :]
            beta[:, t] = ArdalRecombination._logsumexp(scores, axis=1)

        log_norm = float(ArdalRecombination._logsumexp(alpha[:, -1], axis=0))
        log_post = alpha + beta - log_norm
        post = np.exp(log_post)
        post_sum = post.sum(axis=0, keepdims=True)
        post_sum[post_sum == 0.0] = 1.0
        post /= post_sum
        return post


    @staticmethod
    def _merge_short_segments(segments: List[List[int]],
                              min_sites: int) -> List[List[int]]:
        ## merge short segments into neighbors
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
        ## collapse adjacent segments with same state
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


    @staticmethod
    def _build_cooc_weights(encoded_indices: List[int],
                            cooc_dict: Dict,
                            n_loci: int) -> Tuple[np.ndarray, int]:
        """Build per-locus weights from parallel_backend co-occurrence edges using union-find."""
        parent = list(range(n_loci))
        comp_size = [1] * n_loci
        idx_map = {int(encoded_idx): idx for idx, encoded_idx in enumerate(encoded_indices)}

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: int, b: int) -> None:
            ra = find(a)
            rb = find(b)
            if ra == rb:
                return
            if comp_size[ra] < comp_size[rb]:
                ra, rb = rb, ra
            parent[rb] = ra
            comp_size[ra] += comp_size[rb]

        for ref, neighbors in cooc_dict.items():
            ref_idx = idx_map.get(int(ref))
            if ref_idx is None:
                continue
            for neighbor in neighbors:
                neighbor_idx = idx_map.get(int(neighbor))
                if neighbor_idx is not None:
                    union(ref_idx, neighbor_idx)

        roots = [find(idx) for idx in range(n_loci)]
        block_roots = set(roots)
        weights = np.empty(n_loci, dtype=np.float64)
        for idx, root in enumerate(roots):
            weights[idx] = 1.0 / float(comp_size[root])

        return weights, len(block_roots)


    @staticmethod
    def _build_cooc_partitions(chr_slices: Dict[str, Tuple[int, int]],
                               ordered_pos: List[int],
                               mode: str,
                               window_size: Union[int, None]) -> List[Tuple[str, int, int]]:
        """Construct global, chromosome, or genomic-window co-occurrence partitions."""
        if mode == "global":
            return [("global", 0, len(ordered_pos))]

        partitions: List[Tuple[str, int, int]] = []
        for chr_key, (start, end) in chr_slices.items():
            if start >= end:
                continue
            if mode == "chromosome":
                partitions.append((f"chr={chr_key}", start, end))
                continue

            if window_size is None or window_size < 1:
                raise ParameterError("cooc_window_size must be at least 1 when cooc_partition='window'.")

            window_start = start
            while window_start < end:
                window_start_pos = ordered_pos[window_start]
                limit = window_start_pos + window_size
                window_end = window_start + 1
                while window_end < end and ordered_pos[window_end] < limit:
                    window_end += 1
                label = f"chr={chr_key}:{window_start_pos}-{ordered_pos[window_end - 1]}"
                partitions.append((label, window_start, window_end))
                window_start = window_end

        return partitions


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
                       threads : int = 1,
                       max_workers: Union[int, None] = None,
                       parallel_backend: str = "auto",
                       cooc_enabled: bool = False,
                       cooc_threshold: float = 0.90,
                       cooc_partition: str = "global",
                       cooc_window_size: Union[int, None] = None ) -> Dict[str, Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]]]:
        """Detect ancestry tracts by assigning lineage along ordered SNV loci.

        Args:
            guids: Query GUIDs to segment.
            anchor_pops: Mapping of population name -> list of GUIDs.
            allele_id_format: Format string used to decode SNV loci.
            expected_tract_len: Expected tract length in sites (controls switch penalty).
            min_tract_sites: Minimum tract length in sites for reporting.
            pseudocount: Pseudocount for allele frequency smoothing.
            return_loci: When True, return per-locus assignments alongside tracts.
            threads: Thread count for matrix extraction.
            max_workers: Worker count for parallel mode.
            parallel_backend: auto, process, or thread.
            cooc_enabled: When True apply cooc weights to reduce inflation.
            cooc_threshold: Threshold passed to allele cooc computation.
            cooc_partition: Partitioning mode for co-occurrence weighting: global, chromosome, or window.
            cooc_window_size: Genomic window size used when cooc_partition='window'.

        Returns:
            Dict keyed by GUID, each containing "windows" entries (and optional "loci")
            with per-window population assignments, per-population log-likelihoods,
            and tract posteriors derived from chromosome-local state marginals.
        """
        ## validate inputs and defaults
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

        parallel_backend_lower = parallel_backend.lower()
        if parallel_backend_lower not in {"auto", "process", "thread"}:
            raise ParameterError("parallel_backend must be 'auto', 'process', or 'thread'.")
        cooc_partition_lower = cooc_partition.lower()
        if cooc_partition_lower not in {"global", "chromosome", "window"}:
            raise ParameterError("cooc_partition must be 'global', 'chromosome', or 'window'.")
        if cooc_window_size is not None and cooc_window_size < 1:
            raise ParameterError("cooc_window_size must be at least 1.")

        ## build ordered snv loci
        ordered_indices, ordered_pos, ordered_alleles, chr_slices = self._build_ordered_snv_loci(
            allele_id_format=allele_id_format
        )

        allele_idx_arr = np.asarray(ordered_indices, dtype=np.uint64)
        n_loci = allele_idx_arr.size
        pop_names = list(anchor_pops.keys())
        if not pop_names:
            raise ParameterError("anchor_pops must define at least one population.")

        log.info(f"Preparing ancestry model with {len(pop_names)} populations and {n_loci} SNV loci.")

        eps = 1e-6
        logp = np.empty((len(pop_names), n_loci), dtype=np.float64)
        log1p = np.empty((len(pop_names), n_loci), dtype=np.float64)

        ## build emission probabilities per population
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

        ## build cooc weights once and fold them into emissions up front
        cooc_weights = None
        if cooc_enabled:
            partitions = self._build_cooc_partitions(
                chr_slices=chr_slices,
                ordered_pos=ordered_pos,
                mode=cooc_partition_lower,
                window_size=cooc_window_size,
            )
            cooc_weights = np.ones(n_loci, dtype=np.float64)
            total_blocks = 0
            log.info(
                f"Preparing cooc weights using mode='{cooc_partition_lower}' across {len(partitions)} partitions."
            )
            for part_idx, (label, start, end) in enumerate(partitions, start=1):
                part_size = end - start
                if part_size <= 1:
                    log.debug(
                        f"[COOC {part_idx}/{len(partitions)}] Skipping {label}; only {part_size} locus available."
                    )
                    continue

                log.debug(
                    f"[COOC {part_idx}/{len(partitions)}] Starting {label} with {part_size} loci."
                )
                part_indices = ordered_indices[start:end]
                cooc_dict = self._matrix.bitCooccurrence_subset(col_indices=part_indices,
                                                                threshold=cooc_threshold,
                                                                threads=threads)
                part_weights, block_count = self._build_cooc_weights(
                    encoded_indices=part_indices,
                    cooc_dict=cooc_dict,
                    n_loci=part_size,
                )
                cooc_weights[start:end] = part_weights
                total_blocks += block_count

                if log.isEnabledFor(logging.DEBUG):
                    edge_count = sum(len(neighbors) for neighbors in cooc_dict.values())
                    log.debug(
                        f"[COOC {part_idx}/{len(partitions)}] Finished {label}; "
                        f"refs={len(cooc_dict)}, edges={edge_count}, blocks={block_count}."
                    )

            weight_view = cooc_weights[np.newaxis, :]
            logp *= weight_view
            log1p *= weight_view
            log.info(
                f"Applied cooc weights across {len(partitions)} partitions with {total_blocks} blocks."
            )

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
            "cooc_weights": cooc_weights,
        }

        global _RECOMB_SELF, _RECOMB_MODEL
        _RECOMB_SELF = self
        _RECOMB_MODEL = model

        ## choose parallel_backend
        use_pool = len(guids) > 1 and (max_workers is None or max_workers > 1)
        if not use_pool:
            results: Dict[str, Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]]] = {}
            for guid in guids:
                guid_key, out = _recomb_worker(guid, min_tract_sites, return_loci, threads)
                results[guid_key] = out
            return results

        if parallel_backend_lower == "auto":
            parallel_backend_lower = "process" if "fork" in mp.get_all_start_methods() else "thread"

        if parallel_backend_lower == "process" and "fork" not in mp.get_all_start_methods():
            log.warning("Process parallel_backend requires fork; falling back to threads.")
            parallel_backend_lower = "thread"

        results: Dict[str, Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]]] = {}
        if parallel_backend_lower == "process":
            ctx = mp.get_context("fork")
            with ProcessPoolExecutor(mp_context=ctx, max_workers=max_workers) as executor:
                futures = [
                    executor.submit(_recomb_worker, guid, min_tract_sites, return_loci, threads)
                    for guid in guids
                ]
                for future in as_completed(futures):
                    guid_key, out = future.result()
                    results[guid_key] = out
        else:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(_recomb_worker, guid, min_tract_sites, return_loci, threads)
                    for guid in guids
                ]
                for future in as_completed(futures):
                    guid_key, out = future.result()
                    results[guid_key] = out

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
                                parallel_backend: str = "process",
                                cooc_enabled: bool = False,
                                cooc_threshold: float = 0.95,
                                cooc_partition: str = "global",
                                cooc_window_size: Union[int, None] = None
                                ) -> Dict[str, Dict[str, List[Dict[str, Union[str, int, float, Dict[str, float]]]]]]:
        """LEGACY FUNCTION: Parallel recombination tract detection across GUIDs.

        This is a compatibility wrapper for detect_tracts with parallel_backend control.
        """
        ## passthrough to unified detector
        return self.detect_tracts(
            anchor_pops=anchor_pops,
            allele_id_format=allele_id_format,
            guids=guids,
            expected_tract_len=expected_tract_len,
            min_tract_sites=min_tract_sites,
            pseudocount=pseudocount,
            return_loci=return_loci,
            threads=threads,
            max_workers=max_workers,
            parallel_backend=parallel_backend,
            cooc_enabled=cooc_enabled,
            cooc_threshold=cooc_threshold,
            cooc_partition=cooc_partition,
            cooc_window_size=cooc_window_size,
        )
