/*
Copyright 2025 Arthur V. Morris
*/
#include "SnvFallbackBackend.hpp"
#include "utils/PythonLogger.hpp"
#include <algorithm>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils/bitops.hpp"

namespace py = pybind11;

namespace ardal {

namespace {

inline uint32_t compute_snv_distance( const std::vector<uint64_t>& lhs,
                                      const std::vector<uint64_t>& rhs ) {
    std::size_t p = 0;
    std::size_t q = 0;
    uint32_t dist = 0;

    while (p < lhs.size() && q < rhs.size()) {
        const uint64_t li = lhs[p];
        const uint64_t ri = rhs[q];
        const uint32_t locus_l = static_cast<uint32_t>(li >> 3);
        const uint32_t locus_r = static_cast<uint32_t>(ri >> 3);

        if (locus_l == locus_r) {
            const uint32_t base_l = static_cast<uint32_t>(li & 0x7u);
            const uint32_t base_r = static_cast<uint32_t>(ri & 0x7u);
            if (base_l != base_r) ++dist;
            ++p;
            ++q;
        } else if (locus_l < locus_r) {
            ++dist;
            ++p;
        } else {
            ++dist;
            ++q;
        }
    }

    dist += static_cast<uint32_t>((lhs.size() - p) + (rhs.size() - q));
    return dist;
}


inline uint32_t compute_snv_distance_filtered( const std::vector<uint64_t>& lhs,
                                               const std::vector<uint64_t>& rhs,
                                               const std::vector<uint8_t>& locus_mask ) {
    constexpr uint32_t INVALID = std::numeric_limits<uint32_t>::max();
    auto advance = [&](const std::vector<uint64_t>& vec, std::size_t& pos) -> uint32_t {
        while (pos < vec.size()) {
            const uint32_t locus = static_cast<uint32_t>(vec[pos] >> 3);
            if (locus < locus_mask.size() && locus_mask[locus]) {
                return locus;
            }
            ++pos;
        }
        return INVALID;
    };

    auto remaining = [&](const std::vector<uint64_t>& vec, std::size_t pos) -> uint32_t {
        uint32_t count = 0;
        while (pos < vec.size()) {
            const uint32_t locus = static_cast<uint32_t>(vec[pos] >> 3);
            if (locus < locus_mask.size() && locus_mask[locus]) {
                ++count;
            }
            ++pos;
        }
        return count;
    };

    std::size_t p = 0;
    std::size_t q = 0;
    uint32_t dist = 0;

    while (true) {
        const uint32_t locus_l = advance(lhs, p);
        const uint32_t locus_r = advance(rhs, q);
        const bool lhs_done = (locus_l == INVALID);
        const bool rhs_done = (locus_r == INVALID);

        if (lhs_done || rhs_done) {
            if (!lhs_done) dist += remaining(lhs, p);
            if (!rhs_done) dist += remaining(rhs, q);
            break;
        }

        if (locus_l == locus_r) {
            const uint32_t base_l = static_cast<uint32_t>(lhs[p] & 0x7u);
            const uint32_t base_r = static_cast<uint32_t>(rhs[q] & 0x7u);
            if (base_l != base_r) ++dist;
            ++p;
            ++q;
        } else if (locus_l < locus_r) {
            ++dist;
            ++p;
        } else {
            ++dist;
            ++q;
        }
    }

    return dist;
}

} // namespace

SnvFallbackBackend::SnvFallbackBackend( const ardal::detail::FlatMatrix& flat,
                                        std::size_t n_cols_bits,
                                        const std::vector<std::vector<uint32_t>>* missing_mask )
  : flat_(flat),
    n_rows_(flat.n_rows),
    n_cols_bits_(n_cols_bits),
    words_per_row_(flat.wpr),
    missing_mask_(missing_mask),
    has_missing_mask_(missing_mask && !missing_mask->empty()) {}


void SnvFallbackBackend::prepare( py::array_t<uint32_t> allele_to_locus,
                                  py::array_t<uint8_t> allele_to_base ) {
    py::buffer_info locus_buf = allele_to_locus.request();
    py::buffer_info base_buf = allele_to_base.request();

    if (static_cast<std::size_t>(locus_buf.size) != n_cols_bits_) {
        throw std::invalid_argument("allele_to_locus length does not match number of columns");
    }
    if (static_cast<std::size_t>(base_buf.size) != n_cols_bits_) {
        throw std::invalid_argument("allele_to_base length does not match number of columns");
    }

    const auto* locus_ptr = static_cast<const uint32_t*>(locus_buf.ptr);
    const auto* base_ptr = static_cast<const uint8_t*>(base_buf.ptr);

    allele_to_locus_.assign(locus_ptr, locus_ptr + locus_buf.size);
    allele_to_base_.assign(base_ptr, base_ptr + base_buf.size);

    snv_vectors_.clear();
    vectors_ready_ = false;
    snv_entries_.clear();
    entries_ready_ = false;
    lookup_loaded_ = true;
}


/****************************************************************************************************
 * ardal::SnvFallbackBackend::snvHamming
 *
 * Compute the condensed SNV Hamming distance matrix across all rows.
 *
 * INPUT:
 *  threads (int) : Number of OpenMP threads to use.
 *
 * OUTPUT:
 *  py::array_t<uint32_t> : Condensed lower-triangular distance matrix.
 ****************************************************************************************************/
py::array_t<uint32_t> SnvFallbackBackend::snvHamming( int threads, bool mask_missing ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (n_rows_ <= 1) {
        return py::array_t<uint32_t>(static_cast<py::ssize_t>(0));
    }
    if (!lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    ensure_vectors();

    const std::size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<uint32_t> dist_condensed(static_cast<py::ssize_t>(total_pairs));
    auto* dist_ptr = dist_condensed.mutable_data();

    ardal::utils::log_info("Running SnvFallbackBackend::snvHamming.");
    {
        py::gil_scoped_release release;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (std::size_t i = 0; i < n_rows_; ++i) {
            const std::size_t idx_base = (i * (2 * n_rows_ - i - 1)) / 2;
            for (std::size_t j = i + 1; j < n_rows_; ++j) {
                dist_ptr[idx_base + (j - i - 1)] = (mask_missing && has_missing_mask_)
                    ? snvDistanceMasked(i, j, nullptr)
                    : snvDistance(i, j);
            }
        }
    }

    return dist_condensed;
}


/****************************************************************************************************
 * ardal::SnvFallbackBackend::snvHamming_subset
 *
 * Compute an SNV Hamming distance matrix for a subset of rows and columns.
 *
 * INPUT:
 *  row_indices (const std::vector<size_t>&) : Rows to include.
 *  col_indices (const std::vector<size_t>&) : Columns (alleles) to include.
 *
 * PARAMETERS:
 *  threads (int) : Number of OpenMP threads to use.
 *
 * OUTPUT:
 *  py::array_t<uint32_t> : Condensed lower-triangular distance matrix.
 ****************************************************************************************************/
py::array_t<uint32_t> SnvFallbackBackend::snvHamming_subset( const std::vector<size_t>& row_indices,
                                                            const std::vector<size_t>& col_indices,
                                                            int threads,
                                                            bool mask_missing ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (!lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }
    ensure_vectors();

    const size_t sub_rows = row_indices.size();
    if (sub_rows <= 1) {
        return py::array_t<uint32_t>(static_cast<py::ssize_t>(0));
    }

    for (size_t row_idx : row_indices) {
        if (row_idx >= n_rows_) {
            throw std::out_of_range("Row index in row_indices is out of bounds.");
        }
    }

    const size_t total_pairs = sub_rows * (sub_rows - 1) / 2;
    if (col_indices.empty()) {
        py::array_t<uint32_t> zeros(static_cast<py::ssize_t>(total_pairs));
        auto* zero_ptr = zeros.mutable_data();
        std::fill(zero_ptr, zero_ptr + total_pairs, static_cast<uint32_t>(0));
        return zeros;
    }

    const bool use_mask = !isFullColumnSelection(col_indices);
    std::vector<uint8_t> locus_mask;
    if (use_mask) {
        locus_mask = buildLocusMask(col_indices);
        if (locus_mask.empty()) {
            py::array_t<uint32_t> zeros(static_cast<py::ssize_t>(total_pairs));
            auto* zero_ptr = zeros.mutable_data();
            std::fill(zero_ptr, zero_ptr + total_pairs, static_cast<uint32_t>(0));
            return zeros;
        }
    }

    py::array_t<uint32_t> dist_matrix(static_cast<py::ssize_t>(total_pairs));
    auto* dist_ptr = dist_matrix.mutable_data();

    {
        py::gil_scoped_release release;

        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (size_t i = 0; i < sub_rows; ++i) {
            const size_t lhs_idx = row_indices[i];
            const size_t idx_base = (i * (2 * sub_rows - i - 1)) / 2;
            for (size_t j = i + 1; j < sub_rows; ++j) {
                const size_t rhs_idx = row_indices[j];
                uint32_t dist_val;
                if (mask_missing && has_missing_mask_) {
                    dist_val = snvDistanceMasked(lhs_idx, rhs_idx, use_mask ? &locus_mask : nullptr);
                } else if (use_mask) {
                    dist_val = compute_snv_distance_filtered(snv_vectors_[lhs_idx], snv_vectors_[rhs_idx], locus_mask);
                } else {
                    dist_val = snvDistance(lhs_idx, rhs_idx);
                }
                dist_ptr[idx_base + (j - i - 1)] = dist_val;
            }
        }
    }

    return dist_matrix;
}


/****************************************************************************************************
 * ardal::SnvFallbackBackend::snvNeighbourhood
 *
 * Enumerate rows within an SNV Hamming radius of the query row.
 *
 * INPUT:
 *  row_idx (size_t) : Query row index.
 *  epsilon (uint32_t) : Maximum SNV Hamming distance.
 *
 * PARAMETERS:
 *  threads (int) : Number of OpenMP threads to use.
 *
 * OUTPUT:
 *  py::array_t<int64_t> : Array of (row, distance) pairs.
 ****************************************************************************************************/
py::array_t<int64_t> SnvFallbackBackend::snvNeighbourhood( std::size_t row_idx,
                                                           uint32_t epsilon,
                                                           int threads,
                                                           bool mask_missing ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (row_idx >= n_rows_) {
        throw std::runtime_error("Row index out of range.");
    }
    if (!lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    ensure_vectors();

    const int T = std::max(1, threads);
    std::vector<std::vector<std::pair<std::size_t, uint32_t>>> buckets(T);

    ardal::utils::log_info("Running SnvFallbackBackend::snvNeighbourhood.");
    {
        py::gil_scoped_release release;

#if defined(_OPENMP)
        #pragma omp parallel num_threads(T)
#endif
        {
#if defined(_OPENMP)
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            auto& local = buckets[tid];

#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for (std::size_t i = 0; i < n_rows_; ++i) {
                if (i == row_idx) continue;
                const uint32_t dist = (mask_missing && has_missing_mask_)
                    ? snvDistanceMasked(row_idx, i, nullptr)
                    : snvDistance(row_idx, i);
                if (dist <= epsilon) {
                    local.emplace_back(i, dist);
                }
            }
        }
    }

    std::size_t total = 0;
    for (const auto& vec : buckets) {
        total += vec.size();
    }

    py::array_t<int64_t> result({total, static_cast<std::size_t>(2)});
    auto* data = result.mutable_data();

    std::size_t idx = 0;
    for (const auto& vec : buckets) {
        for (const auto& entry : vec) {
            data[idx * 2] = static_cast<int64_t>(entry.first);
            data[idx * 2 + 1] = static_cast<int64_t>(entry.second);
            ++idx;
        }
    }

    return result;
}


/****************************************************************************************************
 * ardal::SnvFallbackBackend::knnSnv
 *
 * Find the k nearest neighbours in SNV Hamming space.
 *
 * INPUT:
 *  row_idx (size_t) : Query row index.
 *  k (uint32_t) : Number of neighbours to return.
 *
 * PARAMETERS:
 *  threads (int) : Number of OpenMP threads to use.
 *
 * OUTPUT:
 *  py::list : List of (row, distance) tuples.
 ****************************************************************************************************/
py::list SnvFallbackBackend::knnSnv( std::size_t row_idx,
                                     uint32_t k,
                                     int threads,
                                     bool mask_missing ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (row_idx >= n_rows_) {
        throw std::runtime_error("Row index out of range.");
    }
    if (!lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    ensure_vectors();

    if (n_rows_ <= 1 || k == 0) {
        return py::list();
    }
    const uint32_t k_eff = std::min<uint32_t>(k, n_rows_ > 1 ? static_cast<uint32_t>(n_rows_ - 1) : 0);
    if (k_eff == 0) {
        return py::list();
    }

    struct Neighbor {
        uint32_t id;
        uint32_t dist;
    };
    struct ByMaxDistance {
        bool operator()(const Neighbor& a, const Neighbor& b) const noexcept {
            if (a.dist == b.dist) return a.id < b.id;
            return a.dist < b.dist;
        }
    };
    auto asc_dist_id = [](const Neighbor& a, const Neighbor& b) noexcept {
        return (a.dist < b.dist) || (a.dist == b.dist && a.id < b.id);
    };

    std::vector<Neighbor> final_neighbors;
    {
        py::gil_scoped_release release;

        const int T = std::max(1, threads);
        std::vector<std::vector<Neighbor>> buckets(T);

#if defined(_OPENMP)
        #pragma omp parallel num_threads(T)
#endif
        {
#if defined(_OPENMP)
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            auto& local = buckets[tid];
            std::priority_queue<Neighbor, std::vector<Neighbor>, ByMaxDistance> heap;

#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for (uint32_t j = 0; j < n_rows_; ++j) {
                if (j == row_idx) continue;

                const uint32_t dist = (mask_missing && has_missing_mask_)
                    ? snvDistanceMasked(row_idx, j, nullptr)
                    : snvDistance(row_idx, j);

                if (heap.size() < k_eff) {
                    heap.push(Neighbor{j, dist});
                } else {
                    const Neighbor& top = heap.top();
                    if (dist < top.dist || (dist == top.dist && j < top.id)) {
                        heap.pop();
                        heap.push(Neighbor{j, dist});
                    }
                }
            }

            local.reserve(heap.size());
            while (!heap.empty()) {
                local.push_back(heap.top());
                heap.pop();
            }
            std::sort(local.begin(), local.end(), asc_dist_id);
        }

        std::size_t total = 0;
        for (const auto& vec : buckets) {
            total += vec.size();
        }
        std::vector<Neighbor> all;
        all.reserve(total);
        for (auto& vec : buckets) {
            all.insert(all.end(), vec.begin(), vec.end());
        }

        if (all.size() > k_eff) {
            std::nth_element(all.begin(), all.begin() + k_eff, all.end(), asc_dist_id);
            all.resize(k_eff);
        }
        std::sort(all.begin(), all.end(), asc_dist_id);
        final_neighbors = std::move(all);
    }

    py::list result;
    for (const auto& nb : final_neighbors) {
        result.append(py::make_tuple(static_cast<int64_t>(nb.id),
                                     static_cast<int64_t>(nb.dist)));
    }
    return result;
}


void SnvFallbackBackend::ensure_vectors() const {
    if (vectors_ready_ && (!has_missing_mask_ || entries_ready_)) return;
    if (!lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    snv_vectors_.assign(n_rows_, {});
    if (has_missing_mask_) {
        snv_entries_.assign(n_rows_, {});
    }

    const uint64_t* base_ptr = flat_.base;
    const bool has_tail = (n_cols_bits_ & 63u) != 0u;
    const uint64_t tail_mask = has_tail ? ardal::tail_mask(n_cols_bits_) : ~0ULL;
    const uint32_t invalid_locus = std::numeric_limits<uint32_t>::max();

    for (std::size_t row = 0; row < n_rows_; ++row) {
        std::vector<std::tuple<uint32_t, uint8_t, uint32_t>> collected;

        for (std::size_t w = 0; w < words_per_row_; ++w) {
            uint64_t word = base_ptr[row * words_per_row_ + w];
            if (has_tail && (w + 1 == words_per_row_)) {
                word &= tail_mask;
            }

            while (word) {
                const unsigned bit = ardal::ctz64_nonzero(word);
                const std::size_t col = (w << 6) | bit;
                word &= (word - 1);
                if (col >= n_cols_bits_) continue;

                const uint32_t locus_id = allele_to_locus_[col];
                if (locus_id == invalid_locus) continue;

                const uint8_t base_code = allele_to_base_[col];
                if (base_code == 0) continue;

                collected.emplace_back(locus_id, base_code, static_cast<uint32_t>(col));
            }
        }

        if (collected.empty()) continue;

        std::sort(collected.begin(), collected.end(),
                  [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });

        std::vector<uint64_t> encoded;
        encoded.reserve(collected.size());
        if (has_missing_mask_) {
            snv_entries_[row].reserve(collected.size());
        }

        std::size_t idx = 0;
        while (idx < collected.size()) {
            std::size_t j = idx + 1;
            bool conflict = false;
            const uint32_t locus_id = std::get<0>(collected[idx]);
            const uint8_t base_code = std::get<1>(collected[idx]);
            const uint32_t col_id = std::get<2>(collected[idx]);
            while (j < collected.size() && std::get<0>(collected[j]) == locus_id) {
                if (std::get<1>(collected[j]) != base_code) {
                    conflict = true;
                }
                ++j;
            }

            if (!conflict) {
                encoded.push_back((static_cast<uint64_t>(locus_id) << 3) |
                                  static_cast<uint64_t>(base_code & 0x7u));
                if (has_missing_mask_) {
                    snv_entries_[row].push_back(SnvEntry{col_id, locus_id, base_code});
                }
            }

            idx = j;
        }

        snv_vectors_[row] = std::move(encoded);
    }

    vectors_ready_ = true;
    if (has_missing_mask_) {
        entries_ready_ = true;
    }
}


uint32_t SnvFallbackBackend::snvDistance( std::size_t i, std::size_t j ) const {
    const auto& lhs = snv_vectors_[i];
    const auto& rhs = snv_vectors_[j];
    return compute_snv_distance(lhs, rhs);
}


bool SnvFallbackBackend::is_missing( const std::vector<uint32_t>* mask, std::size_t col ) const {
    if (!mask || mask->empty()) {
        return false;
    }
    return std::binary_search(mask->begin(), mask->end(), static_cast<uint32_t>(col));
}


bool SnvFallbackBackend::locus_allowed( uint32_t locus, const std::vector<uint8_t>* locus_mask ) const {
    if (!locus_mask || locus_mask->empty()) {
        return true;
    }
    return (locus < locus_mask->size()) && ((*locus_mask)[locus] != 0);
}


void SnvFallbackBackend::collect_encoded_masked( std::size_t row_idx,
                                                 const std::vector<uint32_t>* other_mask,
                                                 const std::vector<uint8_t>* locus_mask,
                                                 std::vector<uint64_t>& out ) const {
    out.clear();
    const auto& entries = snv_entries_[row_idx];
    const std::vector<uint32_t>* self_mask = (has_missing_mask_ && missing_mask_) ? &(*missing_mask_)[row_idx] : nullptr;

    for (const auto& e : entries) {
        if (is_missing(self_mask, e.col) || is_missing(other_mask, e.col)) {
            continue;
        }
        if (!locus_allowed(e.locus, locus_mask)) {
            continue;
        }
        out.push_back((static_cast<uint64_t>(e.locus) << 3) |
                      static_cast<uint64_t>(e.base & 0x7u));
    }
}


uint32_t SnvFallbackBackend::snvDistanceMasked( std::size_t i,
                                                std::size_t j,
                                                const std::vector<uint8_t>* locus_mask ) const {
    if (!has_missing_mask_) {
        return snvDistance(i, j);
    }
    if (!entries_ready_) {
        ensure_vectors();
    }
    const std::vector<uint32_t>* mask_i = (has_missing_mask_ && missing_mask_) ? &(*missing_mask_)[i] : nullptr;
    const std::vector<uint32_t>* mask_j = (has_missing_mask_ && missing_mask_) ? &(*missing_mask_)[j] : nullptr;

    std::vector<uint64_t> lhs_filtered;
    std::vector<uint64_t> rhs_filtered;
    collect_encoded_masked(i, mask_j, locus_mask, lhs_filtered);
    collect_encoded_masked(j, mask_i, locus_mask, rhs_filtered);

    return compute_snv_distance(lhs_filtered, rhs_filtered);
}


bool SnvFallbackBackend::isFullColumnSelection( const std::vector<size_t>& col_indices ) const {
    if (col_indices.size() != n_cols_bits_) {
        return false;
    }
    for (size_t idx = 0; idx < col_indices.size(); ++idx) {
        if (col_indices[idx] != idx) {
            return false;
        }
    }
    return true;
}


std::vector<uint8_t> SnvFallbackBackend::buildLocusMask( const std::vector<size_t>& col_indices ) const {
    const uint32_t invalid_locus = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> loci;
    loci.reserve(col_indices.size());
    uint32_t max_locus = 0;
    bool has_locus = false;

    for (size_t col_idx : col_indices) {
        if (col_idx >= n_cols_bits_) {
            throw std::out_of_range("Column index in col_indices is out of bounds.");
        }
        const uint32_t locus = allele_to_locus_[col_idx];
        if (locus == invalid_locus) continue;
        loci.push_back(locus);
        if (!has_locus || locus > max_locus) {
            max_locus = locus;
            has_locus = true;
        }
    }

    if (!has_locus) {
        return {};
    }

    std::vector<uint8_t> mask(static_cast<std::size_t>(max_locus) + 1, 0u);
    for (uint32_t locus : loci) {
        mask[locus] = 1u;
    }
    return mask;
}

} // namespace ardal
