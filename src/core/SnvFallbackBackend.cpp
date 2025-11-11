/*
Copyright 2025 Arthur V. Morris
*/
#include "SnvFallbackBackend.hpp"
#include "utils/PythonLogger.hpp"
#include <algorithm>
#include <limits>
#include <queue>
#include <stdexcept>
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

} // namespace

SnvFallbackBackend::SnvFallbackBackend( const ardal::detail::FlatMatrix& flat,
                                        std::size_t n_cols_bits )
  : flat_(flat),
    n_rows_(flat.n_rows),
    n_cols_bits_(n_cols_bits),
    words_per_row_(flat.wpr) {}


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
py::array_t<uint32_t> SnvFallbackBackend::snvHamming( int threads ) const {
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
                dist_ptr[idx_base + (j - i - 1)] = snvDistance(i, j);
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
                                                            int threads ) const {
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

    const uint32_t invalid_locus = std::numeric_limits<uint32_t>::max();
    std::unordered_set<uint32_t> allowed_loci;
    allowed_loci.reserve(col_indices.size());
    for (size_t col_idx : col_indices) {
        if (col_idx >= n_cols_bits_) {
            throw std::out_of_range("Column index in col_indices is out of bounds.");
        }
        const uint32_t locus = allele_to_locus_[col_idx];
        if (locus != invalid_locus) {
            allowed_loci.insert(locus);
        }
    }

    if (allowed_loci.empty()) {
        const size_t total_pairs = sub_rows * (sub_rows - 1) / 2;
        py::array_t<uint32_t> zeros(static_cast<py::ssize_t>(total_pairs));
        auto* zero_ptr = zeros.mutable_data();
        std::fill(zero_ptr, zero_ptr + total_pairs, static_cast<uint32_t>(0));
        return zeros;
    }

    for (size_t row_idx : row_indices) {
        if (row_idx >= n_rows_) {
            throw std::out_of_range("Row index in row_indices is out of bounds.");
        }
    }

    std::vector<std::vector<uint64_t>> subset_vectors;
    subset_vectors.reserve(row_indices.size());
    for (size_t row_idx : row_indices) {
        const auto& src = snv_vectors_[row_idx];
        std::vector<uint64_t> filtered;
        filtered.reserve(src.size());
        for (uint64_t entry : src) {
            const uint32_t locus = static_cast<uint32_t>(entry >> 3);
            if (allowed_loci.count(locus) != 0) {
                filtered.push_back(entry);
            }
        }
        subset_vectors.push_back(std::move(filtered));
    }

    const size_t total_pairs = sub_rows * (sub_rows - 1) / 2;
    py::array_t<uint32_t> dist_matrix(static_cast<py::ssize_t>(total_pairs));
    auto* dist_ptr = dist_matrix.mutable_data();

    {
        py::gil_scoped_release release;

        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (size_t i = 0; i < sub_rows; ++i) {
            const auto& lhs = subset_vectors[i];
            const size_t idx_base = (i * (2 * sub_rows - i - 1)) / 2;
            for (size_t j = i + 1; j < sub_rows; ++j) {
                const auto& rhs = subset_vectors[j];
                dist_ptr[idx_base + (j - i - 1)] = compute_snv_distance(lhs, rhs);
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
                                                           int threads ) const {
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
                const uint32_t dist = snvDistance(row_idx, i);
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
                                     int threads ) const {
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

                const uint32_t dist = snvDistance(row_idx, j);

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
    if (vectors_ready_) return;
    if (!lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    snv_vectors_.assign(n_rows_, {});

    const uint64_t* base_ptr = flat_.base;
    const bool has_tail = (n_cols_bits_ & 63u) != 0u;
    const uint64_t tail_mask = has_tail ? ardal::tail_mask(n_cols_bits_) : ~0ULL;
    const uint32_t invalid_locus = std::numeric_limits<uint32_t>::max();

    for (std::size_t row = 0; row < n_rows_; ++row) {
        std::vector<std::pair<uint32_t, uint8_t>> collected;

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

                collected.emplace_back(locus_id, base_code);
            }
        }

        if (collected.empty()) continue;

        std::sort(collected.begin(), collected.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        std::vector<uint64_t> encoded;
        encoded.reserve(collected.size());

        std::size_t idx = 0;
        while (idx < collected.size()) {
            std::size_t j = idx + 1;
            bool conflict = false;
            const uint32_t locus_id = collected[idx].first;
            const uint8_t base_code = collected[idx].second;
            while (j < collected.size() && collected[j].first == locus_id) {
                if (collected[j].second != base_code) {
                    conflict = true;
                }
                ++j;
            }

            if (!conflict) {
                encoded.push_back((static_cast<uint64_t>(locus_id) << 3) |
                                  static_cast<uint64_t>(base_code & 0x7u));
            }

            idx = j;
        }

        snv_vectors_[row] = std::move(encoded);
    }

    vectors_ready_ = true;
}


uint32_t SnvFallbackBackend::snvDistance( std::size_t i, std::size_t j ) const {
    const auto& lhs = snv_vectors_[i];
    const auto& rhs = snv_vectors_[j];
    return compute_snv_distance(lhs, rhs);
}

} // namespace ardal
