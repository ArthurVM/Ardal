/*
Copyright 2025 Arthur V. Morris
*/
#include "RoaringMatrix.hpp"
#include "utils/PythonLogger.hpp"
#include <algorithm>
#include <queue>
#include <limits>
#include <unordered_set>
#include <tuple>

namespace py = pybind11;
namespace ardal {



/****************************************************************************************************
 * ardal::RoaringMatrix::RoaringMatrix
 *
 * Constructor for the RoaringMatrix class.
 *
 * This constructor initializes the RoaringMatrix object by taking a 2D NumPy array of uint8_t
 * (representing a binary matrix) and converting each row into a Roaring bitmap. This
 * representation is particularly efficient for sparse binary data.
 *
 * INPUT:
 *  input_matrix (py::array_t<uint8_t>): A 2D NumPy array containing only binary values (0 or 1).
 *
 * OUTPUT:
 *  None (constructor)
 *
 * EXCEPTIONS:
 *  std::runtime_error:
 *   - If the input matrix is not 2-dimensional.
 *   - If the input matrix contains values other than 0 or 1.
 ****************************************************************************************************/
RoaringMatrix::RoaringMatrix( ardal::detail::FlatMatrix matrix_,
                              std::shared_ptr<const std::vector<int>> row_masses,
                              std::shared_ptr<const std::vector<int>> col_masses,
                              const std::vector<std::vector<uint32_t>>* missing_mask,
                              const MissingRanges* missing_ranges )
  : row_masses_(std::move(row_masses)),
    col_masses_(std::move(col_masses)),
    base_(matrix_.base),
    n_rows_(matrix_.n_rows),
    n_cols_bits_(matrix_.n_cols_bits),
    wpr_(matrix_.wpr)
{   
    if (!matrix_.base) {
        throw std::runtime_error("BitMatrix: null packed_matrix pointer");
    }
    
    // allocate memory for roaring bitmaps
    roaring_bitmap_.resize(n_rows_);

    // allocate memory for col_masses and row_masses
    // NOTE: these are passed from HybridMatrix but are causing a core dump when they are used.
    // This is a temporary fixu until I can work out what is going wrong there.
    // row_masses_.assign(n_rows_, 0);
    // col_masses_.assign(n_cols_bits_, 0);
    ardal::utils::log_info("Constructing Roaring bit-map.");
    // populate roaring bitmaps
    for (size_t i = 0; i < n_rows_; ++i) {
        roaring::Roaring& bitmap = roaring_bitmap_[i];
        for (size_t j = 0; j < wpr_; ++j) {
            uint64_t word = *(base_ + (wpr_ * i) + j);
            for (size_t b = 0; b < 64; ++b) {
                size_t col = j * 64 + b;
                if (col >= n_cols_bits_) break;
                if ((word >> b) & 1) {
                    bitmap.add(col);
                    // row_masses_[i]++;
                    // col_masses_[j]++;
                }
            }
        }
    }

    if (missing_mask && !missing_mask->empty()) {
        if (missing_mask->size() != n_rows_) {
            throw std::runtime_error("RoaringMatrix: missing_mask row count mismatch");
        }
        missing_masks_.resize(n_rows_);
        missing_mask_empty_.resize(n_rows_, true);
        for (size_t i = 0; i < n_rows_; ++i) {
            const auto& cols = (*missing_mask)[i];
            if (cols.empty()) {
                missing_mask_empty_[i] = true;
                continue;
            }
            roaring::Roaring mask;
            for (uint32_t col : cols) {
                if (col < n_cols_bits_) {
                    mask.add(col);
                }
            }
            if (!mask.isEmpty()) {
                missing_masks_[i] = std::move(mask);
                has_missing_mask_ = true;
                missing_mask_empty_[i] = false;
            } else {
                missing_mask_empty_[i] = true;
            }
        }
    } else if (missing_ranges && !missing_ranges->empty()) {
        if (missing_ranges->offsets.size() != n_rows_ + 1) {
            throw std::runtime_error("RoaringMatrix: missing range offsets row count mismatch");
        }
        missing_masks_.resize(n_rows_);
        missing_mask_empty_.resize(n_rows_, true);
        for (size_t i = 0; i < n_rows_; ++i) {
            const uint64_t start_idx = missing_ranges->offsets[i];
            const uint64_t end_idx = missing_ranges->offsets[i + 1];
            if (end_idx < start_idx || end_idx > missing_ranges->ranges.size()) {
                throw std::runtime_error("RoaringMatrix: malformed missing range offsets");
            }
            roaring::Roaring mask;
            for (uint64_t k = start_idx; k < end_idx; ++k) {
                uint32_t left = missing_ranges->ranges[k].first;
                uint32_t right = missing_ranges->ranges[k].second;
                if (right <= left || left >= n_cols_bits_) {
                    continue;
                }
                if (right > n_cols_bits_) {
                    right = static_cast<uint32_t>(n_cols_bits_);
                }
                mask.addRange(left, right);
            }
            if (!mask.isEmpty()) {
                missing_masks_[i] = std::move(mask);
                has_missing_mask_ = true;
                missing_mask_empty_[i] = false;
            }
        }
    } else {
        missing_mask_empty_.assign(n_rows_, true);
    }

    std::stringstream ss;
    ss << "n_rows_ = " << n_rows_ << "; " << "n_cols_bits_ = " << n_cols_bits_ << "; " << "roaring_bitmap_.size() = " << roaring_bitmap_.size();
    ardal::utils::log_debug(static_cast<string>(ss.str()));
}


/****************************************************************************************************
 * ardal::RoaringMatrix::hamming
 *
 * Calculates the Hamming distances between all pairs of rows using Roaring bitmaps.
 *
 * This function calculates the pairwise Hamming distances between all rows of the matrix.
 * The Hamming distance between two rows (Roaring bitmaps) is computed as the cardinality of
 * their symmetric difference (A XOR B).
 *
 * INPUT:
 *  None (operates on the private member roaring_bitmap_)
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<int> : A 1D NumPy array representing the condensed distance matrix containing
 *                     the pairwise Hamming distances. The length of the array is n*(n-1)/2,
 *                     where 'n' is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<uint32_t> RoaringMatrix::hamming( int threads, bool mask_missing ) const {
    ardal::utils::log_info("Running RoaringMatrix::hamming with scalar instructions.");

    const size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<uint32_t> dist_matrix(total_pairs);

    {
        py::gil_scoped_release release;   // kill python

        auto dist_ptr = dist_matrix.mutable_data();

        const bool apply_mask = mask_missing && has_missing_mask_;
        #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
        for (size_t i = 0; i < n_rows_; ++i) {
            for (size_t j = i + 1; j < n_rows_; ++j) {
                size_t idx = (i * (2 * n_rows_ - i - 1)) / 2 + (j - i - 1);
                dist_ptr[idx] = apply_mask ? maskedHammingDistance(i, j)
                                           : hammingDistance(i, j);
            }
        }
    }     
    return dist_matrix;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::hamming_subset
 *
 * Calculate the Hamming distances between pairs of rows from a specified subset of rows and columns.
 *
 * This function first creates a sub-matrix based on the provided row and column indices by
 * intersecting each selected row's bitmap with a mask of the selected columns. It then
 * calculates the pairwise Hamming distances between all rows of this sub-matrix.
 *
 * INPUT:
 *  row_indices (const std::vector<size_t>&) : A vector of row indices to include in the sub-matrix.
 *  col_indices (const std::vector<size_t>&) : A vector of column indices to include in the sub-matrix.
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<uint32_t> : A 1D NumPy array representing the condensed distance matrix.
 ****************************************************************************************************/
py::array_t<uint32_t> RoaringMatrix::hamming_subset( const std::vector<size_t>& row_indices,
                                                     const std::vector<size_t>& col_indices,
                                                     int threads,
                                                     bool mask_missing ) const {
    ardal::utils::log_info("Running RoaringMatrix::hamming_subset with scalar instructions.");

    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }

    for (size_t row_idx : row_indices) {
        if (row_idx >= n_rows_) {
            throw std::out_of_range("Row index in row_indices is out of bounds.");
        }
    }

    const size_t submn_rows_ = row_indices.size();
    if (submn_rows_ <= 1) {
        return py::array_t<uint32_t>(static_cast<py::ssize_t>(0));
    }
    const size_t total_pairs = submn_rows_ * (submn_rows_ - 1) / 2;
    if (col_indices.empty()) {
        py::array_t<uint32_t> zeros(static_cast<py::ssize_t>(total_pairs));
        std::fill_n(zeros.mutable_data(), static_cast<py::ssize_t>(total_pairs), uint32_t{0});
        return zeros;
    }

    const bool use_mask = !isFullColumnSelection(col_indices);
    std::vector<roaring::Roaring> masked_rows;
    std::vector<int> masked_masses;
    std::vector<roaring::Roaring> masked_missing;
    roaring::Roaring column_mask;

    if (use_mask) {
        column_mask = buildColumnMask(col_indices);
        if (column_mask.isEmpty()) {
            py::array_t<uint32_t> zeros(static_cast<py::ssize_t>(total_pairs));
            std::fill_n(zeros.mutable_data(), static_cast<py::ssize_t>(total_pairs), uint32_t{0});
            return zeros;
        }
        masked_rows.resize(submn_rows_);
        masked_masses.resize(submn_rows_);
        if (mask_missing && has_missing_mask_) {
            masked_missing.resize(submn_rows_);
        }
        for (size_t i = 0; i < submn_rows_; ++i) {
            masked_rows[i] = roaring_bitmap_.at(row_indices[i]) & column_mask;
            masked_masses[i] = static_cast<int>(masked_rows[i].cardinality());
            if (!masked_missing.empty()) {
                masked_missing[i] = missing_masks_[row_indices[i]] & column_mask;
            }
        }
    }

    py::array_t<uint32_t> dist_matrix(static_cast<py::ssize_t>(total_pairs));

    py::gil_scoped_release release;   // release GIL for full parallel region
    auto dist_ptr = dist_matrix.mutable_data();

    const bool apply_mask = mask_missing && has_missing_mask_;
    #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t i = 0; i < submn_rows_; ++i) {
        const size_t idx_base = (i * (2 * submn_rows_ - i - 1)) / 2;
        for (size_t j = i + 1; j < submn_rows_; ++j) {
            uint32_t dist_val;
            if (use_mask) {
                if (apply_mask) {
                    roaring::Roaring diff = masked_rows[i] ^ masked_rows[j];
                    if (!masked_missing.empty()) {
                        roaring::Roaring mask_union = masked_missing[i];
                        mask_union |= masked_missing[j];
                        diff -= mask_union;
                    }
                    dist_val = static_cast<uint32_t>(diff.cardinality());
                } else {
                    const int mi = masked_masses[i];
                    const int mj = masked_masses[j];
                    const uint32_t inter = static_cast<uint32_t>(masked_rows[i].and_cardinality(masked_rows[j]));
                    dist_val = static_cast<uint32_t>(mi + mj - (static_cast<int>(inter) << 1));
                }
            } else {
                if (apply_mask) {
                    dist_val = maskedHammingDistance(row_indices[i], row_indices[j]);
                } else {
                    dist_val = hammingDistance(row_indices[i], row_indices[j]);
                }
            }
            dist_ptr[idx_base + (j - i - 1)] = dist_val;
        }
    }

    return dist_matrix;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::hammingDistance
 *
 * Private helper to calculate the Hamming distance between two rows.
 *
 * The Hamming distance between two Roaring bitmaps is the cardinality of their symmetric
 * difference (XOR).
 *
 * INPUT:
 *  i (size_t) : The index of the first row.
 *  j (size_t) : The index of the second row.
 *
 * OUTPUT:
 *  uint32_t : The Hamming distance between the two rows.
 ****************************************************************************************************/
inline uint32_t RoaringMatrix::hammingDistance( size_t i, size_t j ) const {
    // return (roaring_bitmap_[i] ^ roaring_bitmap_[j]).cardinality();
    const int mi = (*row_masses_)[i];
    const int mj = (*row_masses_)[j];
    const uint32_t inter = roaring_bitmap_[i].and_cardinality(roaring_bitmap_[j]);
    return static_cast<uint32_t>(mi + mj - (inter << 1));
}

inline uint32_t RoaringMatrix::maskedHammingDistance( size_t i, size_t j ) const {
    roaring::Roaring diff = roaring_bitmap_[i] ^ roaring_bitmap_[j];
    roaring::Roaring mask_union = missing_masks_[i];
    mask_union |= missing_masks_[j];
    diff -= mask_union;
    return static_cast<uint32_t>(diff.cardinality());
}


/****************************************************************************************************
 * ardal::RoaringMatrix::jaccard
 *
 * Calculates the Jaccard index between all pairs of rows using Roaring bitmaps.
 *
 * This function calculates the pairwise Jaccard index between all rows of the matrix.
 * The Jaccard index between two rows (Roaring bitmaps) is computed as the ratio of
 * the cardinality of their intersection (AND) to the cardinality of their union (OR).
 *              JD = |A ∩ B| / |A ∪ B|
 *
 * INPUT:
 *  None (operates on the private member roaring_bitmap_)
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<double> : A 1D NumPy array representing the condensed distance matrix containing
 *                     the pairwise Jaccard distances. The length of the array is n*(n-1)/2,
 *                     where 'n' is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<double> RoaringMatrix::jaccard( int threads ) const {
    ardal::utils::log_info("Running RoaringMatrix::jaccard with scalar instructions.");
    const size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<double> dist_matrix(total_pairs);
    {
        py::gil_scoped_release release;

        auto dist_ptr = dist_matrix.mutable_data();

        #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
        for (size_t i = 0; i < n_rows_; ++i) {
            for (size_t j = i + 1; j < n_rows_; ++j) {
                size_t idx = (i * (2 * n_rows_ - i - 1)) / 2 + (j - i - 1);
                dist_ptr[idx] = jaccardIndex(i, j);
            }
        }
    }     
    return dist_matrix;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::jaccardIndex
 *
 * Private helper to calculate the Jaccard index between two rows.
 *
 * The Jaccard index is calculated as |A ∩ B| / |A ∪ B|. The union size is derived from
 * the pre-calculated row masses and the intersection size to avoid recomputing the union.
 *
 * INPUT:
 *  i (size_t) : The index of the first row.
 *  j (size_t) : The index of the second row.
 *
 * OUTPUT:
 *  double : The Jaccard index between the two rows. Returns 0.0 if the union is 0.
 *           Note: This returns the Jaccard *Index* (similarity), not distance.
 ****************************************************************************************************/
double RoaringMatrix::jaccardIndex( size_t i, size_t j ) const {
    const auto& row_masses = *row_masses_;

    // native roaring and_cardinality function
    const uint64_t inter = roaring_bitmap_[i].and_cardinality(roaring_bitmap_[j]);

    // union
    const uint64_t uni = row_masses[i] + row_masses[j] - inter;

    if (uni == 0) return 1.0;

    return static_cast<double>(inter) / static_cast<double>(uni);

    // const double intersection_size = (roaring_bitmap_[i] & roaring_bitmap_[j]).cardinality();
    // // const double union_size = row_masses_[i] + row_masses_[j] - intersection_size;
    // // const double union_size = (roaring_bitmap_[i] | roaring_bitmap_[j]).cardinality();
    // const double union_size = row_masses[i] + row_masses[j] - intersection_size;

    // if (union_size == 0) {
    //     // if union is 0, both sets are empty and thus identical. Distance is 0.
    //     return 0.0;
    // }
    // const double jaccard_index = intersection_size / union_size;
    // return jaccard_index;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::neighbourhood
 *
 * Find the epsilon-neighborhood of a row using Hamming distance.
 *
 * This function identifies the rows in the matrix that are within a specified Hamming distance
 * (epsilon) of a given row. It uses pre-calculated row masses for a quick pre-filtering step.
 *
 * INPUT:
 *  row_idx (size_t) : The index of the query row, representing the centroid of the neighbourhood.
 *  epsilon (int)    : The maximum Hamming distance threshold.
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<int64_t> : A 2D NumPy array of shape (N, 2) where each row contains the
 *                         index and Hamming distance of a neighbor.
 ****************************************************************************************************/
py::array_t<int64_t> RoaringMatrix::neighbourhood( size_t row_idx, int epsilon, int threads ) const {
    ardal::utils::log_info("Running RoaringMatrix::neighbourhood with scalar instructions.");

    // do some data cleanliness
    if (epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_idx < 0) {
        throw std::runtime_error("Row index must be non-negative.");
        }
    if (row_idx >= n_rows_) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }
    if (threads <= 0)
        throw std::runtime_error("Number of threads must be positive.");
        
    const auto& row_masses = *row_masses_;
    int q_mass = row_masses[row_idx];    // access pre-calculated row mass

    const roaring::Roaring& q = roaring_bitmap_[row_idx];

    // build candidate list by mass window 
    const int lo = std::max(0, q_mass - epsilon);
    const int hi = q_mass + epsilon;
    std::vector<uint32_t> candidates;
    candidates.reserve(1024);
    for (uint32_t i = 0; i < n_rows_; ++i) {
        if (i == row_idx) continue;
        const int m = row_masses[i];
        if (m >= lo && m <= hi) candidates.push_back(i);
    }

    std::vector<std::vector<std::pair<uint32_t, int>>> thread_results;

    {
        py::gil_scoped_release release;   // release GIL for full parallel region

        #pragma omp parallel num_threads(threads)
        {
            const int thread_id = omp_get_thread_num();

            #pragma omp single
            thread_results.resize(omp_get_num_threads());

            auto& local = thread_results[thread_id];

            #pragma omp for schedule(guided, 2048)
            for (size_t k = 0; k < candidates.size(); ++k) {
                const uint32_t i = candidates[k];
                const int mi = row_masses[i];

                // early impossibility check via required intersection
                // distance <= epsilon  <=>  |A∩B| >= ceil((q_mass + mi - epsilon)/2)
                const int need_num = q_mass + mi - epsilon;
                if (need_num > 0) {
                    const int need = (need_num + 1) >> 1; // ceil(need_num / 2)
                    if (need > std::min(q_mass, mi)) continue; // cannot meet threshold
                }
                
                // due to the restrictions on using a roaring bitmap, early exit is not trivial :(
                int distance = hammingDistance(row_idx, i);

                if (distance <= epsilon)
                    local.emplace_back(i, distance);
            }
        }   // end of parallel region
    }   // GIL reestablished here

    // count total neighbours to preallocate np array
    size_t total_neighbours = 0;
    for (const auto& vec : thread_results) {
        total_neighbours += vec.size();
    }

    // create result np array of shape (total_neighbours, 2)
    py::array_t<int64_t> ep_n({total_neighbours, (size_t)2});
    auto result_ptr = ep_n.mutable_data();

    // populate array with (idex, dist) pairs
    size_t curr_idx = 0;
    for (const auto& vec : thread_results) {
        for (const auto& [idx, dist] : vec) {
            result_ptr[curr_idx * 2 + 0] = idx;
            result_ptr[curr_idx * 2 + 1] = dist;
            ++curr_idx;
        }
    }       
    return ep_n;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::cosineDistance
 *
 * Compute the condensed cosine distance matrix using roaring bitmaps.
 ****************************************************************************************************/
py::array_t<double> RoaringMatrix::cosineDistance( int threads ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (n_rows_ <= 1) {
        return py::array_t<double>(static_cast<py::ssize_t>(0));
    }

    const auto& row_masses = *row_masses_;
    std::vector<double> norms(n_rows_);
    for (size_t i = 0; i < n_rows_; ++i) {
        const int mass = row_masses[i];
        norms[i] = mass > 0 ? std::sqrt(static_cast<double>(mass)) : 0.0;
    }

    const size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<double> dist_condensed(static_cast<py::ssize_t>(total_pairs));
    double* dist_ptr = dist_condensed.mutable_data();

    {
        py::gil_scoped_release release;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (size_t i = 0; i < n_rows_; ++i) {
            const double norm_i = norms[i];
            const size_t idx_base = (i * (2 * n_rows_ - i - 1)) / 2;
            const auto& bitmap_i = roaring_bitmap_[i];

            for (size_t j = i + 1; j < n_rows_; ++j) {
                const double norm_j = norms[j];
                const double denom = norm_i * norm_j;

                double distance;
                if (denom == 0.0) {
                    distance = (norm_i == 0.0 && norm_j == 0.0) ? 0.0 : 1.0;
                } else {
                    const uint64_t dot = bitmap_i.and_cardinality(roaring_bitmap_[j]);
                    double cosine = static_cast<double>(dot) / denom;
                    if (cosine > 1.0) cosine = 1.0;
                    else if (cosine < -1.0) cosine = -1.0;
                    distance = 1.0 - cosine;
                }

                dist_ptr[idx_base + (j - i - 1)] = distance;
            }
        }
    }

    return dist_condensed;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::cosineDistance_subset
 *
 * Compute the condensed cosine distance matrix on a subset of rows/columns.
 ****************************************************************************************************/
py::array_t<double> RoaringMatrix::cosineDistance_subset( const std::vector<size_t>& row_indices,
                                                          const std::vector<size_t>& col_indices,
                                                          int threads ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }

    const size_t submn_rows_ = row_indices.size();
    if (submn_rows_ <= 1) {
        return py::array_t<double>(static_cast<py::ssize_t>(0));
    }

    const size_t total_pairs = submn_rows_ * (submn_rows_ - 1) / 2;
    py::array_t<double> dist_condensed(static_cast<py::ssize_t>(total_pairs));
    double* dist_ptr = dist_condensed.mutable_data();

    if (col_indices.empty()) {
        std::fill(dist_ptr, dist_ptr + total_pairs, 0.0);
        return dist_condensed;
    }

    for (size_t row_idx : row_indices) {
        if (row_idx >= n_rows_) {
            throw std::out_of_range("Row index in row_indices is out of bounds.");
        }
    }

    const bool use_mask = !isFullColumnSelection(col_indices);
    std::vector<double> norms(submn_rows_, 0.0);
    std::vector<roaring::Roaring> masked_rows;

    if (use_mask) {
        roaring::Roaring column_mask = buildColumnMask(col_indices);
        if (column_mask.isEmpty()) {
            std::fill(dist_ptr, dist_ptr + total_pairs, 0.0);
            return dist_condensed;
        }
        masked_rows.resize(submn_rows_);
        for (size_t i = 0; i < submn_rows_; ++i) {
            masked_rows[i] = roaring_bitmap_.at(row_indices[i]) & column_mask;
            const uint64_t mass = masked_rows[i].cardinality();
            norms[i] = mass > 0 ? std::sqrt(static_cast<double>(mass)) : 0.0;
        }
    } else {
        for (size_t i = 0; i < submn_rows_; ++i) {
            const int mass = (*row_masses_)[row_indices[i]];
            norms[i] = mass > 0 ? std::sqrt(static_cast<double>(mass)) : 0.0;
        }
    }


    {
        py::gil_scoped_release release;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (size_t i = 0; i < submn_rows_; ++i) {
            const double norm_i = norms[i];
            const size_t idx_base = (i * (2 * submn_rows_ - i - 1)) / 2;

            for (size_t j = i + 1; j < submn_rows_; ++j) {
                const double norm_j = norms[j];
                const double denom = norm_i * norm_j;

                double distance;
                if (denom == 0.0) {
                    distance = (norm_i == 0.0 && norm_j == 0.0) ? 0.0 : 1.0;
                } else {
                    uint64_t dot;
                    if (use_mask) {
                        dot = masked_rows[i].and_cardinality(masked_rows[j]);
                    } else {
                        dot = roaring_bitmap_[row_indices[i]].and_cardinality(roaring_bitmap_[row_indices[j]]);
                    }
                    double cosine = static_cast<double>(dot) / denom;
                    if (cosine > 1.0) cosine = 1.0;
                    else if (cosine < -1.0) cosine = -1.0;
                    distance = 1.0 - cosine;
                }

                dist_ptr[idx_base + (j - i - 1)] = distance;
            }
        }
    }

    return dist_condensed;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::innerProduct
 *
 * Calculate the inner product between all pairs of rows using Roaring bitmaps.
 *
 * This function calculates the pairwise inner product (number of shared set bits) between all
 * rows of the matrix. For Roaring bitmaps, this is the cardinality of their intersection (AND).
 *
 * INPUT:
 *  None (operates on the private member roaring_bitmap_)
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<int> : A 1D NumPy array representing the condensed inner product matrix.
 ****************************************************************************************************/
py::array_t<int> RoaringMatrix::innerProduct( int threads ) const {
    ardal::utils::log_info("Running RoaringMatrix::innerProduct with scalar instructions.");

    const size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<int> inner_product_matrix(total_pairs);
    {
        py::gil_scoped_release release;

        auto inner_product_ptr = inner_product_matrix.mutable_data();

        #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
        for (size_t i = 0; i < n_rows_; ++i) {
            for (size_t j = i + 1; j < n_rows_; ++j) {
                size_t idx = (i * (2 * n_rows_ - i - 1)) / 2 + (j - i - 1);
                inner_product_ptr[idx] = innerProductRowwise(i, j);
            }
        }   // end of parallel region
    }   // GIL reestablished
    return inner_product_matrix;
}


bool RoaringMatrix::isFullColumnSelection( const std::vector<size_t>& col_indices ) const {
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


roaring::Roaring RoaringMatrix::buildColumnMask( const std::vector<size_t>& col_indices ) const {
    roaring::Roaring mask;
    for (size_t col_idx : col_indices) {
        if (col_idx >= n_cols_bits_) {
            throw std::out_of_range("Column index in col_indices is out of bounds.");
        }
        mask.add(static_cast<uint32_t>(col_idx));
    }
    return mask;
}


std::vector<uint8_t> RoaringMatrix::buildLocusMask( const std::vector<size_t>& col_indices ) const {
    const uint32_t invalid_locus = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> loci;
    loci.reserve(col_indices.size());
    uint32_t max_locus = 0;
    bool has_locus = false;

    for (size_t col_idx : col_indices) {
        if (col_idx >= n_cols_bits_) {
            throw std::out_of_range("Column index in col_indices is out of bounds.");
        }
        const uint32_t locus = allele_to_locus_.empty() ? invalid_locus : allele_to_locus_[col_idx];
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

/****************************************************************************************************
 * ardal::RoaringMatrix::innerProductNeighbourhood
 *
 * Find all rows which share at least `ip_epsilon` set bits with a given row.
 *
 * This function identifies rows in the matrix that have a specified minimum inner product
 * (number of shared bits) with a target query row.
 *
 * INPUT:
 *  row_idx (size_t)    : The index of the query row.
 *  ip_epsilon (int)    : The minimum inner product threshold.
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::list : A Python list of (row_index, inner_product) tuples for all neighbors found.
 *
 ****************************************************************************************************/
py::list RoaringMatrix::innerProductNeighbourhood( size_t row_idx, int ip_epsilon, int threads ) const {
    ardal::utils::log_info("Running RoaringMatrix::innerProductNeighbourhood with scalar instructions.");

    // do some data cleanliness
    if (ip_epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_idx < 0) {
        throw std::runtime_error("Row index must be non-negative.");
        }
    if (row_idx >= n_rows_) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }
    if (threads <= 0)
        throw std::runtime_error("Number of threads must be positive.");
    
    std::vector<std::vector<std::pair<size_t, int>>> thread_results;

    {
        py::gil_scoped_release release;   // release GIL for full parallel region

        #pragma omp parallel num_threads(threads)
        {
            const int thread_id = omp_get_thread_num();

            #pragma omp single
            thread_results.resize(omp_get_num_threads());

            auto& local = thread_results[thread_id];

            #pragma omp for schedule(static)
            for (size_t i = 0; i < n_rows_; ++i) {
                if (i == row_idx) continue;
                
                int distance = innerProductRowwise(row_idx, i);

                if (distance >= ip_epsilon)
                    local.emplace_back(i, distance);
            }
        }   // end of parallel region
    }   // GIL reestablished here

    py::list ipe_n;
    for (const auto& vec : thread_results) {
        for (const auto& [idx, dist] : vec) {
            ipe_n.append(py::make_tuple(idx, dist));
        }
    }
    return ipe_n;
}


py::list RoaringMatrix::knn_hamming( size_t row_idx, int k, int threads ) const {
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

    const size_t n = n_rows_;
    if (row_idx >= n || n == 0 || k <= 0) {
        return py::list();
    }
    const uint32_t k_eff = std::min<uint32_t>(static_cast<uint32_t>(k), n > 1 ? uint32_t(n - 1) : 0);
    if (k_eff == 0) {
        return py::list();
    }

    ardal::utils::log_info("Running RoaringMatrix::knn_hamming.");

    std::vector<Neighbor> final_neighbors;
    {
        py::gil_scoped_release release;

        const int T = std::max(1, threads);
        std::vector<std::vector<Neighbor>> buckets(T);

#if defined(_OPENMP)
        #pragma omp parallel num_threads(T)
#endif
        {
            std::priority_queue<Neighbor, std::vector<Neighbor>, ByMaxDistance> heap;

#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for (uint32_t j = 0; j < n; ++j) {
                if (j == row_idx) continue;

                const uint32_t dist = hammingDistance(row_idx, j);

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

            std::vector<Neighbor> local;
            local.reserve(heap.size());
            while (!heap.empty()) {
                local.push_back(heap.top());
                heap.pop();
            }
            std::sort(local.begin(), local.end(), asc_dist_id);

#if defined(_OPENMP)
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            buckets[tid] = std::move(local);
        }

        size_t total = 0;
        for (const auto& v : buckets) total += v.size();
        std::vector<Neighbor> all;
        all.reserve(total);
        for (auto& v : buckets) {
            all.insert(all.end(), v.begin(), v.end());
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
        result.append(py::make_tuple(static_cast<int64_t>(nb.id), static_cast<int64_t>(nb.dist)));
    }
    return result;
}


py::list RoaringMatrix::knn_inner_product( size_t row_idx, int k, int threads ) const {
    struct Candidate {
        uint32_t id;
        uint32_t score;
    };
    struct MinScore {
        bool operator()(const Candidate& a, const Candidate& b) const noexcept {
            if (a.score == b.score) return a.id > b.id;
            return a.score > b.score;
        }
    };
    struct DescScore {
        bool operator()(const Candidate& a, const Candidate& b) const noexcept {
            if (a.score == b.score) return a.id < b.id;
            return a.score > b.score;
        }
    };

    const size_t n = n_rows_;
    if (row_idx >= n || n == 0 || k <= 0) {
        return py::list();
    }
    const uint32_t k_eff = std::min<uint32_t>(static_cast<uint32_t>(k), n > 1 ? uint32_t(n - 1) : 0);
    if (k_eff == 0) {
        return py::list();
    }

    ardal::utils::log_info("Running RoaringMatrix::knn_inner_product.");

    std::vector<Candidate> final_candidates;
    {
        py::gil_scoped_release release;

        const int T = std::max(1, threads);
        std::vector<std::vector<Candidate>> buckets(T);

#if defined(_OPENMP)
        #pragma omp parallel num_threads(T)
#endif
        {
            std::priority_queue<Candidate, std::vector<Candidate>, MinScore> heap;

#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for (uint32_t j = 0; j < n; ++j) {
                if (j == row_idx) continue;

                const uint32_t score = innerProductRowwise(row_idx, j);

                if (heap.size() < k_eff) {
                    heap.push(Candidate{j, score});
                } else {
                    const Candidate& top = heap.top();
                    if (score > top.score || (score == top.score && j < top.id)) {
                        heap.pop();
                        heap.push(Candidate{j, score});
                    }
                }
            }

            std::vector<Candidate> local;
            local.reserve(heap.size());
            while (!heap.empty()) {
                local.push_back(heap.top());
                heap.pop();
            }
            std::sort(local.begin(), local.end(), DescScore{});

#if defined(_OPENMP)
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            buckets[tid] = std::move(local);
        }

        size_t total = 0;
        for (const auto& v : buckets) total += v.size();
        std::vector<Candidate> all;
        all.reserve(total);
        for (auto& v : buckets) {
            all.insert(all.end(), v.begin(), v.end());
        }

        DescScore desc;
        if (all.size() > k_eff) {
            std::nth_element(all.begin(), all.begin() + k_eff, all.end(), desc);
            all.resize(k_eff);
        }
        std::sort(all.begin(), all.end(), desc);
        final_candidates = std::move(all);
    }

    py::list result;
    for (const auto& cand : final_candidates) {
        result.append(py::make_tuple(static_cast<int64_t>(cand.id), static_cast<int64_t>(cand.score)));
    }
    return result;
}


py::list RoaringMatrix::knn_jaccard( size_t row_idx, int k, int threads ) const {
    struct Candidate {
        uint32_t id;
        double score;
    };
    struct MinScore {
        bool operator()(const Candidate& a, const Candidate& b) const noexcept {
            if (a.score == b.score) return a.id > b.id;
            return a.score > b.score;
        }
    };
    struct DescScore {
        bool operator()(const Candidate& a, const Candidate& b) const noexcept {
            if (a.score == b.score) return a.id < b.id;
            return a.score > b.score;
        }
    };

    const size_t n = n_rows_;
    if (row_idx >= n || n == 0 || k <= 0) {
        return py::list();
    }
    const uint32_t k_eff = std::min<uint32_t>(static_cast<uint32_t>(k), n > 1 ? uint32_t(n - 1) : 0);
    if (k_eff == 0) {
        return py::list();
    }

    ardal::utils::log_info("Running RoaringMatrix::knn_jaccard.");

    const uint32_t mass_query = static_cast<uint32_t>((*row_masses_)[row_idx]);

    std::vector<Candidate> final_candidates;
    {
        py::gil_scoped_release release;

        const int T = std::max(1, threads);
        std::vector<std::vector<Candidate>> buckets(T);

#if defined(_OPENMP)
        #pragma omp parallel num_threads(T)
#endif
        {
            std::priority_queue<Candidate, std::vector<Candidate>, MinScore> heap;

#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for (uint32_t j = 0; j < n; ++j) {
                if (j == row_idx) continue;

                const uint32_t intersection = innerProductRowwise(row_idx, j);
                const uint32_t mass_other = static_cast<uint32_t>((*row_masses_)[j]);
                const uint32_t union_cnt = mass_query + mass_other - intersection;

                double score;
                if (union_cnt == 0) {
                    score = (mass_query == 0 && mass_other == 0) ? 1.0 : 0.0;
                } else {
                    score = static_cast<double>(intersection) / static_cast<double>(union_cnt);
                }

                if (heap.size() < k_eff) {
                    heap.push(Candidate{j, score});
                } else {
                    const Candidate& top = heap.top();
                    if (score > top.score || (score == top.score && j < top.id)) {
                        heap.pop();
                        heap.push(Candidate{j, score});
                    }
                }
            }

            std::vector<Candidate> local;
            local.reserve(heap.size());
            while (!heap.empty()) {
                local.push_back(heap.top());
                heap.pop();
            }
            std::sort(local.begin(), local.end(), DescScore{});

#if defined(_OPENMP)
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            buckets[tid] = std::move(local);
        }

        size_t total = 0;
        for (const auto& v : buckets) total += v.size();
        std::vector<Candidate> all;
        all.reserve(total);
        for (auto& v : buckets) {
            all.insert(all.end(), v.begin(), v.end());
        }

        DescScore desc;
        if (all.size() > k_eff) {
            std::nth_element(all.begin(), all.begin() + k_eff, all.end(), desc);
            all.resize(k_eff);
        }
        std::sort(all.begin(), all.end(), desc);
        final_candidates = std::move(all);
    }

    py::list result;
    for (const auto& cand : final_candidates) {
        result.append(py::make_tuple(static_cast<int64_t>(cand.id), cand.score));
    }
    return result;
}


py::list RoaringMatrix::knn_cosine( size_t row_idx, int k, int threads ) const {
    struct Candidate {
        uint32_t id;
        double score;
    };
    struct MinScore {
        bool operator()(const Candidate& a, const Candidate& b) const noexcept {
            if (a.score == b.score) return a.id > b.id;
            return a.score > b.score;
        }
    };
    struct DescScore {
        bool operator()(const Candidate& a, const Candidate& b) const noexcept {
            if (a.score == b.score) return a.id < b.id;
            return a.score > b.score;
        }
    };

    const size_t n = n_rows_;
    if (row_idx >= n || n == 0 || k <= 0) {
        return py::list();
    }
    const uint32_t k_eff = std::min<uint32_t>(static_cast<uint32_t>(k), n > 1 ? uint32_t(n - 1) : 0);
    if (k_eff == 0) {
        return py::list();
    }

    ardal::utils::log_info("Running RoaringMatrix::knn_cosine.");

    const uint32_t mass_query = static_cast<uint32_t>((*row_masses_)[row_idx]);
    const double norm_query = mass_query > 0 ? std::sqrt(static_cast<double>(mass_query)) : 0.0;

    std::vector<Candidate> final_candidates;
    {
        py::gil_scoped_release release;

        const int T = std::max(1, threads);
        std::vector<std::vector<Candidate>> buckets(T);

#if defined(_OPENMP)
        #pragma omp parallel num_threads(T)
#endif
        {
            std::priority_queue<Candidate, std::vector<Candidate>, MinScore> heap;

#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for (uint32_t j = 0; j < n; ++j) {
                if (j == row_idx) continue;

                const uint32_t intersection = innerProductRowwise(row_idx, j);
                const uint32_t mass_other = static_cast<uint32_t>((*row_masses_)[j]);
                const double norm_other = mass_other > 0 ? std::sqrt(static_cast<double>(mass_other)) : 0.0;
                const double denom = norm_query * norm_other;

                double score = 0.0;
                if (denom > 0.0) {
                    score = static_cast<double>(intersection) / denom;
                    if (score > 1.0) score = 1.0;
                } else if (norm_query == 0.0 && norm_other == 0.0) {
                    score = 0.0;
                }

                if (heap.size() < k_eff) {
                    heap.push(Candidate{j, score});
                } else {
                    const Candidate& top = heap.top();
                    if (score > top.score || (score == top.score && j < top.id)) {
                        heap.pop();
                        heap.push(Candidate{j, score});
                    }
                }
            }

            std::vector<Candidate> local;
            local.reserve(heap.size());
            while (!heap.empty()) {
                local.push_back(heap.top());
                heap.pop();
            }
            std::sort(local.begin(), local.end(), DescScore{});

#if defined(_OPENMP)
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            buckets[tid] = std::move(local);
        }

        size_t total = 0;
        for (const auto& v : buckets) total += v.size();
        std::vector<Candidate> all;
        all.reserve(total);
        for (auto& v : buckets) {
            all.insert(all.end(), v.begin(), v.end());
        }

        DescScore desc;
        if (all.size() > k_eff) {
            std::nth_element(all.begin(), all.begin() + k_eff, all.end(), desc);
            all.resize(k_eff);
        }
        std::sort(all.begin(), all.end(), desc);
        final_candidates = std::move(all);
    }

    py::list result;
    for (const auto& cand : final_candidates) {
        result.append(py::make_tuple(static_cast<int64_t>(cand.id), cand.score));
    }
    return result;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::innerProductRowwise
 *
 * Private helper to calculate the inner product between two rows.
 *
 * The inner product between two Roaring bitmaps is the cardinality of their intersection (AND).
 *
 * INPUT:
 *  i (size_t) : The index of the first row.
 *  j (size_t) : The index of the second row.
 *
 * OUTPUT:
 *  uint32_t : The inner product of the two rows.
 ****************************************************************************************************/
uint32_t RoaringMatrix::innerProductRowwise( size_t i, size_t j ) const {
    return (roaring_bitmap_[i] & roaring_bitmap_[j]).cardinality();
}


/****************************************************************************************************
 * ardal::RoaringMatrix::uniqueSharedBits
 *
 * Finds the indices of bits that are set (1) in ALL specified "ingroup" rows and are NOT set
 * in ANY "outgroup" row (all other rows).
 *
 * This function computes the intersection of all ingroup bitmaps, the union of all outgroup
 * bitmaps, and then finds the difference between these two results.
 *
 * INPUT:
 *  row_indices (const std::vector<size_t>&) : A vector of row indices for the ingroup.
 *
 * OUTPUT:
 *  std::vector<size_t> : A vector containing the column indices of the unique shared bits.
 ****************************************************************************************************/
std::vector<size_t> RoaringMatrix::uniqueSharedBits(const std::vector<size_t>& row_indices) const {
    ardal::utils::log_info("Running RoaringMatrix::uniqueSharedBits with scalar instructions.");

    if (row_indices.empty()) return {};

    // --- Step 1: Intersect ingroup bitmaps ---
    roaring::Roaring ingroup_shared = roaring_bitmap_.at(row_indices[0]);
    for (size_t i = 1; i < row_indices.size(); ++i) {
        ingroup_shared &= roaring_bitmap_.at(row_indices[i]);
    }

    // --- Step 2: Union outgroup bitmaps ---
    roaring::Roaring outgroup_union;
    for (size_t i = 0; i < n_rows_; ++i) {
        if (std::find(row_indices.begin(), row_indices.end(), i) != row_indices.end()) continue;
        outgroup_union |= roaring_bitmap_.at(i);
    }

    // --- Step 3: Compute unique shared bits ---
    roaring::Roaring unique_bits = ingroup_shared - outgroup_union;

    // --- Step 4: Extract to std::vector<size_t> ---
    size_t cardinality = unique_bits.cardinality();
    std::vector<uint32_t> u32_vals(cardinality);
    unique_bits.toUint32Array(u32_vals.data());

    // Cast to size_t for safety (Python side expects size_t)
    return std::vector<size_t>(u32_vals.begin(), u32_vals.end());
}


/****************************************************************************************************
 * ardal::RoaringMatrix::getSetBitIndices
 *
 * Get the indices of set bits for a given row.
 *
 * This function retrieves the column indices of set bits that are present in a specified row of the
 * matrix and returns them as a NumPy array.
 *
 * INPUT:
 *  row_idx (size_t) : The index of the row in the matrix.
 *
 * OUTPUT:
 *  py::array_t<size_t> : A 1D NumPy array containing the column indices for the specified row.
 ****************************************************************************************************/
py::array_t<size_t> RoaringMatrix::getSetBitIndices( size_t row_idx ) const {
    const auto& bitmap = roaring_bitmap_.at(row_idx);

    // Create a numpy array of uint32_t to hold the values.
    py::array_t<uint32_t> values_u32(bitmap.cardinality());

    // Fill the numpy array with values from the roaring bitmap.
    bitmap.toUint32Array(values_u32.mutable_data());

    // Cast the numpy array to size_t, which is what the python side expects.
    return values_u32.cast<py::array_t<size_t>>();
}


/****************************************************************************************************
 * ardal::RoaringMatrix::getRoaringMatrix
 *
 * Get the entire Roaring matrix as a Python list of NumPy arrays.
 *
 * This function iterates through each row of the Roaring matrix and converts it into a NumPy
 * array of set bit indices, returning the result as a list of these arrays.
 *
 * OUTPUT:
 *  py::list : A list where each element is a NumPy array representing the set bits of a row.
 *
 ****************************************************************************************************/
py::list RoaringMatrix::getRoaringMatrix( void ) const {
    if (roaring_bitmap_.size() != n_rows_) {
        std::stringstream ss;
        ss << "roaring_bitmap_ size (" << roaring_bitmap_.size() << ") does not match n_rows_ (" << n_rows_ <<")";
        throw std::runtime_error(static_cast<string>(ss.str()));
    }
    py::list roaring_matrix;
    for (size_t i = 0; i < n_rows_; ++i) {
        // const auto& bitmap = roaring_bitmap_[i];
        roaring_matrix.append(getSetBitIndices(i));
    }
    return roaring_matrix;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::colwiseRoaringFromRowwise
 *
 * Transposes the row-wise Roaring matrix into a column-wise representation.
 *
 * This private helper function creates a new set of Roaring bitmaps where each bitmap
 * represents a column and contains the indices of the rows that have a bit set in that column.
 * This is useful for column-centric operations like `bitCooccurrence`.
 *
 * OUTPUT:
 *  std::vector<roaring::Roaring> : A vector of Roaring bitmaps representing the columns.
 ****************************************************************************************************/
std::vector<roaring::Roaring> RoaringMatrix::colwiseRoaringFromRowwise( void ) const {
    ardal::utils::log_info("Constructing col-wise Roaring bitmap from row-wise.");
    // create one roaring bitmap per column
    std::vector<roaring::Roaring> colwise_roaring(n_cols_bits_);

    // for each row, iterate over its set bits and add the row index to the corresponding column bitmap
    for (size_t row = 0; row < n_rows_; ++row) {
        const auto& row_bitmap = roaring_bitmap_[row];
        // use an iterator to efficiently access set bits
        for (auto it = row_bitmap.begin(); it != row_bitmap.end(); ++it) {
            const size_t col = *it;
            if (col < n_cols_bits_) {
                colwise_roaring[col].add(row);
            }
        }
    }
    return colwise_roaring;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::bitCooccurrence_all
 *
 * Calculates the co-occurrence of bits across all columns in the matrix.
 *
 * This function computes the Jaccard index for every pair of columns. If the index is above
 * a given threshold, the pair is considered co-occurring. The function first transposes the
 * matrix to a column-wise representation for efficiency.
 *
 * INPUT:
 *  threshold (double) : The Jaccard index threshold for co-occurrence (0.0 to 1.0).
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::dict : A dictionary where keys are column indices and values are lists of co-occurring partners.
 ****************************************************************************************************/
py::dict RoaringMatrix::bitCooccurrence_all( double threshold, int threads ) const {
    ardal::utils::log_info("Running RoaringMatrix::bitCooccurrence_all with scalar instructions.");

    // do some input cleanliness
    if (threshold < 0 || threshold > 1) {
        throw std::runtime_error("threshold must be between 0 and 1.");
    }
    if (threads <= 0) {
        throw std::runtime_error("Number of threads must be positive.");
    }
    if (n_rows_ == 0) {
        return py::dict();
    }

    // get the colwise roaring bitmap
    auto colwise_roaring = colwiseRoaringFromRowwise();
    size_t n_cols = colwise_roaring.size();
    std::vector<uint64_t> cardinalities(n_cols, 0);
    for (size_t col = 0; col < n_cols; ++col) {
        cardinalities[col] = colwise_roaring[col].cardinality();
    }

    // sanity checks
    if (roaring_bitmap_.size() != n_rows_) {
        throw std::runtime_error("row-wise roaring size mismatch");
    }
    if (colwise_roaring.size() != n_cols_bits_) {
        throw std::runtime_error("col-wise roaring size mismatch");
    }

    // underflow safety
    if (n_cols < 2) return py::dict();


    // uncomment to debug
    std::stringstream ss;
    ss << "n_rows_ = " << n_rows_ << "; " << "n_cols_bits_ = " << n_cols_bits_ << "; " << "roaring_bitmap_.size() = " << roaring_bitmap_.size() << "; " << "colwise_roaring.size() = " << colwise_roaring.size();
    ardal::utils::log_debug(static_cast<string>(ss.str()));


    std::map<size_t, std::vector<size_t>> global_map;

    {
        py::gil_scoped_release release;   // murder GIL

        #pragma omp parallel num_threads(threads)
        {
            std::map<size_t, std::vector<size_t>> local_map;

            #pragma omp for schedule(static)
            for (size_t i = 0; i < n_cols - 1; ++i) {
                const auto& i_bitmap = colwise_roaring[i];
                const uint64_t i_cardinality = cardinalities[i];
                std::vector<size_t> cooccurring_partners;
                for (size_t j = i + 1; j < n_cols; ++j) {
                    const auto& j_bitmap = colwise_roaring[j];
                    const uint64_t j_cardinality = cardinalities[j];
                    const uint64_t max_cardinality = std::max(i_cardinality, j_cardinality);
                    if (max_cardinality == 0) continue;
                    const double max_possible_jaccard =
                        static_cast<double>(std::min(i_cardinality, j_cardinality)) /
                        static_cast<double>(max_cardinality);
                    if (max_possible_jaccard < threshold) continue;

                    double intersection_size = (i_bitmap & j_bitmap).cardinality();
                    double union_size = static_cast<double>(i_cardinality + j_cardinality) - intersection_size;

                    if (union_size == 0) continue;
                    double jaccard_index = intersection_size / union_size;
                    if (jaccard_index >= threshold) {
                        cooccurring_partners.push_back(j);
                    }
                }
                // std::stringstream ss;
                // ss << "col " << i;
                // ardal::utils::log_info(static_cast<string>(ss.str()));

                // local_map[i] = std::move(cooccurring_partners);
                if (!cooccurring_partners.empty()) {
                    local_map[i] = std::move(cooccurring_partners);
                }
            }

            #pragma omp critical
            {
                for (const auto& [k, v] : local_map) {
                    global_map[k].insert(global_map[k].end(), v.begin(), v.end());
                }
            }
        } // end parallel region
    } // re-establish GIL

    // convert to py dict
    py::dict cooccurrences_py;
    for (const auto& [k, v] : global_map) {
        cooccurrences_py[py::int_(k)] = py::cast(v);
    }
    return cooccurrences_py;
}



/****************************************************************************************************
 * ardal::RoaringMatrix::bitCooccurrence_subset
 *
 * Calculates the co-occurrence of bits for a specified subset of columns.
 *
 * This function computes the Jaccard index for pairs of columns within the provided subset.
 * If the index is above a given threshold, the pair is considered co-occurring.
 *
 * INPUT:
 *  col_indices (const std::vector<size_t>&) : A vector of column indices to analyze.
 *  threshold (double)                       : The Jaccard index threshold for co-occurrence.
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::dict : A dictionary where keys are column indices and values are lists of co-occurring
 *             partners from the subset.
 ****************************************************************************************************/
py::dict RoaringMatrix::bitCooccurrence_subset( const std::vector<size_t>& col_indices, double threshold, int threads ) const {
    ardal::utils::log_info("Running RoaringMatrix::bitCooccurrence_subset with scalar instructions.");

    // do some input cleanliness
    if (threshold < 0 || threshold > 1) {
        throw std::runtime_error("threshold must be between 0 and 1.");
    }
    if (threads <= 0) {
        throw std::runtime_error("Number of threads must be positive.");
    }
    if (n_rows_ == 0) {
        return py::dict();
    }
    if (col_indices.empty() || col_indices.size() == 1) {
        return py::dict();
    }

    std::vector<size_t> subset_cols = col_indices;
    std::sort(subset_cols.begin(), subset_cols.end());
    subset_cols.erase(std::unique(subset_cols.begin(), subset_cols.end()), subset_cols.end());

    // check if col_indices are valid
    for (const auto& idx : subset_cols) {
        if (idx >= n_cols_bits_) {
            throw std::runtime_error("Column index out of bounds.");
        }
    }
    if (subset_cols.size() < 2) {
        return py::dict();
    }

    // get the colwise roaring bitmap
    auto colwise_roaring = colwiseRoaringFromRowwise();
    std::vector<uint64_t> cardinalities(subset_cols.size(), 0);
    for (size_t subset_idx = 0; subset_idx < subset_cols.size(); ++subset_idx) {
        cardinalities[subset_idx] = colwise_roaring[subset_cols[subset_idx]].cardinality();
    }

    std::map<size_t, std::vector<size_t>> global_map;

    {
        py::gil_scoped_release release;   // murder GIL

        #pragma omp parallel num_threads(threads)
        {
            std::map<size_t, std::vector<size_t>> local_map;

            #pragma omp for schedule(static)
            for (size_t subset_i = 0; subset_i < subset_cols.size() - 1; ++subset_i) {
                const size_t i = subset_cols[subset_i];
                std::vector<size_t> cooccurring_partners;
                const auto& i_bitmap = colwise_roaring[i];
                const uint64_t i_cardinality = cardinalities[subset_i];
                for (size_t subset_j = subset_i + 1; subset_j < subset_cols.size(); ++subset_j) {
                    const size_t j = subset_cols[subset_j];
                    const auto& j_bitmap = colwise_roaring[j];
                    const uint64_t j_cardinality = cardinalities[subset_j];
                    const uint64_t max_cardinality = std::max(i_cardinality, j_cardinality);
                    if (max_cardinality == 0) continue;
                    const double max_possible_jaccard =
                        static_cast<double>(std::min(i_cardinality, j_cardinality)) /
                        static_cast<double>(max_cardinality);
                    if (max_possible_jaccard < threshold) continue;

                    double intersection_size = (i_bitmap & j_bitmap).cardinality();
                    double union_size = static_cast<double>(i_cardinality + j_cardinality) - intersection_size;
                    if (union_size == 0) continue;
                    double jaccard_index = intersection_size / union_size;
                    if (jaccard_index >= threshold) {
                        cooccurring_partners.push_back(j);
                    }
                }

                // local_map[i] = std::move(cooccurring_partners);
                if (!cooccurring_partners.empty()) {
                    local_map[i] = std::move(cooccurring_partners);
                }
            }

            #pragma omp critical
            {
                for (const auto& [k, v] : local_map) {
                    global_map[k].insert(global_map[k].end(), v.begin(), v.end());
                }
            }
        } // end parallel region
    } // re-establish GIL

    // convert to py dict
    py::dict cooccurrences_py;
    for (const auto& [k, v] : global_map) {
        cooccurrences_py[py::int_(k)] = py::cast(v);
    }
    return cooccurrences_py;
}


void RoaringMatrix::prepareSnvView( py::array_t<uint32_t> allele_to_locus,
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
    snv_entries_.clear();
    snv_lookup_loaded_ = true;
    snv_vectors_ready_ = false;
    snv_entries_ready_ = false;

    ardal::utils::log_info("RoaringMatrix SNV lookup loaded; vectors will be materialised on demand.");
}


void RoaringMatrix::ensure_snv_vectors() const {
    if (snv_vectors_ready_ && (!has_missing_mask_ || snv_entries_ready_)) return;
    if (!snv_lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    snv_vectors_.assign(n_rows_, {});
    if (has_missing_mask_) {
        snv_entries_.assign(n_rows_, {});
        snv_entries_ready_ = false;
    }
    const uint32_t invalid_locus = std::numeric_limits<uint32_t>::max();

    for (std::size_t row = 0; row < n_rows_; ++row) {
        const auto& bitmap = roaring_bitmap_[row];
        if (bitmap.isEmpty()) continue;

        const roaring::Roaring* self_missing = (has_missing_mask_ && row < missing_masks_.size() && !missing_masks_[row].isEmpty())
                                               ? &missing_masks_[row]
                                               : nullptr;

        std::vector<std::tuple<uint32_t, uint8_t, uint32_t>> collected;
        collected.reserve(bitmap.cardinality());

        for (auto it = bitmap.begin(); it != bitmap.end(); ++it) {
            const uint32_t col = *it;
            if (col >= n_cols_bits_) continue;
            if (self_missing && self_missing->contains(col)) continue;

            const uint32_t locus_id = allele_to_locus_[col];
            if (locus_id == invalid_locus) continue;

            const uint8_t base_code = allele_to_base_[col];
            if (base_code == 0) continue;

            collected.emplace_back(locus_id, base_code, col);
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

    snv_vectors_ready_ = true;
    if (has_missing_mask_) {
        snv_entries_ready_ = true;
    }
}


/****************************************************************************************************
 * ardal::RoaringMatrix::snvDistanceRaw
 *
 * Compute the SNV Hamming distance between two encoded vectors, optionally filtering by locus mask.
 *
 * INPUT:
 *  lhs (const std::vector<uint64_t>&) : Encoded SNV vector for the left row (locus << 3 | base).
 *  rhs (const std::vector<uint64_t>&) : Encoded SNV vector for the right row.
 *  locus_mask (const std::vector<uint8_t>*): Optional locus allow mask; nullptr to disable.
 *
 * OUTPUT:
 *  uint32_t : The SNV Hamming distance between lhs and rhs.
 ****************************************************************************************************/
uint32_t RoaringMatrix::snvDistanceRaw( const std::vector<uint64_t>& lhs,
                                        const std::vector<uint64_t>& rhs,
                                        const std::vector<uint8_t>* locus_mask ) const {
    const bool use_mask = (locus_mask && !locus_mask->empty());
    if (!use_mask) {
        std::size_t p = 0, q = 0;
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
                ++p; ++q;
            } else if (locus_l < locus_r) {
                ++dist; ++p;
            } else {
                ++dist; ++q;
            }
        }
        dist += static_cast<uint32_t>((lhs.size() - p) + (rhs.size() - q));
        return dist;
    }

    constexpr uint32_t INVALID = std::numeric_limits<uint32_t>::max();
    auto advance = [&](const std::vector<uint64_t>& vec, std::size_t& pos) -> uint32_t {
        while (pos < vec.size()) {
            const uint32_t locus = static_cast<uint32_t>(vec[pos] >> 3);
            if (locus < locus_mask->size() && (*locus_mask)[locus]) {
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
            if (locus < locus_mask->size() && (*locus_mask)[locus]) {
                ++count;
            }
            ++pos;
        }
        return count;
    };

    std::size_t p = 0, q = 0;
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
            ++p; ++q;
        } else if (locus_l < locus_r) {
            ++dist; ++p;
        } else {
            ++dist; ++q;
        }
    }
    return dist;
}


/****************************************************************************************************
 * ardal::RoaringMatrix::snvDistanceMaskedEntries
 *
 * Compute SNV distance between two per-row entry lists, applying per-row missing masks and an
 * optional locus mask.
 *
 * INPUT:
 *  lhs (const std::vector<SnvEntry>&) : SNV entries for the left row (col, locus, base).
 *  rhs (const std::vector<SnvEntry>&) : SNV entries for the right row.
 *  mask_l (const roaring::Roaring*)   : Missing mask for the left row (columns to skip).
 *  mask_r (const roaring::Roaring*)   : Missing mask for the right row.
 *  locus_mask (const std::vector<uint8_t>*): Optional locus allow mask; nullptr to disable.
 *
 * OUTPUT:
 *  uint32_t : The SNV Hamming distance with masks applied.
 ****************************************************************************************************/
uint32_t RoaringMatrix::snvDistanceMaskedEntries( const std::vector<SnvEntry>& lhs,
                                                  const std::vector<SnvEntry>& rhs,
                                                  const roaring::Roaring* mask_l,
                                                  const roaring::Roaring* mask_r,
                                                  const std::vector<uint8_t>* locus_mask ) const {
    const bool use_locus_mask = (locus_mask && !locus_mask->empty());
    auto advance_valid = [&](const std::vector<SnvEntry>& vec,
                             std::size_t& idx,
                             const roaring::Roaring* self_mask,
                             const roaring::Roaring* other_mask) -> bool {
        while (idx < vec.size()) {
            const auto& e = vec[idx];
            if (use_locus_mask) {
                if (e.locus >= locus_mask->size() || (*locus_mask)[e.locus] == 0) {
                    ++idx;
                    continue;
                }
            }
            if ((self_mask && self_mask->contains(e.col)) ||
                (other_mask && other_mask->contains(e.col))) {
                ++idx;
                continue;
            }
            return true;
        }
        return false;
    };

    std::size_t p = 0, q = 0;
    uint32_t dist = 0;
    bool has_l = advance_valid(lhs, p, mask_l, mask_r);
    bool has_r = advance_valid(rhs, q, mask_r, mask_l);

    while (has_l || has_r) {
        if (!has_l) {
            ++dist;
            has_r = advance_valid(rhs, q, mask_r, mask_l);
            continue;
        }
        if (!has_r) {
            ++dist;
            has_l = advance_valid(lhs, p, mask_l, mask_r);
            continue;
        }

        const auto& le = lhs[p];
        const auto& re = rhs[q];
        if (le.locus == re.locus) {
            if (le.base != re.base) ++dist;
            has_l = advance_valid(lhs, ++p, mask_l, mask_r);
            has_r = advance_valid(rhs, ++q, mask_r, mask_l);
        } else if (le.locus < re.locus) {
            ++dist;
            has_l = advance_valid(lhs, ++p, mask_l, mask_r);
        } else {
            ++dist;
            has_r = advance_valid(rhs, ++q, mask_r, mask_l);
        }
    }
    return dist;
}


uint32_t RoaringMatrix::snvDistance( size_t i, size_t j ) const {
    const auto& lhs = snv_vectors_[i];
    const auto& rhs = snv_vectors_[j];
    return snvDistanceRaw(lhs, rhs, nullptr);
}


/****************************************************************************************************
 * ardal::RoaringMatrix::snvDistanceMasked
 *
 * Compute SNV distance for two rows, honoring missing masks and an optional locus mask.
 *
 * INPUT:
 *  i (size_t)                         : Row index for the left operand.
 *  j (size_t)                         : Row index for the right operand.
 *  locus_mask (const std::vector<uint8_t>*): Optional locus allow mask; nullptr to disable.
 *
 * OUTPUT:
 *  uint32_t : The SNV Hamming distance with masks applied.
 ****************************************************************************************************/
uint32_t RoaringMatrix::snvDistanceMasked( size_t i,
                                           size_t j,
                                           const std::vector<uint8_t>* locus_mask ) const {
    return snvDistanceRaw(snv_vectors_[i], snv_vectors_[j], (locus_mask && !locus_mask->empty()) ? locus_mask : nullptr);
}


/****************************************************************************************************
 * ardal::RoaringMatrix::snvNeighbourhood
 *
 * Find all rows within a given SNV Hamming radius of the query row.
 *
 * INPUT:
 *  row_idx (size_t) : The index of the query row.
 *  epsilon (uint32_t) : Maximum SNV Hamming distance.
 *
 * PARAMETERS:
 *  threads (int) : Number of threads for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<int64_t> : Array of (row_index, distance) pairs.
 ****************************************************************************************************/
py::array_t<int64_t> RoaringMatrix::snvNeighbourhood( size_t row_idx,
                                                      uint32_t epsilon,
                                                      int threads,
                                                      bool mask_missing ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (row_idx >= n_rows_) {
        throw std::runtime_error("Row index must be non-negative and within bounds.");
    }
    if (!snv_lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    ensure_snv_vectors();
    (void)mask_missing;
    ardal::utils::log_info("Running RoaringMatrix::snvNeighbourhood.");

    const int T = std::max(1, threads);
    std::vector<std::vector<std::pair<size_t, uint32_t>>> buckets(T);

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
 * ardal::RoaringMatrix::knnSnv
 *
 * Compute the k nearest neighbours under the SNV Hamming metric.
 *
 * INPUT:
 *  row_idx (size_t) : The index of the query row.
 *  k (int) : Number of neighbours requested.
 *
 * PARAMETERS:
 *  threads (int) : Number of threads for parallel computation.
 *
 * OUTPUT:
 *  py::list : List of (row_index, distance) tuples.
 ****************************************************************************************************/
py::list RoaringMatrix::knnSnv( size_t row_idx,
                                int k,
                                int threads,
                                bool mask_missing ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (row_idx >= n_rows_) {
        throw std::runtime_error("Row index must be non-negative and within bounds.");
    }
    if (k <= 0) {
        return py::list();
    }
    if (!snv_lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    ensure_snv_vectors();
    (void)mask_missing;
    ardal::utils::log_info("Running RoaringMatrix::knnSnv.");

    const std::size_t n = n_rows_;
    const uint32_t k_eff = std::min<uint32_t>(static_cast<uint32_t>(k),
                                              n > 1 ? static_cast<uint32_t>(n - 1) : 0);
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
            for (uint32_t j = 0; j < n; ++j) {
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


/****************************************************************************************************
 * ardal::RoaringMatrix::snvHamming
 *
 * Compute the condensed SNV Hamming distance matrix across all rows.
 *
 * INPUT:
 *  threads (int) : Number of threads for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<uint32_t> : Condensed lower-triangular distance matrix.
 ****************************************************************************************************/
py::array_t<uint32_t> RoaringMatrix::snvHamming( int threads,
                                                 bool mask_missing ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (n_rows_ <= 1) {
        return py::array_t<uint32_t>(static_cast<py::ssize_t>(0));
    }
    if (!snv_lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    ensure_snv_vectors();

    const std::size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<uint32_t> dist_condensed(static_cast<py::ssize_t>(total_pairs));
    auto* dist_ptr = dist_condensed.mutable_data();

    (void)mask_missing;
    ardal::utils::log_info("Running RoaringMatrix::snvHamming.");
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
 * ardal::RoaringMatrix::snvHamming_subset
 *
 * Compute SNV Hamming distances for a subset of rows and columns.
 *
 * INPUT:
 *  row_indices (const std::vector<size_t>&) : Rows to include.
 *  col_indices (const std::vector<size_t>&) : Columns (alleles) to include.
 *
 * PARAMETERS:
 *  threads (int) : Number of threads for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<uint32_t> : Condensed lower-triangular distance matrix.
 ****************************************************************************************************/
py::array_t<uint32_t> RoaringMatrix::snvHamming_subset( const std::vector<size_t>& row_indices,
                                                       const std::vector<size_t>& col_indices,
                                                       int threads,
                                                       bool mask_missing ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (!snv_lookup_loaded_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }

    ensure_snv_vectors();
    ardal::utils::log_info("Running RoaringMatrix::snvHamming_subset.");

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

    (void)mask_missing;
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
                    uint32_t dist_val = use_mask
                        ? snvDistanceRaw(snv_vectors_[lhs_idx], snv_vectors_[rhs_idx], &locus_mask)
                        : snvDistance(lhs_idx, rhs_idx);
                    dist_ptr[idx_base + (j - i - 1)] = dist_val;
                }
            }
    }

    return dist_matrix;
}

} // namespace ardal
