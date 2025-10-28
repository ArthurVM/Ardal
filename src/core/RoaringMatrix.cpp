/*
Copyright 2025 Arthur V. Morris
*/
#include "RoaringMatrix.hpp"
#include "utils/PythonLogger.hpp"
#include <queue>
#include <limits>

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
                              std::shared_ptr<const std::vector<int>> col_masses )
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
py::array_t<uint32_t> RoaringMatrix::hamming( int threads ) const {
    ardal::utils::log_info("Running RoaringMatrix::hamming with scalar instructions.");

    const size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<uint32_t> dist_matrix(total_pairs);

    {
        py::gil_scoped_release release;   // kill python

        auto dist_ptr = dist_matrix.mutable_data();

        #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
        for (size_t i = 0; i < n_rows_; ++i) {
            for (size_t j = i + 1; j < n_rows_; ++j) {
                size_t idx = (i * (2 * n_rows_ - i - 1)) / 2 + (j - i - 1);
                dist_ptr[idx] = hammingDistance(i, j);
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
                                                     int threads ) const {
    ardal::utils::log_info("Running RoaringMatrix::hamming_subset with scalar instructions.");

    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }

    // create matrix mask
    roaring::Roaring col_mask_bitmap;
    for (size_t col_idx : col_indices) {
        if (col_idx >= n_cols_bits_) {
            throw std::out_of_range("Column index in col_indices is out of bounds.");
        }
        col_mask_bitmap.add(col_idx);
    }

    // create a submatrix using the mask and original matrix
    std::vector<roaring::Roaring> submatrix;
    submatrix.reserve(row_indices.size());
    for (size_t row_idx : row_indices) {
        if (row_idx >= n_rows_) {
            throw std::out_of_range("Row index in row_indices is out of bounds.");
        }
        submatrix.push_back(roaring_bitmap_.at(row_idx) & col_mask_bitmap);
    }

    // hamming distance stuff
    const size_t submn_rows_ = submatrix.size();
    if (submn_rows_ == 0) {
        return py::array_t<uint32_t>(0);
    }
    const size_t total_pairs = submn_rows_ * (submn_rows_ - 1) / 2;
    py::array_t<uint32_t> dist_matrix(total_pairs);

    py::gil_scoped_release release;   // release GIL for full parallel region

    #pragma omp parallel num_threads(threads)
    {
        auto dist_ptr = dist_matrix.mutable_data();

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < submn_rows_; ++i) {
            for (size_t j = i + 1; j < submn_rows_; ++j) {
                size_t idx = (i * (2 * submn_rows_ - i - 1)) / 2 + (j - i - 1);
                dist_ptr[idx] = (submatrix[i] ^ submatrix[j]).cardinality();
            }
        }
    }   // end parallel region
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

    const size_t submn_cols_ = col_indices.size();
    const size_t total_pairs = submn_rows_ * (submn_rows_ - 1) / 2;
    py::array_t<double> dist_condensed(static_cast<py::ssize_t>(total_pairs));
    double* dist_ptr = dist_condensed.mutable_data();

    if (submn_cols_ == 0) {
        std::fill(dist_ptr, dist_ptr + total_pairs, 0.0);
        return dist_condensed;
    }

    roaring::Roaring col_mask_bitmap;
    for (size_t col_idx : col_indices) {
        if (col_idx >= n_cols_bits_) {
            throw std::out_of_range("Column index in col_indices is out of bounds.");
        }
        col_mask_bitmap.add(static_cast<uint32_t>(col_idx));
    }

    std::vector<roaring::Roaring> submatrix;
    submatrix.reserve(submn_rows_);
    for (size_t row_idx : row_indices) {
        if (row_idx >= n_rows_) {
            throw std::out_of_range("Row index in row_indices is out of bounds.");
        }
        roaring::Roaring masked = roaring_bitmap_[row_idx];
        masked &= col_mask_bitmap;
        submatrix.emplace_back(std::move(masked));
    }

    std::vector<double> norms(submn_rows_, 0.0);
    for (size_t i = 0; i < submn_rows_; ++i) {
        uint64_t mass = submatrix[i].cardinality();
        norms[i] = mass > 0 ? std::sqrt(static_cast<double>(mass)) : 0.0;
    }

    {
        py::gil_scoped_release release;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (size_t i = 0; i < submn_rows_; ++i) {
            const double norm_i = norms[i];
            const size_t idx_base = (i * (2 * submn_rows_ - i - 1)) / 2;
            const auto& bitmap_i = submatrix[i];

            for (size_t j = i + 1; j < submn_rows_; ++j) {
                const double norm_j = norms[j];
                const double denom = norm_i * norm_j;

                double distance;
                if (denom == 0.0) {
                    distance = (norm_i == 0.0 && norm_j == 0.0) ? 0.0 : 1.0;
                } else {
                    const uint64_t dot = bitmap_i.and_cardinality(submatrix[j]);
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
                std::vector<size_t> cooccurring_partners;
                for (size_t j = i + 1; j < n_cols; ++j) {
                    const auto& j_bitmap = colwise_roaring[j];
                    double intersection_size = (i_bitmap & j_bitmap).cardinality();
                    double union_size = i_bitmap.cardinality() + j_bitmap.cardinality() - intersection_size;

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
    // check if col_indices are valid
    for (const auto& idx : col_indices) {
        if (idx >= n_cols_bits_) {
            throw std::runtime_error("Column index out of bounds.");
        }
    }

    // get the colwise roaring bitmap
    auto colwise_roaring = colwiseRoaringFromRowwise();
    size_t n_cols = colwise_roaring.size();

    std::map<size_t, std::vector<size_t>> global_map;

    {
        py::gil_scoped_release release;   // murder GIL

        #pragma omp parallel num_threads(threads)
        {
            std::map<size_t, std::vector<size_t>> local_map;

            #pragma omp for schedule(static)
            for (size_t i = 0; i < n_cols - 1; ++i) {
                // check if the column is in the subset
                if (std::find(col_indices.begin(), col_indices.end(), i) == col_indices.end()) {
                    continue; // skip columns not in the subset
                }
                
                std::vector<size_t> cooccurring_partners;
                const auto& i_bitmap = colwise_roaring[i];
                for (size_t j = i + 1; j < n_cols; ++j) {
                    // check if the column is in the subset
                    if (std::find(col_indices.begin(), col_indices.end(), j) == col_indices.end()) {
                        continue; // skip columns not in the subset
                    }
                    const auto& j_bitmap = colwise_roaring[j];
                    double intersection_size = (i_bitmap & j_bitmap).cardinality();
                    double union_size = i_bitmap.cardinality() + j_bitmap.cardinality() - intersection_size;
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

} // namespace ardal
