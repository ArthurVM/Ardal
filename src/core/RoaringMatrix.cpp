/*
Copyright 2025 Arthur V. Morris
*/
#include "RoaringMatrix.hpp"
#include "roaring/roaring.hh"
#include "utils/PythonLogger.hpp"

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
RoaringMatrix::RoaringMatrix( std::shared_ptr<const ardal::detail::WordsVV> packed_matrix,
                              std::shared_ptr<const std::vector<int>> row_masses,
                              std::shared_ptr<const std::vector<int>> col_masses,
                              std::size_t n_rows,
                              std::size_t n_cols_bits )
    : n_rows_(n_rows),
      n_cols_(n_cols_bits) 
{
    // allocate memory for roaring bitmaps
    roaring_bitmap_.resize(n_rows_);

    // allocate memory for col_masses and row_masses
    // NOTE: these are passed from HybridMatrix but are causing a core dump when they are used.
    // This is a temporary fixu until I can work out what is going wrong there.
    row_masses_.assign(n_rows_, 0);
    col_masses_.assign(n_cols_, 0);
    
    ardal::utils::log_info("Constructing Roaring bit-map.");
    // populate roaring bitmaps
    for (size_t i = 0; i < n_rows_; ++i) {
        roaring::Roaring& bitmap = roaring_bitmap_[i];
        for (size_t j = 0; j < n_cols_; ++j) {
            uint64_t word = packed_matrix[i][j];
            for (size_t b = 0; b < 64; ++b) {
                size_t col = j * 64 + b;
                if (col >= n_cols_) break;
                if ((word >> b) & 1) {
                    bitmap.add(col);
                    row_masses_[i]++;
                    col_masses_[j]++;
                }
            }
        }
    }
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
        if (col_idx >= n_cols_) {
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
uint32_t RoaringMatrix::hammingDistance( size_t i, size_t j ) const {
    return (roaring_bitmap_[i] ^ roaring_bitmap_[j]).cardinality();
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
    const double intersection_size = (roaring_bitmap_[i] & roaring_bitmap_[j]).cardinality();
    const double union_size = row_masses_[i] + row_masses_[j] - intersection_size;
    // const double union_size = (roaring_bitmap_[i] | roaring_bitmap_[j]).cardinality();

    if (union_size == 0) {
        // if union is 0, both sets are empty and thus identical. Distance is 0.
        return 0.0;
    }
    const double jaccard_index = intersection_size / union_size;
    return jaccard_index;
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
    
    int q_mass = row_masses_[row_idx];    // access pre-calculated row mass
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

                int mass_d = std::abs(q_mass - row_masses_[i]);
                if (mass_d > epsilon) continue;
                
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
            curr_idx++;
        }
    }       
    return ep_n;
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
    std::vector<roaring::Roaring> colwise_roaring(n_cols_);

    // for each row, iterate over its set bits and add the row index to the corresponding column bitmap
    for (size_t row = 0; row < n_rows_; ++row) {
        const auto& row_bitmap = roaring_bitmap_[row];
        // use an iterator to efficiently access set bits
        for (auto it = row_bitmap.begin(); it != row_bitmap.end(); ++it) {
            const size_t col = *it;
            if (col < n_cols_) {
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
    if (colwise_roaring.size() != n_cols_) {
        throw std::runtime_error("col-wise roaring size mismatch");
    }

    // underflow safety
    if (n_cols < 2) return py::dict();


    // uncomment to debug
    std::stringstream ss;
    ss << "n_rows_ = " << n_rows_ << "; " << "n_cols_ = " << n_cols_ << "; " << "roaring_bitmap_.size() = " << roaring_bitmap_.size() << "; " << "colwise_roaring.size() = " << colwise_roaring.size();
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
        if (idx >= n_cols_) {
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