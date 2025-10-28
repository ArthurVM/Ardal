/*
Copyright 2025 Arthur V. Morris
*/
#ifdef __BMI2__
  #include <immintrin.h>
#endif
#include "BitMatrix.hpp"
#include "utils/simd_utils.hpp"
#include "utils/bitops.hpp"
#include "utils/PythonLogger.hpp"
#include "detail/neighbour.hpp"
using ardal::knn::Neighbor;
using ardal::knn::ByMaxDistance;
using ardal::knn::AscDistanceId;

namespace py = pybind11;
namespace ardal {


/****************************************************************************************************
 * ardal::BitMatrix::BitMatrix
 *
 * Constructor for the BitMatrix class.
 *
 * This constructor initializes the BitMatrix object by taking a 2D NumPy array of uint8_t
 * (representing a binary matrix) and converting it into a memory-efficient bit-packed
 * representation (`matrix_`), where each element (0 or 1) is stored as a single bit.
 * It also performs input validation and captures matrix dimensions.
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
 *   - If the matrix dimensions are too large, potentially causing an overflow.
 *   - If the input matrix contains values other than 0 or 1.
 ****************************************************************************************************/
BitMatrix::BitMatrix( ardal::detail::FlatMatrix matrix_,
                      std::shared_ptr<const std::vector<int>> row_masses,
                      std::shared_ptr<const std::vector<int>> col_masses )
  : matrix_(std::move(matrix_)),
    row_masses_(std::move(row_masses)),
    col_masses_(std::move(col_masses)),
    base_(matrix_.base),
    n_rows_(matrix_.n_rows),
    n_cols_bits_(matrix_.n_cols_bits),
    wpr_(matrix_.wpr)
{
    if (!matrix_.base) {
        throw std::runtime_error("BitMatrix: null packed_matrix pointer");
    }

    // generic sanity checks
    if (wpr_ == 0 && n_cols_bits_ != 0) {
        throw std::runtime_error("BitMatrix: computed zero words_per_row unexpectedly");
    }
    if (wpr_ != (n_cols_bits_ + 63) / 64) {
        throw std::runtime_error("BitMatrix: packed_matrix has wrong words_per_row");
    }

    // light sanity on masses
    if (row_masses_->size() != n_rows_) {
        throw std::runtime_error("BitMatrix: row_masses size mismatch");
    }
    if (col_masses_->size() != n_cols_bits_) {
        throw std::runtime_error("BitMatrix: col_masses size mismatch");
    }

    uint64_t tail_mask_ = ardal::tail_mask(n_cols_bits_);
    
    std::stringstream ss;
    ss << "BitMatrix initialised: rows=" << n_rows_ << " wpr=" << wpr_;
    ardal::utils::log_debug(static_cast<string>(ss.str()));
}



/****************************************************************************************************
 * ardal::BitMatrix::getRowSetBitIndices
 *
 * Get the indices of set bits for a given row.
 *
 * This function retrieves the set of column indices that are present (i.e., have a value of 1)
 * in a specified row of the matrix.
 *
 * INPUT:
 *  row_idx (int) : The index of the row in the matrix.
 *
 * OUTPUT:
 *  std::vector<size_t> : A vector containing the indices of the set bits in the specified row.
 *
 * EXCEPTIONS:
 *  std::out_of_range : If the row index is not between 0 and n_rows_
 ****************************************************************************************************/
std::vector<size_t> BitMatrix::getRowSetBitIndices(size_t row_idx) const {
    // input validation
    if (row_idx >= n_rows_) {
        throw std::out_of_range("Row index out of bounds in getRowSetBitIndices.");
    }
    std::vector<size_t> row_indices;
    row_indices.reserve((*row_masses_)[row_idx]);
    const uint64_t* row_ptr = base_ + row_idx * wpr_;
    for (size_t k = 0; k < wpr_; ++k) {
        uint64_t word = row_ptr[k];
        while (word != 0) {
            int trailing_zeros = __builtin_ctzll(word);
            size_t col_idx = k * 64 + trailing_zeros;
            if (col_idx < n_cols_bits_) {
                row_indices.push_back(col_idx);
            }
            word &= word - 1;   // clear the least significant bit
        }
    }
    return row_indices;
}


/****************************************************************************************************
 * ardal::BitMatrix::getColSetBitIndices
 *
 * Get the indices of set bits for a given column.
 *
 * This function retrieves the set of row indices that are present (i.e., have a value of 1)
 * in a specified column of the matrix.
 *
 * INPUT:
 *  col_idx (size_t) : The index of the column in the matrix.
 *
 * OUTPUT:
 *  std::vector<size_t> : A vector containing the indices of the set bits in the specified column.
 *
 * EXCEPTIONS:
 *  std::out_of_range : If the column index is out of bounds (via getBit).
 ****************************************************************************************************/
std::vector<size_t> BitMatrix::getColSetBitIndices( size_t col_idx ) const {
    const auto& col_masses = *col_masses_;
    std::vector<size_t> col_indices;
    col_indices.reserve(col_masses[col_idx]);

    for (size_t i = 0; i < n_rows_; ++i) {
        if (getBit(i, col_idx)) {
            col_indices.push_back(i);
        }
    }
    return col_indices;
}


// /****************************************************************************************************
//  * ardal::BitMatrix::getColumnVector
//  *
//  * Retrieves the bit vector for a specified column.
//  *
//  * This function efficiently extracts all bits from a single column across all rows and returns
//  * them as a new, bit-packed vector. This format is ideal for performing fast bitwise
//  * operations (AND, OR, XOR) between columns for internal C++ use.
//  *
//  * INPUT:
//  *  col_idx (size_t) : The index of the column to retrieve.
//  *
//  * OUTPUT:
//  *  std::vector<uint64_t> : A bit-packed vector representing the column data. The vector's
//  *                           length is `(n_rows_ + 63) / 64`.
//  *
//  * EXCEPTIONS:
//  *  std::out_of_range : If the column index is out of bounds.
//  ****************************************************************************************************/
// std::vector<uint64_t> BitMatrix::getColumnVector( size_t col_idx ) const {
//     if (col_idx >= n_cols_bits_) {
//         throw std::out_of_range("Column index out of bounds.");
//     }

//     size_t packed_rows = (n_rows_ + 63) / 64;
//     std::vector<uint64_t> packed_col_vector(packed_rows, 0);

//     // pre-calculate the memory location of the bit for the given column index
//     // to avoid repeated calculations inside the loop.
//     const size_t src_chunk_idx = col_idx / 64;
//     const uint64_t src_mask = 1ULL << (col_idx % 64);

//     for (size_t i = 0; i < n_rows_; ++i) {
//         // if the bit is set in the source matrix for the given column...
//         if ((matrix_[i][src_chunk_idx] & src_mask)) {
//             // ...set the corresponding bit in our new column vector.
//             packed_col_vector[i / 64] |= (1ULL << (i % 64));
//         }
//     }

//     // apply mask to the last word if n_rows_ is not a multiple of 64
//     size_t remainder_rows = n_rows_ % 64;
//     if (remainder_rows != 0) {
//         packed_col_vector[packed_rows - 1] &= ((1ULL << remainder_rows) - 1ULL);
//     }

//     return packed_col_vector;
// }


/****************************************************************************************************
 * ardal::BitMatrix::getSetBitIndices
 *
 * Get the set bit indices for a given row.
 *
 * This function retrieves the column indices of set bits that are present in a specified row of the
 * matrix and returns them as a NumPy array.
 *
 * INPUT:
 *  row_idx (int) : The index of the row (GUID) in the matrix.
 *
 * OUTPUT:
 *  py::array_t<size_t> : A 1D NumPy array containing the alleles for the specified row.
 ****************************************************************************************************/
py::array_t<size_t> BitMatrix::getSetBitIndices( size_t row_idx ) const {
    std::vector<size_t> row_indices = getRowSetBitIndices(row_idx);
    return py::array_t<size_t>(row_indices.size(), row_indices.data());
}


/****************************************************************************************************
 * ardal::BitMatrix::getBitMatrix
 *
 * Unpack and return the bit matrix as C-contiguous uint8 (0/1).
 *
 * OUTPUT:
 *   py::array_t<uint8_t> of shape (n_rows_, n_cols_bits_)
 ****************************************************************************************************/
py::array_t<uint8_t> BitMatrix::getBitMatrix() const {
    using std::size_t;

    // allocate C contiguous output
    py::array_t<uint8_t> out({n_rows_, n_cols_bits_});
    auto* out_base = static_cast<uint8_t*>(out.request().ptr);

    if (n_rows_ == 0 || n_cols_bits_ == 0) return out;

    const size_t wpr = wpr_;
    const bool has_tail = (n_cols_bits_ & 63u) != 0u;
    const uint64_t tailmask = has_tail ? ardal::tail_mask(n_cols_bits_) : ~0ULL;

    // 256-entry LUT: for each byte b, LUT[b] is a uint64 whose 8 bytes are 0x00/0x01
    // matching the bits of b (LSB -> first byte)
    static const std::array<uint64_t, 256> EXPAND8 = []() {
        std::array<uint64_t, 256> t{};
        for (int b = 0; b < 256; ++b) {
            uint64_t v = 0;
            for (int bit = 0; bit < 8; ++bit)
                if ((b >> bit) & 1) v |= (uint64_t{1} << (bit * 8));  // set whole byte to 0x01
            t[static_cast<size_t>(b)] = v;
        }
        return t;
    }();

    {
        py::gil_scoped_release release;    // kill py

        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < n_rows_; ++i) {
            const uint64_t* row = base_ + i * wpr;            // packed words
            uint8_t* out_row = out_base + i * n_cols_bits_;   // dest start for row i

            for (size_t w = 0; w < wpr; ++w) {
                const size_t col0 = w * 64;
                const size_t remaining = n_cols_bits_ - col0;
                const size_t take = remaining < 64 ? remaining : 64;

                uint64_t x = row[w];
                if (has_tail && (w + 1 == wpr)) x &= tailmask;

                // write full bytes (8 bits -> 8 bytes) via LUT
                const size_t full_bytes = take >> 3;   // 0..8
                uint64_t tmp = x;
                for (size_t b = 0; b < full_bytes; ++b) {
                    const uint8_t byte = static_cast<uint8_t>(tmp & 0xFFu);
                    const uint64_t expanded = EXPAND8[byte];
                    std::memcpy(out_row + col0 + b * 8, &expanded, sizeof(expanded));
                    tmp >>= 8;
                }

                // leftover bits (<8)
                const size_t leftover = take & 7u;
                for (size_t b = 0; b < leftover; ++b) {
                    out_row[col0 + full_bytes * 8 + b] = static_cast<uint8_t>((tmp >> b) & 1u);
                }
            }
        }
    }   // reanimate py

    return out;
}



/****************************************************************************************************
 * ardal::BitMatrix::getPackedMatrix
 * 
 * Return the packed bit matrix.
 * 
 * INPUT: 
 *  None (operates on the private member matrix_)
 *
 * OUTPUT:
 *  py::array_t<uint8_t> : A 2D numpy array representing a packed matrix.
 ****************************************************************************************************/
py::array_t<uint64_t> BitMatrix::getPackedMatrix( void ) const {
    // allocate NumPy array (C-contiguous)
    py::array_t<uint64_t> out({py::ssize_t(n_rows_), py::ssize_t(wpr_)});
    auto buf = out.request();
    auto* dst = static_cast<uint64_t*>(buf.ptr);

    // copy each row
    for (size_t r = 0; r < n_rows_; ++r) {
        const uint64_t* src = base_ + r * wpr_;
        std::memcpy(dst + r * wpr_, src, wpr_ * sizeof(uint64_t));
    }
    return out;
}


/****************************************************************************************************
 * ardal::BitMatrix::hamming
 * 
 * Calculate the Hamming distances between all pairs of rows. SIMD optimised for vectorised distance
 * calculation, with the option of scalar computation for testing and fallback in instances where 
 * CPU does not support AVX2.
 *
 * This function calculates the pairwise Hamming distances between all rows of the allele matrix.
 * The Hamming distance between two rows is the number of positions at which the corresponding
 * elements differ.  The results are returned as a condensed distance matrix.
 *
 * INPUT: 
 *  None (operates on the private member matrix_)
 * 
 * PARAMETERS: 
 *  use_simd (bool)   : A boolean flag to specify whether to use SIMD for vectorised distance
 *                      calculation.
 *  fill_cache (bool) : A boolean flag to specify whether to fill the cache (CURRENTLY INACTIVE).
 *  threads (int)     : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<int>  : A 1D NumPy array representing the condensed distance matrix containing
 *                      the pairwise Hamming distances.  The length of the array is n*(n-1)/2,
 *                      where 'n' is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<uint32_t> BitMatrix::hamming( bool fill_cache,
                                          bool use_simd,
                                          int threads ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }

    const size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<uint32_t> dist_matrix(total_pairs);

    ardal::utils::log_debug("Starting BitMatrix::hamming.");
        
    // dispatch to backend with/without SIMD
    simd_dispatchers::HammingFn hamming_func = simd_dispatchers::hamming_dispatcher(use_simd);
    {
        // open mp thread stuff
        py::gil_scoped_release release;   // release GIL

        auto dist_ptr = dist_matrix.mutable_data();

        #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
        for (size_t i = 0; i < n_rows_; ++i) {
            // get row i
            const uint64_t* __restrict row_i = base_ + i * wpr_;

            // compute the index base once per i
            size_t base = (i * (2 * n_rows_ - i - 1)) / 2;

            for (size_t j = i + 1; j < n_rows_; ++j) {

                // get row j
                const uint64_t* __restrict row_j = base_ + j * wpr_;

                // compute hamming disct
                const uint32_t dist = hamming_func(row_i, row_j, wpr_, n_cols_bits_);
                
                // use the base to construct the index
                dist_ptr[base + (j - i - 1)] = dist;

                // if (fill_cache) {
                //     // currently not in action to save some headaches
                //     #pragma omp critical
                //     _hamming_cache[{i, j}] = dist;
                // }
            }
        }
    }   // GIL reestablished
    return dist_matrix;
}


/****************************************************************************************************
 * ardal::BitMatrix::hamming_subset
 *
 * Calculate the Hamming distances between pairs of rows from a specified subset of rows and columns.
 *
 * This function first creates a sub-matrix based on the provided row and column indices. It then
 * calculates the pairwise Hamming distances between all rows of this sub-matrix. The results
 * are returned as a condensed distance matrix.
 *
 * INPUT:
 *  row_indices (const std::vector<size_t>&) : A vector of row indices to include in the sub-matrix.
 *  col_indices (const std::vector<size_t>&) : A vector of column indices to include in the sub-matrix.
 *
 * PARAMETERS:
 *  use_simd (bool) : A boolean flag to specify whether to use SIMD for vectorised distance calculation.
 *  threads (int)   : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<uint32_t> : A 1D NumPy array representing the condensed distance matrix.
 *
 * EXCEPTIONS:
 *  std::runtime_error : If the thread count is not positive.
 ****************************************************************************************************/
py::array_t<uint32_t> BitMatrix::hamming_subset( const std::vector<size_t>& row_indices,
                                                 const std::vector<size_t>& col_indices,
                                                 bool use_simd,
                                                 int threads ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }

    const size_t submn_rows_ = row_indices.size();
    const size_t submn_cols_ = col_indices.size();
    const size_t submwpr_ = (submn_cols_ + 63) / 64;

    // subset the matrix
    std::vector<std::vector<uint64_t>> submatrix = subsetPackedMatrix(row_indices, col_indices, threads);

    const size_t total_pairs = submn_rows_ * (submn_rows_ - 1) / 2;
    py::array_t<uint32_t> dist_matrix(total_pairs);

    ardal::utils::log_debug("Starting BitMatrix::hamming_subset.");

    // dispatch to backend with/without SIMD
    simd_dispatchers::HammingFn hamming_func = simd_dispatchers::hamming_dispatcher(use_simd);
    {
        // open mp thread stuff
        py::gil_scoped_release release;   // release GIL for full parallel region

        auto dist_ptr = dist_matrix.mutable_data();

        #pragma omp parallel num_threads(threads)
        {
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < submn_rows_; ++i) {
                for (size_t j = i + 1; j < submn_rows_; ++j) {
                    size_t idx = (i * (2 * submn_rows_ - i - 1)) / 2 + (j - i - 1);

                    uint32_t dist = hamming_func(&submatrix[i][0], &submatrix[j][0], submwpr_, n_cols_bits_);

                    dist_ptr[idx] = dist;
                }
            }
        }
    }   // GIL reestablished
    return dist_matrix;
}


/****************************************************************************************************
 * ardal::BitMatrix::innerProduct
 * 
 * Calculate the Inner Product between all pairs of rows. SIMD optimised for vectorised distance
 * calculation, with the option of scalar computation for testing and fallback in instances where 
 * CPU does not support AVX2.
 *
 * This function calculates the pairwise Inner Product between all rows of the bit matrix.
 * The Inner Product between two rows is the number of positions at which the corresponding
 * elements are both set (1). The results are returned as a condensed distance matrix.
 *
 * INPUT: 
 *  None (operates on the private member matrix_)
 * 
 * PARAMETERS: 
 *  use_simd (bool)   : A boolean flag to specify whether to use SIMD for vectorised distance
 *                      calculation.
 *  fill_cache (bool) : A boolean flag to specify whether to fill the cache (CURRENTLY INACTIVE).
 *  threads (int)     : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<int>  : A 1D NumPy array representing the condensed distance matrix containing
 *                      the pairwise Inner Products.  The length of the array is n*(n-1)/2,
 *                      where 'n' is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<int> BitMatrix::innerProduct( bool fill_cache,
                                          bool use_simd,
                                          int threads ) const {
    if (threads <= 0)
        throw std::runtime_error("Thread count must be positive.");
    
    const size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<int> inner_product_matrix(total_pairs);

    ardal::utils::log_debug("Starting BitMatrix::innerProduct.");
        
    // dispatch to backend with/without SIMD
    simd_dispatchers::InnerProductFn inner_product_func = simd_dispatchers::inner_product_dispatcher(use_simd);
    {
        // open mp threading
        py::gil_scoped_release release;   // release GIL for full parallel region

        auto inner_product_ptr = inner_product_matrix.mutable_data();

        #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
        for (size_t i = 0; i < n_rows_; ++i) {
            const uint64_t* __restrict row_i = base_ + i * wpr_;

            // construct the base for the ip results matrix
            size_t ip_ptr_base = (i * (2 * n_rows_ - i - 1)) / 2;

            for (size_t j = i + 1; j < n_rows_; ++j) {    
                const uint64_t* __restrict row_j = base_ + j * wpr_;
                int inner_product = inner_product_func(row_i, row_j, wpr_, n_cols_bits_);

                inner_product_ptr[ip_ptr_base + (j - i - 1)] = inner_product;

                // if (fill_cache) {
                //     // currently not in action to save some headaches
                //     #pragma omp critical
                //     _hamming_cache[{i, j}] = inner_product;
                // }
            }
        }
    }   // GIL reestablished here
    return inner_product_matrix;
}


/****************************************************************************************************
 * ardal::BitMatrix::cosineDistanceAll
 *
 * Compute the all-vs-all cosine distance matrix between rows. Cosine distance is defined as
 * 1 - (dot(x, y) / (||x|| * ||y||)), where rows are treated as binary vectors. Distances are
 * symmetric and the diagonal is zero.
 *
 * INPUT:
 *  None (operates on packed matrix state)
 *
 * PARAMETERS:
 *  use_simd (bool) : Whether to use SIMD kernels for the inner product.
 *  threads (int)   : Number of threads for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<double> : Condensed (n_rows * (n_rows - 1) / 2) vector of pairwise distances.
 ****************************************************************************************************/
py::array_t<double> BitMatrix::cosineDistanceAll( bool use_simd,
                                                  int threads ) const {
    if (threads <= 0) {
        throw std::runtime_error("Thread count must be positive.");
    }
    if (n_rows_ <= 1) {
        return py::array_t<double>(static_cast<py::ssize_t>(0));
    }

    ardal::utils::log_debug("Starting BitMatrix::cosineDistanceAll.");

    const auto& row_masses = *row_masses_;
    std::vector<double> norms(n_rows_);
    for (size_t i = 0; i < n_rows_; ++i) {
        norms[i] = row_masses[i] > 0 ? std::sqrt(static_cast<double>(row_masses[i])) : 0.0;
    }

    const size_t total_pairs = n_rows_ * (n_rows_ - 1) / 2;
    py::array_t<double> dist_condensed(static_cast<py::ssize_t>(total_pairs));
    double* dist_ptr = dist_condensed.mutable_data();

    simd_dispatchers::InnerProductFn inner_product_func = simd_dispatchers::inner_product_dispatcher(use_simd);

    {
        py::gil_scoped_release release;

        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (std::size_t i = 0; i < n_rows_; ++i) {
            const uint64_t* __restrict row_i = base_ + i * wpr_;
            const double norm_i = norms[i];
            const std::size_t idx_base = (i * (2 * n_rows_ - i - 1)) / 2;

            for (std::size_t j = i + 1; j < n_rows_; ++j) {
                const uint64_t* __restrict row_j = base_ + j * wpr_;
                const double norm_j = norms[j];
                const double denom = norm_i * norm_j;

                double distance;
                if (denom == 0.0) {
                    distance = (norm_i == 0.0 && norm_j == 0.0) ? 0.0 : 1.0;
                } else {
                    const int dot = inner_product_func(row_i, row_j, wpr_, n_cols_bits_);
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
 * ardal::BitMatrix::cosineDistance_subset
 *
 * Compute the condensed cosine distance vector for a subset defined by specific rows and columns.
 ****************************************************************************************************/
py::array_t<double> BitMatrix::cosineDistance_subset( const std::vector<size_t>& row_indices,
                                                      const std::vector<size_t>& col_indices,
                                                      bool use_simd,
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

    std::vector<std::vector<uint64_t>> submatrix = subsetPackedMatrix(row_indices, col_indices, threads);
    const size_t submwpr_ = (submn_cols_ + 63) / 64;

    std::vector<double> norms(submn_rows_, 0.0);
    for (size_t i = 0; i < submn_rows_; ++i) {
        uint64_t count = 0;
        const auto& row = submatrix[i];
        for (size_t w = 0; w < submwpr_; ++w) {
            count += ardal::popcnt64(row[w]);
        }
        norms[i] = count > 0 ? std::sqrt(static_cast<double>(count)) : 0.0;
    }

    simd_dispatchers::InnerProductFn inner_product_func = simd_dispatchers::inner_product_dispatcher(use_simd);

    {
        py::gil_scoped_release release;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (size_t i = 0; i < submn_rows_; ++i) {
            const uint64_t* __restrict row_i = submatrix[i].data();
            const double norm_i = norms[i];
            const size_t idx_base = (i * (2 * submn_rows_ - i - 1)) / 2;

            for (size_t j = i + 1; j < submn_rows_; ++j) {
                const uint64_t* __restrict row_j = submatrix[j].data();
                const double norm_j = norms[j];
                const double denom = norm_i * norm_j;

                double distance;
                if (denom == 0.0) {
                    distance = (norm_i == 0.0 && norm_j == 0.0) ? 0.0 : 1.0;
                } else {
                    const int dot = inner_product_func(row_i, row_j, submwpr_, submn_cols_);
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
 * ardal::BitMatrix::neighbourhood
 * 
 * Find the epsilon-neighborhood of a row using Hamming distance.
 *
 * This function identifies the rows in the matrix that are within a specified Hamming distance
 * (epsilon) of a given row.
 *
 * INPUT:
 *  row_coord (size_t) : The index of the query GUID, representing the centroid of the neighbourhood.
 *  epsilon (int)      : The maximum Hamming distance threshold.
 *
 * PARAMETERS:
 *  use_simd (bool) : A boolean flag to specify whether to use SIMD for vectorised distance
 *                      calculation.
 *  threads (int)   : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::list<tuple(row, dist)> : A 1D NumPy array containing the indices of the rows that are within 
 *                               the epsilon-neighborhood of the target row.
 *
 * EXCEPTIONS:
 *  std::runtime_error : If row_coord is out of range.
 ****************************************************************************************************/
py::array_t<int64_t> BitMatrix::neighbourhood( size_t row_idx,
                                               int epsilon,
                                               bool use_simd,
                                               int threads ) const {
    // do some input cleanliness
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
    std::vector<std::vector<std::pair<size_t, int>>> thread_results;
    
    ardal::utils::log_debug("Starting BitMatrix::neighbourhood.");
        
    // dispatch to backend with/without SIMD
    simd_dispatchers::HammingEpsFn epsilon_neighbourhood_func = simd_dispatchers::epsilon_hamming_dispatcher(use_simd);
    {
        // open mp threading
        py::gil_scoped_release release;   // release GIL for full parallel region

        #pragma omp parallel num_threads(threads)
        {
            const int tid = omp_get_thread_num();

            #pragma omp single
            thread_results.resize(omp_get_num_threads());

            // thread local vectors
            std::vector<std::pair<std::size_t,int>> local;
            local.reserve(1024);   // heuristic, dont really need so can be removed

            // get query row
            const uint64_t* __restrict row_q = base_ + row_idx * wpr_;

            #pragma omp for schedule(static)
            for (size_t i = 0; i < n_rows_; ++i) {
                if (i == row_idx) continue;

                int mass_d = std::abs(q_mass - row_masses[i]);
                if (mass_d > epsilon) continue;

                const uint64_t* __restrict row_i = base_ + i * wpr_;

                int distance = epsilon_neighbourhood_func(row_q, row_i, wpr_, n_cols_bits_, epsilon);

                if (distance <= epsilon)
                    local.emplace_back(i, distance);
            }
            thread_results[tid] = std::move(local);
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
 * ardal::BitMatrix::innerProductNeighbourhood
 *
 * Find all rows which share at least ip_epsilon alleles with row_idx.
 *
 * This function identifies GUIDs (rows) in the allele matrix that have a specified minimum
 * number of alleles in common with a target query GUID.
 *
 * INPUT:
 *  row_idx (size_t) : The index of the query GUID.
 *  ip_epsilon (int) : The minimum inner product threshold.
 *
 * OUTPUT:
 *   py::list<tuple(row, ip)> : A 1D NumPy array containing the indices of the rows that are within the 
 *                              IP-neighbourhood of the target row.
 *
 * EXCEPTIONS:
 *  std::runtime_error : If query_guid_index is out of range or k_alleles is negative.
 ****************************************************************************************************/
py::list BitMatrix::innerProductNeighbourhood( size_t row_idx,
                                               int ip_epsilon,
                                               bool use_simd ) const {
     // do some data cleanliness
    if (ip_epsilon < 0) {
        throw std::runtime_error("ip_epsilon must be non-negative.");
        }
    if (row_idx < 0) {  
        throw std::runtime_error("Row index must be non-negative.");
        }
    if (row_idx >= n_rows_) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }

    py::list ipe_n;   // results list

    const auto& row_masses = *row_masses_;

    // if query row has lower mass than ip_epsilon then no matches are possible
    if (static_cast<int>(row_masses[row_idx]) < ip_epsilon) {
        return ipe_n;
    }

    ardal::utils::log_debug("Starting BitMatrix::innerProductNeighbourhood.");
            
    // dispatch to backend with/without SIMD
    simd_dispatchers::InnerProductFn ip_neighbourhood_func = simd_dispatchers::inner_product_dispatcher(use_simd);
    {
        // get query row
        const uint64_t* __restrict row_q = base_ + row_idx * wpr_;

        // iterate through other rows
        for (size_t i = 0; i < n_rows_; ++i) {
            if (i == row_idx) {
                continue;
            }

            // mass optimisation
            if (row_masses[i] < ip_epsilon) {
                continue;
            }

            const uint64_t* __restrict row_i = base_ + i * wpr_;

            int inner_product = ip_neighbourhood_func(row_q, row_i, wpr_, n_cols_bits_);
            
            if (inner_product >= ip_epsilon) {
                ipe_n.append(py::make_tuple(i, inner_product));
            }
        }
    }
    return ipe_n;
}


/****************************************************************************************************
 * ardal::BitMatrix::knn
 *
 * Find the k-nearest-neighbours.
 * NOTE: Currently this does not handle ties.
 *
 * INPUT:
 *  row_idx (size_t) : The index of the query GUID.
 *  k (uint32_t) : The number of neighbours to return.
 *
 * OUTPUT:
 *   py::list<tuple(row, ip)> : A 1D NumPy array containing the indices of the rows which represent
 *                              the nearest neighbours.
 *
 * EXCEPTIONS:
 *  std::runtime_error : If query_guid_index is out of range or k_alleles is negative.
 ****************************************************************************************************/
py::list
BitMatrix::knn_hamming(size_t row_idx, int k, bool use_simd, int threads) const {
    using ardal::knn::Neighbor;
    using ardal::knn::ByMaxDistance;
    ardal::knn::AscDistanceId asc_dist_id{};

    const size_t n = n_rows_;
    if (row_idx >= n || n == 0 || k <= 0 || wpr_ == 0) {
        return py::list();
    }
    const uint32_t k_eff = std::min<uint32_t>(static_cast<uint32_t>(k), n > 1 ? uint32_t(n - 1) : 0);
    if (k_eff == 0) {
        return py::list();
    }

    ardal::utils::log_debug("Starting BitMatrix::knn_hamming.");

    // dispatch to backend with/without SIMD
    auto hamming_eps = ardal::simd_dispatchers::epsilon_hamming_dispatcher(use_simd);

    std::vector<Neighbor> final_neighbors;
    {
        py::gil_scoped_release release;  // release GIL

        const uint64_t* A   = base_ + row_idx * wpr_;
        const uint32_t  a_m = static_cast<uint32_t>((*row_masses_)[row_idx]);

        std::atomic<uint32_t> d_k_global{ std::numeric_limits<uint32_t>::max() };

        const int T = std::max(1, threads);
        std::vector<std::vector<Neighbor>> buckets(T);

        if (T == 1) {
            // serial path
            std::priority_queue<Neighbor, std::vector<Neighbor>, ByMaxDistance> heap;
            uint32_t d_k_local = std::numeric_limits<uint32_t>::max();

            for (uint32_t j = 0; j < n; ++j) {
                if (j == row_idx) continue;

                uint32_t prune_dk = std::numeric_limits<uint32_t>::max();
                if (heap.size() == k_eff) {
                    prune_dk = d_k_local;
                    const uint32_t b_m = static_cast<uint32_t>((*row_masses_)[j]);
                    const uint32_t lbm = (a_m > b_m) ? (a_m - b_m) : (b_m - a_m);
                    if (lbm >= prune_dk) continue;
                }

                const uint64_t* B = base_ + size_t(j) * wpr_;
                const uint32_t eps = (heap.size() == k_eff && prune_dk > 0)
                                   ? (prune_dk - 1)
                                   : (std::numeric_limits<uint32_t>::max() / 2);

                const uint32_t d = hamming_eps(A, B, wpr_, n_cols_bits_, eps);

                if (heap.size() < k_eff) {
                    heap.push({j, d});
                    if (heap.size() == k_eff) d_k_local = heap.top().d;
                } else if (d < d_k_local) {
                    heap.pop(); heap.push({j, d}); d_k_local = heap.top().d;
                }
            }

            std::vector<Neighbor> local; local.reserve(heap.size());
            while (!heap.empty()) { local.push_back(heap.top()); heap.pop(); }
            std::sort(local.begin(), local.end(), asc_dist_id);
            buckets[0] = std::move(local);
        }
        else {
            // parallel path
            #pragma omp parallel num_threads(T)
            {
                std::priority_queue<Neighbor, std::vector<Neighbor>, ByMaxDistance> heap;
                uint32_t d_k_local = std::numeric_limits<uint32_t>::max();

                #pragma omp for schedule(static)
                for (uint32_t j = 0; j < n; ++j) {
                    if (j == row_idx) continue;

                    uint32_t prune_dk = std::numeric_limits<uint32_t>::max();
                    if (heap.size() == k_eff) {
                        const uint32_t g = d_k_global.load(std::memory_order_relaxed);
                        prune_dk = (d_k_local < g) ? d_k_local : g;

                        const uint32_t b_m = static_cast<uint32_t>((*row_masses_)[j]);
                        const uint32_t lbm = (a_m > b_m) ? (a_m - b_m) : (b_m - a_m);
                        if (lbm >= prune_dk) continue;
                    }

                    const uint64_t* B = base_ + size_t(j) * wpr_;
                    const uint32_t eps = (heap.size() == k_eff && prune_dk > 0)
                                       ? (prune_dk - 1)
                                       : (std::numeric_limits<uint32_t>::max() / 2);

                    const uint32_t d = hamming_eps(A, B, wpr_, n_cols_bits_, eps);

                    if (heap.size() < k_eff) {
                        heap.push({j, d});
                        if (heap.size() == k_eff) {
                            d_k_local = heap.top().d;
                            uint32_t cur = d_k_global.load(std::memory_order_relaxed);
                            while (cur > d_k_local &&
                                   !d_k_global.compare_exchange_weak(cur, d_k_local, std::memory_order_relaxed)) {}
                        }
                    } else if (d < d_k_local) {
                        heap.pop(); heap.push({j, d}); d_k_local = heap.top().d;
                        uint32_t cur = d_k_global.load(std::memory_order_relaxed);
                        while (cur > d_k_local &&
                               !d_k_global.compare_exchange_weak(cur, d_k_local, std::memory_order_relaxed)) {}
                    }
                }

                std::vector<Neighbor> local; local.reserve(heap.size());
                while (!heap.empty()) { local.push_back(heap.top()); heap.pop(); }
                std::sort(local.begin(), local.end(), asc_dist_id);

                const int tid =
                #ifdef _OPENMP
                    omp_get_thread_num();
                #else
                    0;
                #endif
                buckets[tid] = std::move(local);
            } // parallel
        }

        // merge from threads
        size_t total = 0; for (auto& v : buckets) total += v.size();
        std::vector<Neighbor> all; all.reserve(total);
        for (auto& v : buckets) all.insert(all.end(), v.begin(), v.end());

        if (all.size() > k_eff) {
            std::nth_element(all.begin(), all.begin() + k_eff, all.end(), asc_dist_id);
            all.resize(k_eff);
        }
        std::sort(all.begin(), all.end(), asc_dist_id);

        final_neighbors = std::move(all);
    } // GIL reacquired

    py::list result;
    for (const auto& nb : final_neighbors) {
        result.append(py::make_tuple(static_cast<int64_t>(nb.id),
                                     static_cast<int64_t>(nb.d)));
    }
    return result;
}


py::list
BitMatrix::knn_inner_product(size_t row_idx, int k, bool use_simd, int threads) const {
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
    if (row_idx >= n || n == 0 || k <= 0 || wpr_ == 0) {
        return py::list();
    }
    const uint32_t k_eff = std::min<uint32_t>(static_cast<uint32_t>(k), n > 1 ? uint32_t(n - 1) : 0);
    if (k_eff == 0) {
        return py::list();
    }

    ardal::utils::log_debug("Starting BitMatrix::knn_inner_product.");

    auto inner_product = ardal::simd_dispatchers::inner_product_dispatcher(use_simd);

    std::vector<Candidate> final_candidates;
    {
        py::gil_scoped_release release;

        const uint64_t* A = base_ + row_idx * wpr_;

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

                const uint64_t* B = base_ + size_t(j) * wpr_;
                const uint32_t score = inner_product(A, B, wpr_, n_cols_bits_);

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
        for (auto& v : buckets) total += v.size();
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
        result.append(py::make_tuple(static_cast<int64_t>(cand.id),
                                     static_cast<int64_t>(cand.score)));
    }
    return result;
}


py::list
BitMatrix::knn_jaccard(size_t row_idx, int k, bool use_simd, int threads) const {
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
    if (row_idx >= n || n == 0 || k <= 0 || wpr_ == 0) {
        return py::list();
    }
    const uint32_t k_eff = std::min<uint32_t>(static_cast<uint32_t>(k), n > 1 ? uint32_t(n - 1) : 0);
    if (k_eff == 0) {
        return py::list();
    }

    ardal::utils::log_debug("Starting BitMatrix::knn_jaccard.");

    auto inner_product = ardal::simd_dispatchers::inner_product_dispatcher(use_simd);

    const uint32_t mass_query = static_cast<uint32_t>((*row_masses_)[row_idx]);

    std::vector<Candidate> final_candidates;
    {
        py::gil_scoped_release release;

        const uint64_t* A = base_ + row_idx * wpr_;

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

                const uint64_t* B = base_ + size_t(j) * wpr_;
                const uint32_t intersection = inner_product(A, B, wpr_, n_cols_bits_);
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
        for (auto& v : buckets) total += v.size();
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


py::list
BitMatrix::knn_cosine(size_t row_idx, int k, bool use_simd, int threads) const {
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
    if (row_idx >= n || n == 0 || k <= 0 || wpr_ == 0) {
        return py::list();
    }
    const uint32_t k_eff = std::min<uint32_t>(static_cast<uint32_t>(k), n > 1 ? uint32_t(n - 1) : 0);
    if (k_eff == 0) {
        return py::list();
    }

    ardal::utils::log_debug("Starting BitMatrix::knn_cosine.");

    auto inner_product = ardal::simd_dispatchers::inner_product_dispatcher(use_simd);

    const uint32_t mass_query = static_cast<uint32_t>((*row_masses_)[row_idx]);
    const double norm_query = mass_query > 0 ? std::sqrt(static_cast<double>(mass_query)) : 0.0;

    std::vector<Candidate> final_candidates;
    {
        py::gil_scoped_release release;

        const uint64_t* A = base_ + row_idx * wpr_;

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

                const uint64_t* B = base_ + size_t(j) * wpr_;
                const uint32_t intersection = inner_product(A, B, wpr_, n_cols_bits_);
                const uint32_t mass_other = static_cast<uint32_t>((*row_masses_)[j]);
                const double norm_other = mass_other > 0 ? std::sqrt(static_cast<double>(mass_other)) : 0.0;

                double score = 0.0;
                const double denom = norm_query * norm_other;
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
        for (auto& v : buckets) total += v.size();
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
 * ardal::BitMatrix::uniqueSharedBits
 *
 * Finds the indices of bits that are set (1) in ALL specified rows and are NOT set in ANY other row.
 *
 * INPUT:
 *  row_indices (const std::vector<size_t>&) : A vector of row indices (sample IDs) to consider as the group.
 * 
 * PARAMETERS:
 *  use_simd (bool) : A boolean flag to specify whether to use SIMD for vectorised bitwise operations.
 *
 * OUTPUT:
 *  std::vector<size_t> : A vector containing the column indices of the unique shared bits.
 ****************************************************************************************************/
std::vector<size_t> BitMatrix::uniqueSharedBits( const std::vector<size_t>& row_indices,
                                                 bool use_simd ) const {
    if (row_indices.empty()) return {};
 
    const size_t ingroup_size = row_indices.size();
    std::vector<uint64_t> group_and(wpr_, ~0ULL);   // initialize with all bits sets

    // TODO: this is a bit of a mess
    if (use_simd) {
        ardal::utils::log_debug("Running BitMatrix::uniqueSharedBits with AVX2 instructions.");
        // SIMD ingroup AND
        size_t k = 0;
        for (; k + 4 <= wpr_; k += 4) {
            // get query row
            const uint64_t* __restrict row_k = base_ + k * wpr_;
            __m256i acc = _mm256_loadu_si256((__m256i const*)row_k);
            for (size_t idx = 1; idx < ingroup_size; ++idx) {
                const uint64_t* __restrict row_idx = base_ + idx * wpr_;
                __m256i row = _mm256_loadu_si256((__m256i const*)row_idx);
                acc = _mm256_and_si256(acc, row);
            }
            _mm256_storeu_si256((__m256i*)&group_and[k], acc);
        }
        // remainder loop for ingroup AND
        for (; k < wpr_; ++k) {
            uint64_t acc = *(base_ + k * wpr_);
            for (size_t idx = 1; idx < ingroup_size; ++idx) {
                acc &= *(base_ + idx * wpr_);
            }
            group_and[k] = acc;
        }
    } else {
        ardal::utils::log_debug("Running BitMatrix::uniqueSharedBits with scalar instructions.");
        // scalar ingroup AND
        for (size_t k = 0; k < wpr_; ++k) {
            uint64_t acc = *(base_ + k * wpr_);
            for (size_t idx = 1; idx < ingroup_size; ++idx) {
                acc &= *(base_ + idx * wpr_);
            }
            group_and[k] = acc;
        }
    }
    
    // check column mass for uniqueness
    // this checks each column mass (number of set bits) against the number of ingroup rows
    // in order for that column to be unique and shared by the in group, these two values must be identical
    const auto& col_masses = *col_masses_;
    std::vector<size_t> unique_bits;
    for (size_t k = 0; k < wpr_; ++k) {
        uint64_t chunk = group_and[k];
        if (chunk == 0) continue;   // optimisation: skip empty chunks

        while (chunk != 0) {
            int trailing_zeros = __builtin_ctzll(chunk);
            size_t col_idx = k * 64 + trailing_zeros;
            if (col_idx < n_cols_bits_) {
                if (col_idx < n_cols_bits_ && col_masses[col_idx] == ingroup_size) {
                    unique_bits.push_back(col_idx);
                }
            }
            chunk &= chunk - 1;   // clear the LSB
        }
    }
    return unique_bits;
}


/****************************************************************************************************
 * ardal::BitMatrix::density
 *
 * Calculates the density of the matrix, defined as the proportion of set bits.
 *
 * INPUT:
 *  None
 *
 * OUTPUT:
 *  double : The density of the matrix.
 ****************************************************************************************************/
double BitMatrix::density( void ) const {
    if (n_rows_ == 0 || n_cols_bits_ == 0) {
        return 0.0;
    }
    const auto& col_masses = *col_masses_;
    double total_set_bits = std::accumulate(col_masses.begin(), col_masses.end(), 0.0);
    return total_set_bits / (static_cast<double>(n_rows_) * n_cols_bits_);
}


/****************************************************************************************************
 * ardal::BitMatrix::colFrequency
 *
 * Calculates the frequency of set bits for each column within a specified subset of rows.
 *
 * INPUT:
 *  row_indices (std::vector<size_t>&) : A vector of row indices to consider for frequency calculation.
 *
 * OUTPUT:
 *  py::array_t<double> : A 1D NumPy array of row frequencies, one for each column.
 ****************************************************************************************************/
py::array_t<double> BitMatrix::colFrequency( std::vector<size_t>& row_indices ) const {
    py::array_t<double> frequencies(n_cols_bits_);
    auto frequencies_ptr = frequencies.mutable_data();
    std::fill(frequencies_ptr, frequencies_ptr + n_cols_bits_, 0.0);

    const size_t n = row_indices.size();
    if (n_cols_bits_ == 0) return frequencies;

    const auto& col_masses = *col_masses_;
    // if no indices provided or all rows requested, return normalized column masses
    if (n == 0 || n == n_rows_) {
        for (size_t col = 0; col < n_cols_bits_; ++col) {
            frequencies_ptr[col] = static_cast<double>(col_masses[col]) / static_cast<double>(n_rows_);
        }
        return frequencies;
    }

    {
        py::gil_scoped_release release;

        for (size_t row : row_indices) {
            if (row >= n_rows_) continue;
            const uint64_t* row_ptr = base_ + row * wpr_;
            for (size_t k = 0; k < wpr_; ++k) {
                uint64_t chunk = row_ptr[k];
                while (chunk) {
                    int tz = __builtin_ctzll(chunk);
                    const size_t col = k * 64 + static_cast<size_t>(tz);
                    if (col < n_cols_bits_) frequencies_ptr[col] += 1.0;
                    chunk &= chunk - 1;
                }
            }
        }
    } // GIL restored

    const double denom = static_cast<double>(n);
    for (size_t col = 0; col < n_cols_bits_; ++col) {
        frequencies_ptr[col] /= denom;
    }

    return frequencies;
}


/****************************************************************************************************
 * ardal::BitMatrix::bitCooccurrence
 * 
 * NOTE: This function is much slower than the RoaringMatrix version, so it is currently not used.
 * It may be revisited in the future if the RoaringMatrix version is not suitable for some reason.
 *
 * Calculates the co-occurrence of bits across columns in the matrix, returning a dictionary of
 * column indices and their co-occurring partners.
 *
 * INPUT:
 *  threshold (double)   : The threshold for co-occurrence.
 *  use_simd (bool)      : A boolean flag to specify whether to use SIMD for vectorised distance
 *                         calculation.
 *  threads (int)        : The number of threads to use for parallel computation.
 *  cache_bytes (size_t) : The size of the cache in bytes for storing column vectors.
 *
 * OUTPUT:
 *  py::dict : A dictionary where keys are column indices and values are lists
 *             of co-occurring column indices.
 ****************************************************************************************************/
// py::dict BitMatrix::bitCooccurrence(double threshold, bool use_simd, int threads, size_t cache_bytes) const {
//     // do some input cleanliness
//     if (threshold < 0) {
//         throw std::runtime_error("threshold must be non-negative.");
//     }
//     if (threads <= 0) {
//         throw std::runtime_error("Number of threads must be positive.");
//     }
//     if (n_rows_ == 0) {
//         return py::dict();
//     }

//     // construct column cache
//     size_t packed_rows = (n_rows_ + 63) / 64;
//     const double max_disagreements = (1.0 - threshold) * static_cast<double>(n_rows_);

//     // calculate cache size
//     if (cache_bytes == 0) {
//         cache_bytes = packed_rows * n_cols_bits_ * sizeof(uint64_t);
//     } else if (cache_bytes < packed_rows * n_cols_bits_ * sizeof(uint64_t)) {
//         throw std::runtime_error("Cache size is too small for the number of columns and rows.");
//     }
//     // distribute the cache size evenly across threads
//     if (cache_bytes % threads != 0) {
//         cache_bytes -= cache_bytes % threads;  // make it divisible by threads
//     }
//     int thread_cache_bytes = cache_bytes / threads;

//     // COMMENT OUT FOR DEBUG
//     std::cout << "Calculating bit co-occurrences with threshold: " << threshold
//               << " (missing = " << max_disagreements << ")" << std::endl;

//     std::map<size_t, std::vector<size_t>> global_map;
//     std::unordered_set<size_t> global_visited;

//     auto start_time = std::chrono::high_resolution_clock::now();

//     #pragma omp parallel num_threads(threads)
//     {
//         ColumnCache local_cache(thread_cache_bytes, packed_rows);
//         std::map<size_t, std::vector<size_t>> local_map;
//         std::unordered_set<size_t> local_visited;

//         #pragma omp for schedule(dynamic)
//         for (size_t i = 0; i < n_cols_bits_; ++i) {
//             if (local_visited.count(i)) {
//                 // std::cout << "0. Skipping column " << i << " as it has already been visited." << std::endl;
//                 continue;
//             }

//             // std::cout << "Processing column " << i << " with mass " << col_masses_[i] << std::endl;

//             const std::vector<uint64_t>& i_vec = local_cache.get(i, [this](size_t idx) {
//                 return getColumnVector(idx);
//             });

//             std::vector<size_t> i_ref_cooccur_vec;
//             size_t i_mass = col_masses_[i];
//             size_t valid_tail = n_rows_ % 64;

//             for (size_t j = i + 1; j < n_cols_bits_; ++j) {
//                 if (local_visited.count(j)) {
//                     // std::cout << "1. Skipping column " << j << " as it has already been visited." << std::endl;
//                     continue;
//                 }

//                 // col mass optimisation
//                 if (std::abs(static_cast<long>(i_mass) - static_cast<long>(col_masses_[j])) > max_disagreements) {
//                     // std::cout << "2. Skipping pair (" << i << ", " << j << ") due to mass difference." << std::endl;
//                     continue;
//                 }

//                 const std::vector<uint64_t>& j_vec = local_cache.get(j, [this](size_t idx) {
//                     return getColumnVector(idx);
//                 });

//                 // calculate jaccard
//                 size_t intersection_size = 0;
//                 size_t union_size = 0;

//                 popcount_pairwise(i_vec, j_vec, valid_tail, use_simd,
//                                   intersection_size, union_size);

//                 if ((union_size - intersection_size) > max_disagreements) {
//                     // std::cout << "3. Skipping pair (" << i << ", " << j << ") due to early exit." << std::endl;
//                     continue;
//                 }

//                 if (union_size == 0) {
//                     // std::cout << "4. Skipping pair (" << i << ", " << j << ") due to zero union size." << std::endl;
//                     continue;
//                 }

//                 // std::cout << "5. Intersection: " << intersection_size << ", Union size: " << union_size << std::endl;

//                 double jaccard = static_cast<double>(intersection_size) / union_size;
//                 if (jaccard >= threshold) {
//                     i_ref_cooccur_vec.push_back(j);
//                     local_visited.insert(j);
//                 }
//             }

//             if (!i_ref_cooccur_vec.empty()) {
//                 local_map[i] = std::move(i_ref_cooccur_vec);
//                 local_visited.insert(i);
//             } else {
//                 // std::cout << "6. No co-occurrences found for column " << i << "." << std::endl;
//             }

//             local_cache.evict(i);
//         }

//         #pragma omp critical
//         {
//             for (auto& [k, v] : local_map) global_map[k] = std::move(v);
//             global_visited.insert(local_visited.begin(), local_visited.end());
//         }
//     }

//     auto end_time = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed = end_time - start_time;
//     std::cout << "[bitCooccurrence] Time elapsed: " << elapsed.count() << " seconds" << std::endl;

//     // Convert to Python dictionary
//     py::dict cooccurrences_py;
//     for (auto& [k, v] : global_map) {
//         cooccurrences_py[py::int_(k)] = py::cast(v);
//     }

//     return cooccurrences_py;
// }



/****************************************************************************************************
 * ardal::BitMatrix::columnEntropy
 *
 * Calculates the Shannon entropy for each column in the matrix.
 * Entropy H(X) = -p * log2(p) - (1-p) * log2(1-p), where p is the bit frequency.
 *
 * INPUT:
 *  None
 *
 * OUTPUT:
 *  py::array_t<double> : A 1D NumPy array of entropy values, one for each column.
 ****************************************************************************************************/
py::array_t<double> BitMatrix::columnEntropy( void ) const {
    ardal::utils::log_debug("Running BitMatrix::columnEntropy.");

    const auto& col_masses = *col_masses_;
    py::array_t<double> entropies(n_cols_bits_);
    auto entropies_ptr = entropies.mutable_data();

    for (size_t j = 0; j < n_cols_bits_; ++j) {
        if (n_rows_ == 0) {
            entropies_ptr[j] = 0.0;
            continue;
        }
        double p = static_cast<double>(col_masses[j]) / n_rows_;

        if (p == 0.0 || p == 1.0) {
            entropies_ptr[j] = 0.0;
        } else {
            entropies_ptr[j] = -p * std::log2(p) - (1.0 - p) * std::log2(1.0 - p);
        }
    }
    return entropies;
}


/****************************************************************************************************
 * ardal::BitMatrix::klDivergence
 *
 * Calculates the Kullback-Leibler (KL) divergence for each column, measuring how well it
 * distinguishes an "ingroup" of rows from the "outgroup" (all other rows).
 *
 * INPUT:
 *  ingroup_indices (const std::vector<size_t>&) : A vector of row indices for the ingroup.
 *
 * OUTPUT:
 *  py::array_t<double> : A 1D NumPy array of KL divergence scores, one for each column.
 ****************************************************************************************************/
py::array_t<double> BitMatrix::klDivergence( const std::vector<size_t>& ingroup_indices ) const {
    ardal::utils::log_debug("Running BitMatrix::klDivergence.");

    const size_t ingroup_size = ingroup_indices.size();
    if (ingroup_size == 0 || ingroup_size == n_rows_) {
        // if ingroup is empty or all rows then discrimination is meaningless
        py::array_t<double> scores(n_cols_bits_);
        std::fill(scores.mutable_data(), scores.mutable_data() + n_cols_bits_, 0.0);
        return scores;
    }
    const size_t outgroup_size = n_rows_ - ingroup_size;

    // calculate column masses for the ingroup
    // could probably be optimised at some point
    std::vector<int> ingroupcol_masses_(n_cols_bits_, 0);
    for (const auto& row_idx : ingroup_indices) {
        if (row_idx >= n_rows_) continue;
        const uint64_t* row_ptr = base_ + row_idx * wpr_;
        for (size_t k = 0; k < wpr_; ++k) {
            uint64_t chunk = row_ptr[k];
            while (chunk != 0) {
                int trailing_zeros = __builtin_ctzll(chunk);
                size_t col_idx = k * 64 + trailing_zeros;
                if (col_idx < n_cols_bits_) {
                    ingroupcol_masses_[col_idx]++;
                }
                chunk &= chunk - 1;   // clear the LSB
            }
        }
    }

    py::array_t<double> scores(n_cols_bits_);
    auto scores_ptr = scores.mutable_data();
    const auto& col_masses = *col_masses_;

    for (size_t j = 0; j < n_cols_bits_; ++j) {
        double ingroup_mass = static_cast<double>(ingroupcol_masses_[j]);
        double outgroup_mass = static_cast<double>(col_masses[j] - ingroupcol_masses_[j]);

        // Laplace smoothing to avoid log(0) or division by 0 evil
        double p_ingroup = (ingroup_mass + 1.0) / (ingroup_size + 2.0);
        double p_outgroup = (outgroup_mass + 1.0) / (outgroup_size + 2.0);

        double q_ingroup = 1.0 - p_ingroup;
        double q_outgroup = 1.0 - p_outgroup;

        double dkl = 0.0;
        if (p_ingroup > 0.0 && p_outgroup > 0.0) {
            dkl += p_ingroup * std::log2(p_ingroup / p_outgroup);
        }
        if (q_ingroup > 0.0 && q_outgroup > 0.0) {
            dkl += q_ingroup * std::log2(q_ingroup / q_outgroup);
        }
        scores_ptr[j] = dkl;
    }
    return scores;
}


/****************************************************************************************************
 * ardal::BitMatrix::jsDivergence
 *
 * Calculates the Jensen-Shannon divergence for each column, measuring how well it distinguishes an
 * "ingroup" of rows from the "outgroup" (all other rows).
 *
 * INPUT:
 *  ingroup_indices (const std::vector<size_t>&) : A vector of row indices for the ingroup.
 *
 * OUTPUT:
 *  py::array_t<double> : A 1D NumPy array of Jensen-Shannon divergence scores, one for each column.
 ****************************************************************************************************/
py::array_t<double> BitMatrix::jsDivergence( const std::vector<size_t>& ingroup_indices ) const {
    ardal::utils::log_debug("Running BitMatrix::jsDivergence.");

    const size_t ingroup_size = ingroup_indices.size();
    if (ingroup_size == 0 || ingroup_size == n_rows_) {
        py::array_t<double> scores(n_cols_bits_);
        std::fill(scores.mutable_data(), scores.mutable_data() + n_cols_bits_, 0.0);
        return scores;
    }
    const size_t outgroup_size = n_rows_ - ingroup_size;

    // calculate column masses for the ingroup
    // could probably be optimised at some point
    std::vector<int> ingroupcol_masses_(n_cols_bits_, 0);
    for (const auto& row_idx : ingroup_indices) {
        if (row_idx >= n_rows_) continue;
        const uint64_t* row_ptr = base_ + row_idx * wpr_;
        for (size_t k = 0; k < wpr_; ++k) {
            uint64_t chunk = row_ptr[k];
            while (chunk) {
                int tz = __builtin_ctzll(chunk);
                size_t col = k * 64 + tz;
                if (col < n_cols_bits_) {
                    ingroupcol_masses_[col]++;
                }
                chunk &= chunk - 1;    // clear the LSB
            }
        }
    }

    py::array_t<double> scores(n_cols_bits_);
    auto scores_ptr = scores.mutable_data();
    const auto& col_masses = *col_masses_;

    for (size_t j = 0; j < n_cols_bits_; ++j) {
        double ig_mass = static_cast<double>(ingroupcol_masses_[j]);
        double og_mass = static_cast<double>(col_masses[j] - ingroupcol_masses_[j]);

        // Laplace smoothing
        double p = (ig_mass + 1.0) / (ingroup_size + 2.0);
        double q = (og_mass + 1.0) / (outgroup_size + 2.0);

        double m = 0.5 * (p + q);
        double js = 0.0;

        if (p > 0.0) js += 0.5 * p * std::log2(p / m);
        if (q > 0.0) js += 0.5 * q * std::log2(q / m);

        double mp = 1.0 - p, mq = 1.0 - q, mm = 1.0 - m;
        if (mp > 0.0) js += 0.5 * mp * std::log2(mp / mm);
        if (mq > 0.0) js += 0.5 * mq * std::log2(mq / mm);

        scores_ptr[j] = js;
    }
    return scores;
}


/****************************************************************************************************
 * ardal::BitMatrix::informationGain
 *
 * Calculates the information gain for each column in the matrix with respect to a specified ingroup.
 * Information gain is calculated as the difference between the entropy of the class and the conditional
 * entropy given the presence of a specific allele (column).
 *
 * INPUT:
 *  ingroup_indices (const std::vector<size_t>&) : A vector of row indices for the ingroup.
 *
 * OUTPUT:
 *  py::array_t<double> : A 1D NumPy array of information gain scores, one for each column.
 ****************************************************************************************************/
py::array_t<double> BitMatrix::informationGain( const std::vector<size_t>& ingroup_indices ) const {
    ardal::utils::log_debug("Running BitMatrix::informationGain.");

    const size_t ingroup_size = ingroup_indices.size();
    if (ingroup_size == 0 || ingroup_size == n_rows_) {
        py::array_t<double> scores(n_cols_bits_);
        std::fill(scores.mutable_data(), scores.mutable_data() + n_cols_bits_, 0.0);
        return scores;
    }

    // calculate column masses for the ingroup
    // could probably be optimised at some point
    std::vector<int> ingroupcol_masses_(n_cols_bits_, 0);
    for (const auto& row_idx : ingroup_indices) {
        if (row_idx >= n_rows_) continue;
        const uint64_t* row_ptr = base_ + row_idx * wpr_;
        for (size_t k = 0; k < wpr_; ++k) {
            uint64_t chunk = row_ptr[k];
            while (chunk) {
                int tz = __builtin_ctzll(chunk);
                size_t col = k * 64 + tz;
                if (col < n_cols_bits_) {
                    ingroupcol_masses_[col]++;
                }
                chunk &= chunk - 1;    // clear the LSB
            }
        }
    }

    const size_t total_size = n_rows_;
    const double p_class = static_cast<double>(ingroup_size) / total_size;

    auto entropy = [](double p) -> double {
        return (p > 0.0 && p < 1.0) ? -p * std::log2(p) - (1.0 - p) * std::log2(1.0 - p) : 0.0;
    };

    double H_C = entropy(p_class);

    py::array_t<double> scores(n_cols_bits_);
    auto scores_ptr = scores.mutable_data();
    const auto& col_masses = *col_masses_;

    for (size_t j = 0; j < n_cols_bits_; ++j) {
        size_t n_11 = static_cast<size_t>(ingroupcol_masses_[j]);
        size_t n_10 = static_cast<size_t>(col_masses[j] - ingroupcol_masses_[j]);
        size_t n_01 = ingroup_size - n_11;
        size_t n_00 = (n_rows_ - ingroup_size) - n_10;

        double H_C_given_snp = 0.0;

        double n0 = n_00 + n_01;
        if (n0 > 0) {
            double p0 = static_cast<double>(n_01) / n0;
            H_C_given_snp += (n0 / static_cast<double>(total_size)) * entropy(p0);
        }

        double n1 = n_10 + n_11;
        if (n1 > 0) {
            double p1 = static_cast<double>(n_11) / n1;
            H_C_given_snp += (n1 / static_cast<double>(total_size)) * entropy(p1);
        }

        scores_ptr[j] = H_C - H_C_given_snp;
    }
    return scores;
}


/****************************************************************************************************
 * ardal::BitMatrix::subsetPackedMatrix
 *
 * Creates a new bit-packed matrix from a subset of rows and columns.
 *
 * This is a helper function that extracts specified rows and columns from the main `matrix_`
 * and constructs a new, smaller, bit-packed matrix. This is used by functions that operate on
 * a subset of the data, like `hamming_subset`.
 *
 * INPUT:
 *  row_indices (const std::vector<size_t>&) : A vector of original row indices to include.
 *  col_indices (const std::vector<size_t>&) : A vector of original column indices to include.
 *
 * OUTPUT:
 *  std::vector<std::vector<uint64_t>> : The new, smaller bit-packed matrix.
 ****************************************************************************************************/
#include <vector>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <algorithm>

std::vector<std::vector<uint64_t>>
BitMatrix::subsetPackedMatrix(const std::vector<size_t>& row_indices,
                              const std::vector<size_t>& col_indices,
                              int threads) const
{
    // Assume: base_ points to a row-major packed bit matrix:
    //   nrows_ rows, wpr_ 64-bit words per row (so ncols = 64*wpr_).
    // Each row r lives at base_ + r*wpr_.

    const size_t R = row_indices.size();
    const size_t C = col_indices.size();
    const size_t W = (C + 63u) / 64u;     // words per output row

    std::vector<std::vector<uint64_t>> out(R, std::vector<uint64_t>(W, 0ULL));
    if (R == 0 || C == 0) return out;

    const size_t ncols = 64u * wpr_;

    // --- light runtime validation (cheap and prevents OOB) ---
    for (size_t r : row_indices) {
        if (r >= n_rows_) throw std::out_of_range("row_indices contains out-of-range row");
    }
    for (size_t c : col_indices) {
        if (c >= ncols) throw std::out_of_range("col_indices contains out-of-range column");
    }

    // Precompute source/dest word/mask per selected column
    std::vector<uint32_t> src_w(C), dst_w(C);
    std::vector<uint64_t> src_mask(C), dst_mask(C);
    for (size_t j = 0; j < C; ++j) {
        const size_t c = col_indices[j];
        src_w[j]    = static_cast<uint32_t>(c >> 6);          // which input word
        src_mask[j] = (1ULL << (c & 63u));                    // bit in that word
        dst_w[j]    = static_cast<uint32_t>(j >> 6);          // which output word
        dst_mask[j] = (1ULL << (j & 63u));                    // bit in that word
    }

    const bool has_tail = (C & 63u) != 0u;
    const uint64_t tailmask = has_tail ? ((1ULL << (C & 63u)) - 1ULL) : ~0ULL;

    // threads <= 0 -> single-thread
    const int nthreads = (threads > 0 ? threads : 1);

    // If you expose this via pybind11, you can add GIL release here.
    // pybind11::gil_scoped_release release;

#if defined(_OPENMP)
    #pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
    for (size_t ri = 0; ri < R; ++ri) {
        const size_t r = row_indices[ri];
        const uint64_t* __restrict src = base_ + r * wpr_;
        uint64_t* __restrict dst = out[ri].data();

        // Generic path: test & scatter
        for (size_t j = 0; j < C; ++j) {
            if (src[src_w[j]] & src_mask[j]) {
                dst[dst_w[j]] |= dst_mask[j];
            }
        }

        // Zero the slack bits in the final word
        if (has_tail) {
            dst[W - 1] &= tailmask;
        }
    }

    return out;
}

 
double BitMatrix::getDensity( void ) const {
    return density_;
}


size_t BitMatrix::getNCols( void ) const {
    return n_cols_bits_;
}


size_t BitMatrix::getNRows( void ) const {
    return n_rows_;
}


/****************************************************************************************************
 * ardal::BitMatrix::getRowMasses
 *
 * Get the popcount of each row in the matrix as a NumPy array.
 *
 * This function returns the pre-calculated popcount (mass) of each row in the matrix. The popcount
 * of a row represents the number of set bits (1s) present in that row.
 *
 * INPUT:
 *  None (operates on the private member row_masses_)
 *
 * OUTPUT:
 *  py::array_t<int> : A NumPy array containing the popcount of each row in the matrix.
 ****************************************************************************************************/
py::array_t<int> BitMatrix::getRowMasses( void ) const {
    const auto& row_masses = *row_masses_;
    return py::array_t<int>(row_masses.size(), row_masses.data());
}


/****************************************************************************************************
 * ardal::BitMatrix::getRowMassesVector
 *
 * Get the popcount of each row in the matrix as a std::vector.
 *
 * This function returns a const reference to the internal vector storing the pre-calculated
 * popcount (mass) of each row. This is primarily for efficient internal use, such as passing
 * data to the RoaringMatrix backend without extra copies.
 *
 * INPUT:
 *  None (operates on the private member row_masses_)
 *
 * OUTPUT:
 *  const std::vector<int>& : A const reference to the vector of row masses.
 ****************************************************************************************************/
const std::vector<int>& BitMatrix::getRowMassesVector( void ) const {
    return *row_masses_;
}


/****************************************************************************************************
 * ardal::BitMatrix::getColumnMasses
 *
 * Get the popcount of each column in the matrix.
 *
 * This function returns the pre-calculated popcount of each column in the matrix. The popcount of a
 * column represents the number of set bits (1s) present in that column.
 *
 * INPUT:
 *  None (operates on the private member col_masses_)
 *
 * OUTPUT:
 *  py::array_t<int> : A NumPy array containing the popcount of each column in the matrix.
 ****************************************************************************************************/
py::array_t<int> BitMatrix::getColumnMasses( void ) const {
    const auto& col_masses = *col_masses_;
    return py::array_t<int>(col_masses.size(), col_masses.data());
}


/****************************************************************************************************
 * ardal::BitMatrix::getSubsetMatrix_rows
 *
 * Get a subset of the packed matrix as a NumPy array. This function is executed if only row ids
 * are provided, since this can be done rapidly.
 *
 *
 * INPUT:
 *  row_indices (const std::vector<size_t>&) : A vector of row indices for the subset.
 *
 * OUTPUT:
 *  py::array_t<uint64_t> : A 2D NumPy array representing the packed sub-matrix.
 ****************************************************************************************************/
py::array_t<uint64_t> BitMatrix::getSubsetMatrix_rows( const std::vector<size_t>& row_indices,
                                                       int threads) const {
    ardal::utils::log_debug("BitMatrix::getSubsetMatrix_rows.");
    const size_t R = row_indices.size();
    const size_t W = wpr_;                 // 64-bit words per row

    py::array_t<uint64_t> out({R, W});
    if (R == 0 || W == 0) return out;

    // basic bounds checks
    #ifndef NDEBUG
    for (size_t r : row_indices) {
        if (r >= n_rows_) throw std::out_of_range("row_indices contains out-of-range row");
    }
    #endif

    uint64_t* out_ptr = out.mutable_data();

    const bool has_tail = (n_cols_bits_ & 63u) != 0u;
    const uint64_t tail = has_tail ? ((uint64_t{1} << (n_cols_bits_ & 63u)) - 1u) : ~uint64_t{0};

    // release GIL
    py::gil_scoped_release release;

    const int nthreads = threads > 0 ? threads : 1;

    #if defined(_OPENMP)
    #pragma omp parallel for schedule(static) num_threads(nthreads)
    #endif
    for (size_t i = 0; i < R; ++i) {
        const size_t r = row_indices[i];
        const uint64_t* src = base_ + r * W;
        uint64_t* dst = out_ptr + i * W;
        std::memcpy(dst, src, W * sizeof(uint64_t));
        if (has_tail) dst[W - 1] &= tail; // ensure slack bits are zeroed
    }

    return out;
}



/****************************************************************************************************
 * ardal::BitMatrix::getSubsetMatrix
 *
 * Get a subset of the packed matrix as a NumPy array.
 *
 * This function creates a sub-matrix from the specified row and column indices and returns it
 * as a 2D NumPy array of uint64_t, preserving the bit-packed format. This is mainly for
 * debugging or specialized external use.
 *
 * INPUT:
 *  row_indices (const std::vector<size_t>&) : A vector of row indices for the subset.
 *  col_indices (const std::vector<size_t>&) : A vector of column indices for the subset.
 *
 * OUTPUT:
 *  py::array_t<uint64_t> : A 2D NumPy array representing the packed sub-matrix.
 ****************************************************************************************************/
py::array_t<uint64_t> BitMatrix::getSubsetMatrix( const std::vector<size_t>& row_indices, 
                                                  const std::vector<size_t>& col_indices,
                                                  const int threads ) const {
    ardal::utils::log_debug("BitMatrix::getSubsetMatrix.");
    const std::vector<std::vector<uint64_t>> submat = subsetPackedMatrix(row_indices, col_indices, threads);
    const size_t nrows = row_indices.size();
    const size_t ncols = col_indices.size();
    const size_t wpr = (ncols + 63) / 64;

    py::array_t<uint64_t> packed_matrix({py::ssize_t(nrows), py::ssize_t(wpr)});
    if (nrows == 0 || wpr == 0) {
        return packed_matrix;
    }

    uint64_t* out_ptr = packed_matrix.mutable_data();
    for (size_t i = 0; i < nrows; ++i) {
        const std::vector<uint64_t>& row = submat[i];
        std::memcpy(out_ptr + i * wpr, row.data(), wpr * sizeof(uint64_t));
    }

    return packed_matrix;
}

} // namespace ardal
