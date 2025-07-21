/*
Copyright 2025 Arthur V. Morris
*/
#include "BitMatrix.hpp"
#include "ColumnCache.hpp"
#include "simd_popcount.hpp"
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <iostream>

namespace py = pybind11;
namespace _ardal {


/****************************************************************************************************
 * ardal::BitMatrix::BitMatrix
 *
 * Constructor for the BitMatrix class.
 *
 * This constructor initializes the BitMatrix object by taking a 2D NumPy array of uint8_t
 * (representing a binary matrix) and converting it into a memory-efficient bit-packed
 * representation (`_packed_matrix`), where each element (0 or 1) is stored as a single bit.
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
BitMatrix::BitMatrix(py::array_t<uint8_t> input_matrix) {
    auto buf = input_matrix.request();

    if (buf.ndim != 2) {
        throw std::runtime_error("Input matrix must be 2-dimensional");
    }

    // capture matrix dimensions
    _n_rows = buf.shape[0];
    _n_cols = buf.shape[1];

    // calculate packed matrix columns
    _packed_cols = (_n_cols + 63) / 64;

    // size check
    if (_n_rows > std::numeric_limits<size_t>::max() / _n_cols) {
        throw std::runtime_error("Matrix dimensions are too large, potential overflow.");
    }

    // allocate memory for packed matrix
    auto ptr = static_cast<uint8_t*>(buf.ptr);
    _packed_matrix.resize(_n_rows, std::vector<uint64_t>(_packed_cols, 0));
    
    // get the mass (popcount) of each row and column
    _row_masses.resize(_n_rows, 0);
    _col_masses.resize(_n_cols, 0);

    // do some packing
    for (size_t i = 0; i < _n_rows; ++i) {
        for (size_t j = 0; j < _n_cols; ++j) {
            uint8_t val = ptr[i * _n_cols + j];
            if (val != 0 && val != 1) {
                throw std::runtime_error("Input matrix must only contain binary values (0 or 1)");
            }
            if (val) {
                _packed_matrix[i][j / 64] |= (1ULL << (j % 64));
                // populate mass counters
                _row_masses[i]++;
                _col_masses[j]++;
            }
        }
    }
    // calculate matrix density
    _density = density();
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
 *  std::out_of_range : If the row index is not between 0 and _n_rows
 ****************************************************************************************************/
std::vector<size_t> BitMatrix::getRowSetBitIndices( size_t row_idx ) const {
    // input validation
    if (row_idx >= _n_rows || row_idx < 0) {
        throw std::out_of_range("Row index out of bounds in getRowSetBitIndices.");
    }
    std::vector<size_t> row_indices;
    row_indices.reserve(_row_masses[row_idx]);
    for (size_t k = 0; k < _packed_cols; ++k) {
        uint64_t chunk = _packed_matrix[row_idx][k];
        while (chunk != 0) {
            int trailing_zeros = __builtin_ctzll(chunk);
            size_t col_idx = k * 64 + trailing_zeros;
            if (col_idx < _n_cols) {
                row_indices.push_back(col_idx);
            }
            chunk &= chunk - 1;   // clear the least significant bit
        }
    }
    return row_indices;
}


std::vector<size_t> BitMatrix::getColSetBitIndices( size_t col_idx ) const {
    std::vector<size_t> col_indices;
    col_indices.reserve(_col_masses[col_idx]);

    for (size_t i = 0; i < _n_rows; ++i) {
        if (getBit(i, col_idx)) {
            col_indices.push_back(i);
        }
    }
    return col_indices;
}


/****************************************************************************************************
 * ardal::BitMatrix::getColumnVector
 *
 * Retrieves the bit vector for a specified column.
 *
 * This function efficiently extracts all bits from a single column across all rows and returns
 * them as a new, bit-packed vector. This format is ideal for performing fast bitwise
 * operations (AND, OR, XOR) between columns for internal C++ use.
 *
 * INPUT:
 *  col_idx (size_t) : The index of the column to retrieve.
 *
 * OUTPUT:
 *  std::vector<uint64_t> : A bit-packed vector representing the column data. The vector's
 *                           length is `(_n_rows + 63) / 64`.
 *
 * EXCEPTIONS:
 *  std::out_of_range : If the column index is out of bounds.
 ****************************************************************************************************/
std::vector<uint64_t> BitMatrix::getColumnVector(size_t col_idx) const {
    if (col_idx >= _n_cols) {
        throw std::out_of_range("Column index out of bounds.");
    }

    size_t packed_rows = (_n_rows + 63) / 64;
    std::vector<uint64_t> packed_col_vector(packed_rows, 0);

    // Pre-calculate the memory location of the bit for the given column index
    // to avoid repeated calculations inside the loop.
    const size_t src_chunk_idx = col_idx / 64;
    const uint64_t src_mask = 1ULL << (col_idx % 64);

    for (size_t i = 0; i < _n_rows; ++i) {
        // If the bit is set in the source matrix for the given column...
        if ((_packed_matrix[i][src_chunk_idx] & src_mask)) {
            // ...set the corresponding bit in our new column vector.
            packed_col_vector[i / 64] |= (1ULL << (i % 64));
        }
    }
    return packed_col_vector;
}


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
 * Unpack and return the bit matrix.
 * 
 * INPUT: 
 *  None (operates on the private member _packed_matrix)
 *
 * OUTPUT:
 *  py::array_t<uint8_t> : A 2D numpy array representing a binary matrix.
 ****************************************************************************************************/
py::array_t<uint8_t> BitMatrix::getBitMatrix( void ) const {
    py::array_t<uint8_t> unpacked_matrix({_n_rows, _n_cols});
    auto unpacked_ptr = static_cast<uint8_t*>(unpacked_matrix.request().ptr);

    for (size_t i = 0; i < _n_rows; ++i) {
        for (size_t j = 0; j < _n_cols; ++j) {
            unpacked_ptr[i * _n_cols + j] = getBit(i, j);
        }
    }
    return unpacked_matrix; 
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
 *  None (operates on the private member _packed_matrix)
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
py::array_t<uint32_t> BitMatrix::hamming( bool fill_cache, bool use_simd, int threads ) const {
    if (threads <= 0)
        throw std::runtime_error("Thread count must be positive.");
    
    const size_t total_pairs = _n_rows * (_n_rows - 1) / 2;
    py::array_t<uint32_t> dist_matrix(total_pairs);
    
    {
        py::gil_scoped_release release;   // release GIL for full parallel region
        omp_set_num_threads(threads);

        auto dist_ptr = dist_matrix.mutable_data();

        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < _n_rows; ++i) {
            for (size_t j = i + 1; j < _n_rows; ++j) {
                size_t idx = (i * (2 * _n_rows - i - 1)) / 2 + (j - i - 1);
                uint32_t dist = use_simd
                    ? hammingDistanceSIMD(i, j)
                    : hammingDistanceScalar(i, j);

                dist_ptr[idx] = dist;

                // if (fill_cache) {
                //     // currently not in action to save some headaches
                //     #pragma omp critical
                //     _hamming_cache[{i, j}] = dist;
                // }
            }
        }
    }   // GIL reestablished here
    return dist_matrix;
}



/****************************************************************************************************
 * ardal::BitMatrix::hammingDistanceScalar
 *
 * Calculates the Hamming distance between two rows of the bit-packed matrix using scalar operations.
 * The Hamming distance is the number of positions at which the corresponding bits differ.
 * This function iterates through each byte of the packed rows, performs a bitwise XOR, and then
 * counts the set bits (popcount) in the result.
 *
 * This is a private helper function used by `hamming`.
 *
 * INPUT:
 *  i (size_t) : The row index of the first allele sequence.
 *  j (size_t) : The row index of the second allele sequence.
 *
 * OUTPUT:
 *  int : The Hamming distance between the two specified rows.
 * 
 * EXCEPTIONS:
 *  std::out_of_range : If i or j are not smaller than _n_rows
 ****************************************************************************************************/
uint32_t BitMatrix::hammingDistanceScalar( size_t i, size_t j ) const {
    // input validation in the name of paranoia/robustness
    if (i >= _n_rows || j >= _n_rows) {
        throw std::out_of_range("Row index out of bounds in hammingDistanceScalar.");
    }

    uint32_t dist = 0;
    for (size_t k = 0; k < _packed_cols; ++k) {
        uint64_t xor_byte = _packed_matrix[i][k] ^ _packed_matrix[j][k];
        dist += __builtin_popcountll(xor_byte);
    }
    return dist;
}


/****************************************************************************************************
 * ardal::BitMatrix::hammingDistanceSIMD
 *
 * Calculates the Hamming distance between two rows of the bit-packed matrix using SIMD (AVX2)
 * intrinsics for optimized performance.
 * The function processes the packed bytes in chunks of 32 bytes (256 bits) using AVX2 instructions
 * to perform bitwise XOR and then uses `__builtin_popcount` on 32-bit lanes to sum the differing bits.
 * A remainder loop handles any bytes not fitting into 32-byte chunks.
 *
 * This is a private helper function used by `hamming`.
 *
 * INPUT:
 *  i (size_t) : The row index of the first allele sequence.
 *  j (size_t): The row index of the second allele sequence.
 *
 * OUTPUT:
 *  int : The Hamming distance between the two specified rows.
 * 
 * EXCEPTIONS:
 *  std::out_of_range : If i or j are not smaller than _n_rows
 ****************************************************************************************************/
uint32_t BitMatrix::hammingDistanceSIMD( size_t i, size_t j ) const {
    // input valiation for robustness
    if (i >= _n_rows || j >= _n_rows) {
        throw std::out_of_range("Row index out of bounds in hammingDistanceSIMD.");
    }

    uint32_t dist = 0;
    size_t k = 0;

    // SIMD block
    for (; k + 4 <= _packed_cols; k += 4) {
        __m256i a = _mm256_loadu_si256((__m256i const*)&_packed_matrix[i][k]);
        __m256i b = _mm256_loadu_si256((__m256i const*)&_packed_matrix[j][k]);
        __m256i xor_result = _mm256_xor_si256(a, b);

        // horizontal popcount on 4x 64-bit integers
        uint64_t chunks[4];
        _mm256_storeu_si256((__m256i*)chunks, xor_result);
        dist += __builtin_popcountll(chunks[0]) + __builtin_popcountll(chunks[1]) +
                __builtin_popcountll(chunks[2]) + __builtin_popcountll(chunks[3]);
    }

    // cover the non 64 chunk tail
    for (; k < _packed_cols; ++k) {
        uint64_t xor_chunk = _packed_matrix[i][k] ^ _packed_matrix[j][k];
        dist += __builtin_popcountll(xor_chunk);
    }
    return dist;
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
 *  None (operates on the private member _packed_matrix)
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
py::array_t<int> BitMatrix::innerProduct( bool fill_cache, bool use_simd, int threads ) const {
    if (threads <= 0)
        throw std::runtime_error("Thread count must be positive.");
    
    const size_t total_pairs = _n_rows * (_n_rows - 1) / 2;
    py::array_t<int> inner_product_matrix(total_pairs);
    
    {
        py::gil_scoped_release release;   // release GIL for full parallel region

        omp_set_num_threads(threads);

        auto inner_product_ptr = inner_product_matrix.mutable_data();

        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < _n_rows; ++i) {
            for (size_t j = i + 1; j < _n_rows; ++j) {
                size_t idx = (i * (2 * _n_rows - i - 1)) / 2 + (j - i - 1);
                int inner_product = use_simd
                    ? innerProductSIMD(i, j)
                    : innerProductScalar(i, j);

                inner_product_ptr[idx] = inner_product;

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
py::array_t<int64_t> BitMatrix::neighbourhood( size_t row_idx, int epsilon, bool use_simd, int threads ) const {
    // do some input cleanliness
    if (epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_idx < 0) {
        throw std::runtime_error("Row index must be non-negative.");
        }
    if (row_idx >= _n_rows) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }
    if (threads <= 0)
        throw std::runtime_error("Number of threads must be positive.");
    
    int q_mass = _row_masses[row_idx];    // access pre-calculated row mass
    std::vector<std::vector<std::pair<size_t, int>>> thread_results;

    {
        py::gil_scoped_release release;   // release GIL for full parallel region
        omp_set_num_threads(threads);     // Explicitly control number of threads

        #pragma omp parallel
        {
            const int thread_id = omp_get_thread_num();

            #pragma omp single
            thread_results.resize(omp_get_num_threads());

            auto& local = thread_results[thread_id];

            #pragma omp for schedule(static)
            for (size_t i = 0; i < _n_rows; ++i) {
                if (i == row_idx) continue;

                int mass_d = std::abs(q_mass - _row_masses[i]);
                if (mass_d > epsilon) continue;

                int distance = use_simd
                    ? epsilonNeighbourhoodSIMD(row_idx, i, epsilon)
                    : epsilonNeighbourhoodScalar(row_idx, i, epsilon);

                if (distance <= epsilon)
                    local.emplace_back(i, distance);
            }
        }
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
 * ardal::BitMatrix::epsilonNeighbourhoodScalar
 *
 * Calculates the Hamming distance between two rows of the bit-packed matrix and assesses whether
 * they exist in the same neighbourhood, performing an early exit if epsilon is exceeded. 
 * Scalar implementation for for testing and fallback in instances where CPU does not support AVX2.
 * 
 * The function processes the packed bytes in chunks of 32 bytes (256 bits) using AVX2 instructions
 * to perform bitwise XOR and then uses `__builtin_popcount` on 32-bit lanes to sum the differing bits.
 * A remainder loop handles any bytes not fitting into 32-byte chunks.
 *
 * This is a private helper function used by `neighbourhood`.
 *
 * INPUT:
 *  i (size_t)       : The row index of the first allele sequence.
 *  j (size_t)       : The row index of the second allele sequence.
 *  epsilon (int)    : The size of the neighbourhood (in hamming distance).
 *
 * OUTPUT:
 *  int : The Hamming distance between the two specified rows.
 * 
 * EXCEPTIONS:
 *  std::out_of_range : If i or j are not smaller than _n_rows
 ****************************************************************************************************/
int BitMatrix::epsilonNeighbourhoodScalar( size_t i, size_t j, int epsilon ) const {
    // input valiation for robustness
    if (i >= _n_rows || j >= _n_rows) {
        throw std::out_of_range("Row index out of bounds in hammingDistanceSIMD.");
    }

    int dist = 0;
    for (size_t k = 0; k < _packed_cols; ++k) {
        uint64_t xor_chunk = _packed_matrix[i][k] ^ _packed_matrix[j][k];
        dist += __builtin_popcountll(xor_chunk);
        if (dist > epsilon) {
            return dist;  // early break where epsilon is exceeded
        }
    }
    return dist;
}


/****************************************************************************************************
 * ardal::BitMatrix::epsilonNeighbourhoodSIMD
 *
 * Calculates the Hamming distance between two rows of the bit-packed matrix and assesses whether
 * they exist in the same neighbourhood, performing an early exit if epsilon is exceeded. 
 * Vectorised with SIMD (AVX2) intrinsics for optimized performance.
 * The function processes the packed bytes in chunks of 32 bytes (256 bits) using AVX2 instructions
 * to perform bitwise XOR and then uses `__builtin_popcount` on 32-bit lanes to sum the differing bits.
 * A remainder loop handles any bytes not fitting into 32-byte chunks.
 *
 * This is a private helper function used by `neighbourhood`.
 *
 * INPUT:
 *  i (size_t)       : The row index of the first allele sequence.
 *  j (size_t)       : The row index of the second allele sequence.
 *  epsilon (size_t) : The size of the neighbourhood (in hamming distance).
 *
 * OUTPUT:
 *  int : The Hamming distance between the two specified rows.
 * 
 * EXCEPTIONS:
 *  std::out_of_range : If i or j are not smaller than _n_rows
 ****************************************************************************************************/
int BitMatrix::epsilonNeighbourhoodSIMD( size_t i, size_t j, int epsilon ) const {
    // input valiation for robustness
    if (i >= _n_rows || j >= _n_rows) {
        throw std::out_of_range("Row index out of bounds in hammingDistanceSIMD.");
    }

    int dist = 0;
    size_t k = 0;

    // SIMD block
    for (; k + 4 <= _packed_cols; k += 4) {
        __m256i a = _mm256_loadu_si256((__m256i const*)&_packed_matrix[i][k]);
        __m256i b = _mm256_loadu_si256((__m256i const*)&_packed_matrix[j][k]);
        __m256i xor_result = _mm256_xor_si256(a, b);

        // horizontal popcount on 4x 64bit integers
        uint64_t chunks[4];
        _mm256_storeu_si256((__m256i*)chunks, xor_result);
        dist += __builtin_popcountll(chunks[0]) + __builtin_popcountll(chunks[1]) +
                __builtin_popcountll(chunks[2]) + __builtin_popcountll(chunks[3]);

        if (dist > epsilon) {
            return dist;   // early break where epsilon is exceeded
        }
    }

    // cover the non 32 chunk tail
    for (; k < _packed_cols; ++k) {
        uint64_t xor_chunk = _packed_matrix[i][k] ^ _packed_matrix[j][k];
        dist += __builtin_popcountll(xor_chunk);
    }
    return dist;
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
py::list BitMatrix::innerProductNeighbourhood( size_t row_idx, int ip_epsilon, bool use_simd ) const {
     // do some data cleanliness
    if (ip_epsilon < 0) {
        throw std::runtime_error("ip_epsilon must be non-negative.");
        }
    if (row_idx < 0) {
        throw std::runtime_error("Row index must be non-negative.");
        }
    if (row_idx >= _n_rows) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }

    py::list ipe_n;   // results list

    // if query row has lower mass than ip_epsilon then no matches are possible
    if (static_cast<int>(_row_masses[row_idx]) < ip_epsilon) {
        return ipe_n;
    }

    // iterate through other rows
    for (size_t i = 0; i < _n_rows; ++i) {
        if (i == row_idx) {
            continue;
        }

        // mass optimisation
        if (_row_masses[i] < ip_epsilon) {
            continue;
        }

        int inner_product;   // initialise inner product
        if (use_simd) {
            inner_product = innerProductSIMD(row_idx, i);
        } else {
            inner_product = innerProductScalar(row_idx, i);
        }
        
        if (inner_product >= ip_epsilon) {
            ipe_n.append(py::make_tuple(i, inner_product));
        }
    }
    return ipe_n;
}


/****************************************************************************************************
 * ardal::BitMatrix::innerProductScalar
 *
 * Calculates the inner product (number of shared set bits) between two rows of the bit-packed matrix
 * using scalar operations.
 * This function iterates through each byte of the packed rows, performs a bitwise AND, and then
 * counts the set bits (popcount) in the result.
 *
 * This is a private helper function used by `innerProductNeighbourhood`.
 *
 * INPUT:
 *  i (size_t) : The row index of the first allele sequence.
 *  j (size_t) : The row index of the second allele sequence.
 *
 * OUTPUT:
 *  int : The inner product between the two specified rows.
 *
 * EXCEPTIONS:
 *  std::out_of_range : If i or j are not smaller than _n_rows
 ****************************************************************************************************/
 
int BitMatrix::innerProductScalar( size_t i, size_t j ) const {
    int ip = 0;
    for (size_t k = 0; k < _packed_cols; ++k) {
        uint64_t and_chunk = _packed_matrix[i][k] & _packed_matrix[j][k];
        ip += __builtin_popcountll(and_chunk);
    }
    return ip;
}


/****************************************************************************************************
 * ardal::BitMatrix::innerProductSIMD
 *
 * Calculates the inner product (number of shared set bits) between two rows of the bit-packed matrix
 * using SIMD (AVX2) intrinsics for optimized performance.
 * The function processes the packed bytes in chunks of 32 bytes (256 bits) using AVX2 instructions
 * to perform bitwise AND and then uses `__builtin_popcount` on 32-bit lanes to sum the common bits.
 * A remainder loop handles any bytes not fitting into 32-byte chunks.
 *
 * This is a private helper function used by `innerProductNeighbourhood`.
 *
 * INPUT:
 *  i (size_t) : The row index of the first allele sequence.
 *  j (size_t) : The row index of the second allele sequence.
 *
 * OUTPUT:
 *  int : The inner product between the two specified rows.
 *
 * EXCEPTIONS:
 *  std::out_of_range : If i or j are not smaller than _n_rows
 ****************************************************************************************************/
 
int BitMatrix::innerProductSIMD( size_t i, size_t j ) const {
    int ip = 0;
    size_t k = 0;

    // SIMD block
    for (; k + 4 <= _packed_cols; k += 4) {
        __m256i a = _mm256_loadu_si256((__m256i const*)&_packed_matrix[i][k]);
        __m256i b = _mm256_loadu_si256((__m256i const*)&_packed_matrix[j][k]);
        __m256i and_result = _mm256_and_si256(a, b);

        // horizontal popcount on 4x 64-bit integers
        uint64_t chunks[4];
        _mm256_storeu_si256((__m256i*)chunks, and_result);
        ip += __builtin_popcountll(chunks[0]) + __builtin_popcountll(chunks[1]) +
              __builtin_popcountll(chunks[2]) + __builtin_popcountll(chunks[3]);
    }

    // cover the non 32 chunk tail
    for (; k < _packed_cols; ++k) {
        uint64_t and_chunk = _packed_matrix[i][k] & _packed_matrix[j][k];
        ip += __builtin_popcountll(and_chunk);
    }
    return ip;
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
std::vector<size_t> BitMatrix::uniqueSharedBits(const std::vector<size_t>& row_indices, bool use_simd) const {
    if (row_indices.empty()) return {};
 
    const size_t ingroup_size = row_indices.size();
    std::vector<uint64_t> group_and(_packed_cols, ~0ULL);   // initialize with all bits set

    if (use_simd) {
        // SIMD ingroup AND
        size_t k = 0;
        for (; k + 4 <= _packed_cols; k += 4) {
            __m256i acc = _mm256_loadu_si256((__m256i const*)&_packed_matrix[row_indices[0]][k]);
            for (size_t idx = 1; idx < ingroup_size; ++idx) {
                __m256i row = _mm256_loadu_si256((__m256i const*)&_packed_matrix[row_indices[idx]][k]);
                acc = _mm256_and_si256(acc, row);
            }
            _mm256_storeu_si256((__m256i*)&group_and[k], acc);
        }
        // remainder loop for ingroup AND
        for (; k < _packed_cols; ++k) {
            uint64_t acc = _packed_matrix[row_indices[0]][k];
            for (size_t idx = 1; idx < ingroup_size; ++idx) {
                acc &= _packed_matrix[row_indices[idx]][k];
            }
            group_and[k] = acc;
        }
    } else {
        // scalar ingroup AND
        for (size_t k = 0; k < _packed_cols; ++k) {
            uint64_t acc = _packed_matrix[row_indices[0]][k];
            for (size_t idx = 1; idx < ingroup_size; ++idx) {
                acc &= _packed_matrix[row_indices[idx]][k];
            }
            group_and[k] = acc;
        }
    }

    // check column mass for uniqueness
    // this checks each column mass (number of set bits) against the number of ingroup rows
    // in order for that column to be unique and shared by the in group, these two values must be identical
    std::vector<size_t> unique_bits;
    for (size_t k = 0; k < _packed_cols; ++k) {
        uint64_t chunk = group_and[k];
        if (chunk == 0) continue;   // optimisation: skip empty chunks

        while (chunk != 0) {
            int trailing_zeros = __builtin_ctzll(chunk);
            size_t col_idx = k * 64 + trailing_zeros;
            if (col_idx < _n_cols) {
                if (col_idx < _n_cols && _col_masses[col_idx] == ingroup_size) {
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
    if (_n_rows == 0 || _n_cols == 0) {
        return 0.0;
    }
    double total_set_bits = std::accumulate(_col_masses.begin(), _col_masses.end(), 0.0);
    return total_set_bits / (static_cast<double>(_n_rows) * _n_cols);
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
py::array_t<double> BitMatrix::colFrequency(std::vector<size_t>& row_indices) const {
    py::array_t<double> frequencies(_n_cols);
    auto frequencies_ptr = frequencies.mutable_data();

    std::fill(frequencies_ptr, frequencies_ptr + _n_cols, 0.0);

    size_t n = row_indices.size();
    if (_n_cols == 0) {
        return frequencies;
    }

    // if no or all indices are provided then quickly return normalised column masses
    if (n == _n_cols || n == 0) {
        for (size_t col = 0; col < _n_cols; ++col) {
            frequencies_ptr[col] = _col_masses[col] / static_cast<double>(_n_rows);
        }
        return frequencies;
    }

    for (size_t row : row_indices) {
        if (row >= _n_rows) continue;
        for (size_t col = 0; col < _n_cols; ++col) {
            if (getBit(row, col)) {
                frequencies_ptr[col] += 1.0;
            }
        }
    }

    for (size_t col = 0; col < _n_cols; ++col) {
        frequencies_ptr[col] /= static_cast<double>(n);
    }

    return frequencies;
}



py::dict BitMatrix::bitCooccurrence(double threshold, bool use_simd, int threads, size_t cache_bytes) const {
    if (threshold < 0) {
        throw std::runtime_error("threshold must be non-negative.");
    }
    if (threads <= 0) {
        throw std::runtime_error("Number of threads must be positive.");
    }
    if (_n_rows == 0) {
        return py::dict();
    }

    // construct column cache
    size_t packed_rows = (_n_rows + 63) / 64;
    const double max_disagreements = (1.0 - threshold) * static_cast<double>(_n_rows);
    ColumnCache col_cache(cache_bytes, packed_rows);

    // COMMENT OUT FOR DEBUG
    std::cout << "Calculating bit co-occurrences with threshold: " << threshold
              << " (missing = " << max_disagreements << ")" << std::endl;

    std::map<size_t, std::vector<size_t>> global_map;
    std::unordered_set<size_t> global_visited;

    auto start_time = std::chrono::high_resolution_clock::now();

    #pragma omp parallel num_threads(threads)
    {
        std::map<size_t, std::vector<size_t>> local_map;
        std::unordered_set<size_t> local_visited;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < _n_cols; ++i) {
            if (local_visited.count(i)) {
                // std::cout << "0. Skipping column " << i << " as it has already been visited." << std::endl;
                continue;
            }

            std::cout << "Processing column " << i << " with mass " << _col_masses[i] << std::endl;

            const std::vector<uint64_t>& i_vec = col_cache.get(i, [this](size_t idx) {
                return getColumnVector(idx);
            });

            std::vector<size_t> i_ref_cooccur_vec;
            size_t i_mass = _col_masses[i];
            size_t valid_tail = _n_rows % 64;

            for (size_t j = i + 1; j < _n_cols; ++j) {
                if (local_visited.count(j)) {
                    // std::cout << "1. Skipping column " << j << " as it has already been visited." << std::endl;
                    continue;
                }

                // col mass optimisation
                if (std::abs(static_cast<long>(i_mass) - static_cast<long>(_col_masses[j])) > max_disagreements) {
                    // std::cout << "2. Skipping pair (" << i << ", " << j << ") due to mass difference." << std::endl;
                    continue;
                }

                const std::vector<uint64_t>& j_vec = col_cache.get(j, [this](size_t idx) {
                    return getColumnVector(idx);
                });

                // Calculate co-occurrence (Jaccard index)
                size_t intersection_size = 0;
                size_t union_size = 0;

                popcount_pairwise(i_vec, j_vec, valid_tail, use_simd,
                                  intersection_size, union_size);

                if ((union_size - intersection_size) > max_disagreements) {
                    // std::cout << "3. Skipping pair (" << i << ", " << j << ") due to early exit." << std::endl;
                    continue;
                }

                if (union_size == 0) {
                    // std::cout << "4. Skipping pair (" << i << ", " << j << ") due to zero union size." << std::endl;
                    continue;
                }

                // std::cout << "5. Intersection: " << intersection_size << ", Union size: " << union_size << std::endl;

                double jaccard = static_cast<double>(intersection_size) / union_size;
                if (jaccard >= threshold) {
                    i_ref_cooccur_vec.push_back(j);
                    local_visited.insert(j);
                }
            }

            if (!i_ref_cooccur_vec.empty()) {
                local_map[i] = std::move(i_ref_cooccur_vec);
                local_visited.insert(i);
            } else {
                // std::cout << "6. No co-occurrences found for column " << i << "." << std::endl;
            }

            col_cache.evict(i);
        }

        #pragma omp critical
        {
            for (auto& [k, v] : local_map) global_map[k] = std::move(v);
            global_visited.insert(local_visited.begin(), local_visited.end());
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "[bitCooccurrence] Time elapsed: " << elapsed.count() << " seconds" << std::endl;

    // Convert to Python dictionary
    py::dict cooccurrences_py;
    for (auto& [k, v] : global_map) {
        cooccurrences_py[py::int_(k)] = py::cast(v);
    }

    return cooccurrences_py;
}



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
py::array_t<double> BitMatrix::columnEntropy() const {
    py::array_t<double> entropies(_n_cols);
    auto entropies_ptr = entropies.mutable_data();

    for (size_t j = 0; j < _n_cols; ++j) {
        if (_n_rows == 0) {
            entropies_ptr[j] = 0.0;
            continue;
        }
        double p = static_cast<double>(_col_masses[j]) / _n_rows;

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
py::array_t<double> BitMatrix::klDivergence(const std::vector<size_t>& ingroup_indices) const {
    const size_t ingroup_size = ingroup_indices.size();
    if (ingroup_size == 0 || ingroup_size == _n_rows) {
        // if ingroup is empty or all rows then discrimination is meaningless
        py::array_t<double> scores(_n_cols);
        std::fill(scores.mutable_data(), scores.mutable_data() + _n_cols, 0.0);
        return scores;
    }
    const size_t outgroup_size = _n_rows - ingroup_size;

    // calculate column masses for the ingroup
    // could probably be optimised at some point
    std::vector<int> ingroup_col_masses(_n_cols, 0);
    for (const auto& row_idx : ingroup_indices) {
        if (row_idx >= _n_rows) continue;   // safety check
        for (size_t k = 0; k < _packed_cols; ++k) {
            uint64_t chunk = _packed_matrix[row_idx][k];
            while (chunk != 0) {
                int trailing_zeros = __builtin_ctzll(chunk);
                size_t col_idx = k * 64 + trailing_zeros;
                if (col_idx < _n_cols) {
                    ingroup_col_masses[col_idx]++;
                }
                chunk &= chunk - 1;   // clear the LSB
            }
        }
    }

    py::array_t<double> scores(_n_cols);
    auto scores_ptr = scores.mutable_data();

    for (size_t j = 0; j < _n_cols; ++j) {
        double ingroup_mass = static_cast<double>(ingroup_col_masses[j]);
        double outgroup_mass = static_cast<double>(_col_masses[j] - ingroup_col_masses[j]);

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
 * ardal::BitMatrix::getDensity
 *
 * Get the density of the matrix.
 *
 * This function returns the pre-calculated density of the matrix. The density is defined as the
 * proportion of set bits (1s) in the entire matrix.
 *
 * INPUT: 
 *  None (operates on the private member _density)
 *
 * OUTPUT:
 *  double : The density of the matrix.
 ****************************************************************************************************/
 
double BitMatrix::getDensity( void ) const {
    return _density;
}


/****************************************************************************************************
 * ardal::BitMatrix::getRowMasses
 *
 * Get the popcount of each row in the matrix.
 *
 * This function returns the pre-calculated popcount of each row in the matrix. The popcount of a row
 * represents the number of set bits (1s) present in that row.
 *
 * INPUT: 
 *  None (operates on the private member _row_masses)
 *
 * OUTPUT:
 *  std::vector<int> : A vector containing the popcount of each row in the matrix.
 ****************************************************************************************************/
py::array_t<int> BitMatrix::getRowMasses( void ) {
    return py::array_t<int>(_row_masses.size(), _row_masses.data());
}


const std::vector<int>& BitMatrix::getRowMassesVector( void ) {
    return _row_masses;
}


/****************************************************************************************************
 * ardal::BitMatrix::getColumnMasses
 *
 * Get the popcount of each column in the matrix.
 *
 * This function returns the pre-calculated popcount of each column in the matrix. The popcount of a
 * columnn represents the number of set bits (1s) present in that column.
 *
 * INPUT: 
 *  None (operates on the private member _col_masses)
 *
 * OUTPUT:
 *  std::vector<int> : A vector containing the popcount of each column in the matrix.
 ****************************************************************************************************/
py::array_t<int> BitMatrix::getColumnMasses( void ) {
    return py::array_t<int>(_col_masses.size(), _col_masses.data());
}

} // namespace _ardal