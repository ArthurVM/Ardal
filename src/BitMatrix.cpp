/*
Copyright 2025 Arthur V. Morris
*/
#include "BitMatrix.hpp"
#include <stdexcept>
#include <cmath>
#include <cstring>

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
    _packed_cols = (_n_cols + 7) / 8;

    // size check
    if (_n_rows > std::numeric_limits<size_t>::max() / _n_cols) {
        throw std::runtime_error("Matrix dimensions are too large, potential overflow.");
    }

    // allocate memory for packed matrix
    auto ptr = static_cast<uint8_t*>(buf.ptr);
    _packed_matrix.resize(_n_rows, std::vector<uint8_t>(_packed_cols, 0));
    
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
                _packed_matrix[i][j / 8] |= (1 << (j % 8));
                // populate mass counters
                _row_masses[i]++;
                _col_masses[j]++;
            }
        }
    }
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
    for (size_t k = 0; k < _n_cols; ++k) {
        if (getBit(row_idx, k)) {
            row_indices.push_back(k);
        }
    }
    return row_indices;
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
 * ardal::BitMatrix::getMatrix
 * 
 * Unpack and return the bit matrix.
 * 
 * INPUT: 
 *  None (operates on the private member _packed_matrix)
 *
 * OUTPUT:
 *  py::array_t<uint8_t> : A 2D numpy array representing a binary matrix.
 ****************************************************************************************************/
py::array_t<uint8_t> BitMatrix::getMatrix( void ) const {
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
int BitMatrix::hammingDistanceScalar( size_t i, size_t j ) const {
    // input validation in the name of paranoia/robustness
    if (i >= _n_rows || j >= _n_rows) {
        throw std::out_of_range("Row index out of bounds in hammingDistanceScalar.");
    }

    int dist = 0;
    for (size_t k = 0; k < _packed_cols; ++k) {
        uint8_t xor_byte = _packed_matrix[i][k] ^ _packed_matrix[j][k];
        dist += __builtin_popcount(xor_byte);
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
int BitMatrix::hammingDistanceSIMD( size_t i, size_t j ) const {
    // input valiation for robustness
    if (i >= _n_rows || j >= _n_rows) {
        throw std::out_of_range("Row index out of bounds in hammingDistanceSIMD.");
    }

    int dist = 0;
    size_t k = 0;

    // SIMD block
    for (; k + 31 < _packed_cols; k += 32) {
        __m256i a = _mm256_loadu_si256((__m256i const*)&_packed_matrix[i][k]);
        __m256i b = _mm256_loadu_si256((__m256i const*)&_packed_matrix[j][k]);
        __m256i xor_result = _mm256_xor_si256(a, b);

        // horizontal popcount using 8 lanes of 32-bit integers
        for (int lane = 0; lane < 8; ++lane) {
            uint32_t chunk = _mm256_extract_epi32(xor_result, lane);
            dist += __builtin_popcount(chunk);
        }
    }

    // cover the non 32 chunk tail
    for (; k < _packed_cols; ++k) {
        uint8_t xor_byte = _packed_matrix[i][k] ^ _packed_matrix[j][k];
        dist += __builtin_popcount(xor_byte);
    }
    return dist;
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
 *
 * OUTPUT:
 *  py::array_t<int>  : A 1D NumPy array representing the condensed distance matrix containing
 *                      the pairwise Hamming distances.  The length of the array is n*(n-1)/2,
 *                      where 'n' is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<int> BitMatrix::hamming( bool fill_cache, bool use_simd ) const {
    size_t total_pairs = _n_rows * (_n_rows - 1) / 2;
    py::array_t<int> dist_matrix(total_pairs);
    auto dist_ptr = dist_matrix.mutable_data();

    size_t idx = 0;
    for (size_t i = 0; i < _n_rows; ++i) {
        for (size_t j = i + 1; j < _n_rows; ++j) {
            // run as SIMD
            if (use_simd)
                dist_ptr[idx++] = hammingDistanceSIMD(i, j);
            // run as scalar
            else
                dist_ptr[idx++] = hammingDistanceScalar(i, j);
        }
    }
    return dist_matrix;
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
    for (; k + 31 < _packed_cols; k += 32) {
        __m256i a = _mm256_loadu_si256((__m256i const*)&_packed_matrix[i][k]);
        __m256i b = _mm256_loadu_si256((__m256i const*)&_packed_matrix[j][k]);
        __m256i xor_result = _mm256_xor_si256(a, b);

        // horizontal popcount using 8 lanes of 32-bit integers
        for (int lane = 0; lane < 8; ++lane) {
            uint32_t chunk = _mm256_extract_epi32(xor_result, lane);
            dist += __builtin_popcount(chunk);
        }
        if (dist > epsilon) {
            return dist;   // early break where epsilon is exceeded
        }
    }

    // cover the non 32 chunk tail
    for (; k < _packed_cols; ++k) {
        uint8_t xor_byte = _packed_matrix[i][k] ^ _packed_matrix[j][k];
        dist += __builtin_popcount(xor_byte);
    }
    return dist;
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
        uint8_t xor_byte = _packed_matrix[i][k] ^ _packed_matrix[j][k];
        dist += __builtin_popcount(xor_byte);
        if (dist > epsilon) {
            return dist;  // early break where epsilon is exceeded
        }
    }
    return dist;
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
 * OUTPUT:
 *  py::list<tuple(row, dist)> : A 1D NumPy array containing the indices of the rows that are within 
 *                               the epsilon-neighborhood of the target row.
 *
 * EXCEPTIONS:
 *  std::runtime_error : If row_coord is out of range.
 ****************************************************************************************************/
py::list BitMatrix::neighbourhood( size_t row_idx, int epsilon, bool use_simd ) const {
    // do some data cleanliness
    if (epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_idx < 0) {
        throw std::runtime_error("Row index must be non-negative.");
        }
    if (row_idx >= _n_rows) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }

    py::list ep_n;   // store results
    int q_mass = _row_masses[row_idx];   // access pre-calculated row mass

    for (size_t i = 0; i < _n_rows; ++i) {
        if (i != row_idx) {
            // row mass filter
            int mass_d = std::abs(q_mass - _row_masses[i]);
            if (mass_d > epsilon) {
                continue;   // skip neighbourhood evaluation
            }

            int distance;   // initialise distance
            if (use_simd) {
                distance = epsilonNeighbourhoodSIMD(row_idx, i, epsilon);
            } else {
                distance = epsilonNeighbourhoodScalar(row_idx, i, epsilon);
            }
            if (distance <= epsilon) {
                ep_n.append(py::make_tuple(i, distance));
            }
        }
    }
    return ep_n;
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
        uint8_t and_byte = _packed_matrix[i][k] & _packed_matrix[j][k];
        ip += __builtin_popcount(and_byte);
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
    for (; k + 31 < _packed_cols; k += 32) {
        __m256i a = _mm256_loadu_si256((__m256i const*)&_packed_matrix[i][k]);
        __m256i b = _mm256_loadu_si256((__m256i const*)&_packed_matrix[j][k]);
        __m256i and_result = _mm256_and_si256(a, b);

        // horizontal popcount using 8 lanes of 32-bit integers
        for (int lane = 0; lane < 8; ++lane) {
            uint32_t chunk = _mm256_extract_epi32(and_result, lane);
            ip += __builtin_popcount(chunk);
        }
    }

    // cover the non 32 chunk tail
    for (; k < _packed_cols; ++k) {
        uint8_t and_byte = _packed_matrix[i][k] & _packed_matrix[j][k];
        ip += __builtin_popcount(and_byte);
    }
    return ip;
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
    if (row_idx >= _n_rows) {
        throw std::runtime_error("Query row index out of range.");
    }
    if (ip_epsilon < 0) {
        throw std::runtime_error("ip_epsilon must be non-negative.");
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


// Pybind methods
PYBIND11_MODULE(_ardal, m) {  // _ardal module and method bindings
    py::class_<_ardal::BitMatrix>(m, "BitMatrix")
        .def(py::init<py::array_t<uint8_t>&>())
        .def("hamming", &_ardal::BitMatrix::hamming, py::arg("fill_cache") = false, py::arg("use_simd") = true)
        .def("neighbourhood", &_ardal::BitMatrix::neighbourhood, py::arg("row_idx"), py::arg("epsilon"), py::arg("use_simd") = true)
        .def("innerProductNeighbourhood", &_ardal::BitMatrix::innerProductNeighbourhood, py::arg("row_idx"), py::arg("ip_epsilon"), py::arg("use_simd") = true)
        .def("getSetBitIndices", &_ardal::BitMatrix::getSetBitIndices, py::arg("row_idx"))
        .def("getRowMasses", &_ardal::BitMatrix::getRowMasses)
        .def("getColumnMasses", &_ardal::BitMatrix::getColumnMasses)
        .def("getMatrix", &_ardal::BitMatrix::getMatrix);
}