/*
Copyright 2025 Arthur V. Morris
*/
#include "RoaringMatrix.hpp"
#include <stdexcept>
#include <cmath>
#include <cstring>

namespace py = pybind11;
namespace _ardal {


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
RoaringMatrix::RoaringMatrix(py::array_t<uint8_t> input_matrix) {
    auto buf = input_matrix.request();

    if (buf.ndim != 2) {
        throw std::runtime_error("Input matrix must be 2-dimensional");
    }

    // capture matrix dimensions
    _n_rows = buf.shape[0];
    _n_cols = buf.shape[1];

    // allocate memory for roaring bitmaps
    _roaring_matrix.resize(_n_rows);
    
    // populate roaring bitmaps
    auto ptr = static_cast<uint8_t*>(buf.ptr);
    for (size_t i = 0; i < _n_rows; ++i) {
        for (size_t j = 0; j < _n_cols; ++j) {
            uint8_t val = ptr[i * _n_cols + j];
            if (val != 0 && val != 1) {
                throw std::runtime_error("Input matrix must only contain binary values (0 or 1)");
            }
            if (val) {
                _roaring_matrix[i].push_back(j);
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
 *  None (operates on the private member _roaring_matrix)
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<int> : A 1D NumPy array representing the condensed distance matrix containing
 *                     the pairwise Hamming distances. The length of the array is n*(n-1)/2,
 *                     where 'n' is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<int> RoaringMatrix::hamming( bool use_simd, int threads ) const {
    const size_t total_pairs = _n_rows * (_n_rows - 1) / 2;
    py::array_t<int> dist_matrix(total_pairs);
    auto dist_ptr = dist_matrix.mutable_data();
    
    for (size_t i = 0; i < _n_rows; ++i) {
        for (size_t j = i + 1; j < _n_rows; ++j) {
            size_t idx = (i * (2 * _n_rows - i - 1)) / 2 + (j - i - 1);
            int dist = use_simd
                    ? hammingDistanceSIMD(i, j)
                    : hammingDistanceScalar(i, j);
            dist_ptr[idx] = dist;
        }
    }
    return dist_matrix;
}


int RoaringMatrix::hammingDistanceScalar( size_t i, size_t j ) const {
    std::vector<size_t> row_i = _roaring_matrix[i];
    std::vector<size_t> row_j = _roaring_matrix[j];
    int hamming_dist = (row_i ^ row_j).cardinality();
    return hamming_dist;
}


int RoaringMatrix::hammingDistanceSIMD( size_t i, size_t j ) const {
    std::vector<size_t> row_i = _roaring_matrix[i];
    std::vector<size_t> row_j = _roaring_matrix[j];
    int hamming_dist = (row_i ^ row_j).cardinality();
    return hamming_dist;
}


py::array_t<int> RoaringMatrix::innerProduct( bool use_simd, int threads ) const {
    py::array_t<int> inner_product_matrix;
    return inner_product_matrix;
}


py::array_t<int64_t> RoaringMatrix::neighbourhood( size_t row_idx, int epsilon, bool use_simd, int threads ) const {
    py::array_t<int64_t> ep_n;
    return ep_n;
}


py::list RoaringMatrix::innerProductNeighbourhood( size_t row_idx, int ip_epsilon, bool use_simd ) const {
    py::list ip_n_list;
    return ip_n_list;
}


std::vector<size_t> RoaringMatrix::uniqueSharedBits( const std::vector<size_t>& row_indices, bool use_simd ) const {
    std::vector<size_t> shared_bits;
    return shared_bits;
}


int RoaringMatrix::epsilonNeighbourhoodScalar( size_t i, size_t j ) const {
    int dist = 0;
    return dist;
}


int RoaringMatrix::epsilonNeighbourhoodSIMD( size_t i, size_t j ) const {
    int dist = 0;
    return dist;
}


py::array_t<size_t> RoaringMatrix::getSetBitIndices( size_t row_idx ) const {
    const auto& row_vec = _roaring_matrix[row_idx];
    return py::array_t<size_t>(row_vec.size(), row_vec.data());
}


py::list RoaringMatrix::getRoaringMatrix( void ) const {
    py::list roaring_matrix;
    for (size_t i = 0; i < _n_rows; ++i) {
        const auto& row_vec = _roaring_matrix[i];
        roaring_matrix.append(py::array_t<size_t>(row_vec.size(), row_vec.data()));
    }
    return roaring_matrix;
}

} // namespace _ardal