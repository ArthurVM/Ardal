// BitMatrix.hpp
#pragma once

#include <vector>
#include <unordered_map>
#include <utility>
#include <cstdint>
#include <set>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <immintrin.h>
#include <omp.h>

namespace py = pybind11;
namespace _ardal {

// custom hash function for std::pair
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};


class BitMatrix {
public:
    // constructor: takes a NumPy matrix
    BitMatrix( py::array_t<uint8_t> input_matrix );

    // distance methods
    py::array_t<int> hamming( bool fill_cache = false, bool use_simd = true, int threads = 1 ) const;
    py::array_t<int> innerProduct( bool fill_cache = false, bool use_simd = true, int threads = 1 ) const;
 
    // neighbourhood methods
    py::array_t<int64_t> neighbourhood( size_t row_idx, int epsilon, bool use_simd = true, int threads = 1 )  const;
    py::list innerProductNeighbourhood( size_t row_idx, int ip_epsilon, bool use_simd = true ) const;

    // set operation methods
    std::vector<size_t> uniqueSharedBits(const std::vector<size_t>& row_indices, bool use_simd = true) const;

    // statistical methods
    py::array_t<double> columnEntropy() const;
    py::array_t<double> klDivergence(const std::vector<size_t>& ingroup_indices) const;

    // get methods
    py::array_t<size_t> getSetBitIndices( size_t row_idx ) const;
    py::array_t<int> getRowMasses( void );
    py::array_t<int> getColumnMasses( void );
    py::array_t<uint8_t> getMatrix( void ) const;

private:
    // bit-packed matrix
    std::vector<std::vector<uint64_t>> _packed_matrix;

    // attributes
    size_t _n_rows;            // the number of rows (guids)
    size_t _n_cols;            // the number of columns (alleles)
    size_t _packed_cols;       // the number of packed columns
    std::vector<int> _row_masses;   // the mass of each row
    std::vector<int> _col_masses;   // the mass of each column

    // access methods
    std::vector<size_t> getRowSetBitIndices( size_t row_idx ) const;

    // distance methods
    int hammingDistanceScalar( size_t i, size_t j ) const;
    int hammingDistanceSIMD( size_t i, size_t j ) const;

    // neighbourhood methods
    int epsilonNeighbourhoodScalar( size_t i, size_t j, int epsilon ) const;
    int epsilonNeighbourhoodSIMD( size_t i, size_t j, int epsilon ) const;
    int innerProductScalar( size_t i, size_t j ) const;
    int innerProductSIMD( size_t i, size_t j ) const;

    // hamming distance cache
    std::unordered_map<std::pair<size_t, size_t>, int, pair_hash> _hamming_cache;

    // internal helper: bit unpacking
    inline bool getBit( size_t row, size_t col ) const {
        uint64_t byte = _packed_matrix[row][col / 64];
        return (byte >> (col % 64)) & 1;
    }
};

} // namespace _ardal