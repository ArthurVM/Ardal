// AlleleMatrix.hpp (Adjusted to Original Code Structure)
#pragma once

#include <vector>
#include <unordered_map>
#include <utility>
#include <cstdint>
#include <set>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <immintrin.h>

namespace py = pybind11;
namespace _ardal {

// Custom hash function for std::pair
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};


class AlleleMatrix {
public:
    // constructor: takes a NumPy matrix
    AlleleMatrix( py::array_t<uint8_t> input_matrix );

    // distance methods
    py::array_t<int> hamming( bool fill_cache = false, bool use_simd = true ) const;
 
    // neighbourhood methods
    py::list neighbourhood( size_t row_idx, int epsilon, bool use_simd = true )  const;
    py::list innerProductNeighbourhood( size_t row_idx, int ip_epsilon, bool use_simd = true ) const;

    // get methods
    py::array_t<size_t> getAlleles( size_t row_idx ) const;
    py::array_t<int> getMass( void );
    py::array_t<uint8_t> getMatrix( void ) const;

private:
    // bit-packed matrix
    std::vector<std::vector<uint8_t>> _packed_matrix;

    // attributes
    size_t _n_rows;            // the number of rows (guids)
    size_t _n_cols;            // the number of columns (alleles)
    size_t _packed_cols;       // the number of packed columns
    std::vector<int> _rmass;   // the mass of each row

    // access methods
    std::vector<size_t> accessGUID( size_t row_idx ) const;

    // distance methods
    int hammingDistanceScalar( size_t i, size_t j ) const;
    int hammingDistanceSIMD( size_t i, size_t j ) const;

    // neighbourhood methods
    int epsilonNeighbourhoodScalar( size_t i, size_t j, int epsilon ) const;
    int epsilonNeighbourhoodSIMD( size_t i, size_t j, int epsilon ) const;
    int innerProductScalar( size_t i, size_t j ) const;
    int innerProductSIMD( size_t i, size_t j ) const;

    // row mass methods
    std::vector<int> mass( void ) const;
    int rowMass( size_t row_idx ) const;

    // hamming distance cache
    std::unordered_map<std::pair<size_t, size_t>, int, pair_hash> _hamming_cache;

    // internal helper: bit unpacking
    inline bool get_allele( size_t row, size_t col ) const {
        uint8_t byte = _packed_matrix[row][col / 8];
        return (byte >> (col % 8)) & 1;
    }
};

} // namespace _ardal