// BitMatrix.hpp
/*
Copyright 2025 Arthur V. Morris
*/
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

    // distance functions
    py::array_t<uint32_t> hamming( bool fill_cache = false, bool use_simd = true, int threads = 1 ) const;
    py::array_t<int> innerProduct( bool fill_cache = false, bool use_simd = true, int threads = 1 ) const;
 
    // neighbourhood functions
    py::array_t<int64_t> neighbourhood( size_t row_idx, int epsilon, bool use_simd = true, int threads = 1 ) const;
    py::list innerProductNeighbourhood( size_t row_idx, int ip_epsilon, bool use_simd = true ) const;

    // set operation functions
    std::vector<size_t> uniqueSharedBits( const std::vector<size_t>& row_indices, bool use_simd = true ) const;

    // statistical functions
    py::dict bitCooccurrence( double threshold = 0.95, int threads = 1 ) const;
    py::array_t<double> colFrequency( std::vector<size_t>& row_indices ) const;
    py::array_t<double> columnEntropy( void ) const;
    py::array_t<double> klDivergence( const std::vector<size_t>& ingroup_indices ) const;
    py::array_t<double> jsDivergence( const std::vector<size_t>& ingroup_indices ) const;
    py::array_t<double> informationGain( const std::vector<size_t>& ingroup_indices ) const;

    // get functions
    py::array_t<size_t> getSetBitIndices( size_t row_idx ) const;
    py::array_t<int> getRowMasses( void );
    const std::vector<int>& getRowMassesVector( void );
    py::array_t<int> getColumnMasses( void );
    double getDensity( void ) const;
    py::array_t<uint8_t> getBitMatrix( void ) const;

private:
    // bit-packed matrix
    std::vector<std::vector<uint64_t>> _packed_matrix;

    // attributes
    size_t _n_rows;                 // the number of rows (guids)
    size_t _n_cols;                 // the number of columns (alleles)
    size_t _packed_cols;            // the number of packed columns
    std::vector<int> _row_masses;   // the mass of each row
    std::vector<int> _col_masses;   // the mass of each column
    double _density;                // the density of the matrix

    // access functions
    std::vector<size_t> getRowSetBitIndices( size_t row_idx ) const;
    std::vector<uint64_t> getColumnVector(size_t col_idx) const;
    std::vector<size_t> getColSetBitIndices( size_t col_idx ) const;

    // distance functions
    uint32_t hammingDistance_scalar( size_t i, size_t j ) const;
    uint32_t hammingDistance_SIMD( size_t i, size_t j ) const;

    // neighbourhood functions
    int epsilonNeighbourhood_scalar( size_t i, size_t j, int epsilon ) const;
    int epsilonNeighbourhood_SIMD( size_t i, size_t j, int epsilon ) const;
    int innerProduct_scalar( size_t i, size_t j ) const;
    int innerProduct_SIMD( size_t i, size_t j ) const;

    // statistics helper functions
    double density( void ) const;

    // hamming distance cache
    std::unordered_map<std::pair<size_t, size_t>, int, pair_hash> _hamming_cache;

    // internal helpers
    // bit unpacking
    inline bool getBit( size_t row, size_t col ) const {
        uint64_t byte = _packed_matrix[row][col / 64];
        return (byte >> (col % 64)) & 1;
    }
};

} // namespace _ardal