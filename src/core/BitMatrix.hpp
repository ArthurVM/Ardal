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
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <iostream>
#include "core/types.hpp"


namespace py = pybind11;
namespace ardal {

    
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
    BitMatrix( std::shared_ptr<const ardal::detail::WordsVV> vv,
               std::shared_ptr<const std::vector<int>> row_masses,
               std::shared_ptr<const std::vector<int>> col_masses,
               size_t n_rows,
               size_t n_cols_bits );

    // distance functions
    py::array_t<uint32_t> hamming( bool fill_cache = false,
                                   bool use_simd = true,
                                   int threads = 1 ) const;
    py::array_t<uint32_t> hamming_subset( const std::vector<size_t>& row_indices,
                                          const std::vector<size_t>& col_indices,
                                          bool use_simd,
                                          int threads ) const;
    py::array_t<int> innerProduct( bool fill_cache = false,
                                   bool use_simd = true,
                                   int threads = 1 ) const;
 
    // neighbourhood functions
    py::array_t<int64_t> neighbourhood( size_t row_idx,
                                        int epsilon,
                                        bool use_simd = true,
                                        int threads = 1 ) const;
    py::list innerProductNeighbourhood( size_t row_idx,
                                        int ip_epsilon,
                                        bool use_simd = true ) const;

    // set operation functions
    std::vector<size_t> uniqueSharedBits( const std::vector<size_t>& row_indices,
                                          bool use_simd = true ) const;

    // statistical functions
    // py::dict bitCooccurrence( double threshold = 0.95,
    //                           int threads = 1 ) const;
    py::array_t<double> colFrequency( std::vector<size_t>& row_indices ) const;
    py::array_t<double> columnEntropy( void ) const;
    py::array_t<double> klDivergence( const std::vector<size_t>& ingroup_indices ) const;
    py::array_t<double> jsDivergence( const std::vector<size_t>& ingroup_indices ) const;
    py::array_t<double> informationGain( const std::vector<size_t>& ingroup_indices ) const;

    // get functions
    py::array_t<size_t> getSetBitIndices( size_t row_idx ) const;
    py::array_t<int> getRowMasses( void ) const;
    const std::vector<int>& getRowMassesVector( void ) const;
    py::array_t<int> getColumnMasses( void ) const;
    double getDensity( void ) const;
    size_t getNCols( void ) const;
    size_t getNRows( void ) const;
    py::array_t<uint8_t> getBitMatrix( void ) const;
    py::array_t<uint64_t> getSubsetMatrix( const std::vector<size_t>& row_indices, 
                                           const std::vector<size_t>& col_indices,
                                           const int threads = 1 ) const;
    py::array_t<uint64_t> getPackedMatrix( void ) const;

private:
    // bit-packed matrix
    std::shared_ptr<const ardal::detail::WordsVV> packed_matrix_;       // packed matrix
    std::vector<const std::uint64_t*> row_ptrs_;         // row pointers

    // attributes
    size_t n_rows_;                 // the number of rows (guids)
    size_t n_cols_;                 // the number of columns (alleles)
    size_t n_packed_cols_;          // the number of packed columns
    std::shared_ptr<const std::vector<int>> row_masses_;   // the mass of each row
    std::shared_ptr<const std::vector<int>> col_masses_;   // the mass of each column
    double density_;                // the density of the matrix
    uint64_t tail_mask_ = ~0ULL;

    // access functions
    std::vector<size_t> getRowSetBitIndices( size_t row_idx ) const;
    std::vector<uint64_t> getColumnVector(size_t col_idx) const;
    std::vector<size_t> getColSetBitIndices( size_t col_idx ) const;
    std::vector<std::vector<uint64_t>> subsetPackedMatrix( const std::vector<size_t>& row_indices,
                                                           const std::vector<size_t>& col_indices,
                                                           const int threads = 1 ) const;

    // statistics helper functions
    double density( void ) const;

    // hamming distance cache
    std::unordered_map<std::pair<size_t, size_t>, int, pair_hash> _hamming_cache;

    // internal helpers
    // bit unpacking
    inline bool getBit(std::size_t row, std::size_t col) const {
        if (row >= n_rows_ || col >= n_cols_) return false;
        const std::size_t w = col >> 6;
        const unsigned b = static_cast<unsigned>(col & 63);
        uint64_t word = row_ptrs_[row][w];
        if (w + 1 == n_packed_cols_) word &= tail_mask_;
        return (word >> b) & 1ULL;
    }

    inline uint64_t getWord(std::size_t row, std::size_t w) const {
        uint64_t x = row_ptrs_[row][w];
        if (w + 1 == n_packed_cols_) x &= tail_mask_;
        return x;
    }
};

} // namespace ardal