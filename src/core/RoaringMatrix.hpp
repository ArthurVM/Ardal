// RoaringMatrix.hpp
/*
Copyright 2025 Arthur V. Morris
*/
#pragma once

#include <vector>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <cstdint>
#include <set>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <immintrin.h>
#include <omp.h>
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <sstream>
#include "core/types.hpp"


namespace py = pybind11;
namespace ardal {

class RoaringMatrix {
public:
    // constructor: takes a NumPy matrix
    RoaringMatrix( std::shared_ptr<const ardal::detail::WordsVV> vv,
                   std::shared_ptr<const std::vector<int>> row_masses,
                   std::shared_ptr<const std::vector<int>> col_masses,
                   std::size_t n_rows,
                   std::size_t n_cols_bits );

    // distance functions
    py::array_t<uint32_t> hamming( int threads = 1 ) const;
    py::array_t<uint32_t> hamming_subset( const std::vector<size_t>& row_indices,
                                          const std::vector<size_t>& col_indices,
                                          int threads = 1 ) const;
    py::array_t<double> jaccard( int threads = 1 ) const;
    py::array_t<int> innerProduct( int threads = 1 ) const;
 
    // neighbourhood functions
    py::array_t<int64_t> neighbourhood( size_t row_idx, int epsilon, int threads = 1 ) const;
    py::list innerProductNeighbourhood( size_t row_idx, int ip_epsilon, int threads = 1 ) const;

    // set operation functions
    std::vector<size_t> uniqueSharedBits( const std::vector<size_t>& row_indices ) const;

    // get functions
    py::array_t<size_t> getSetBitIndices( size_t row_idx ) const;
    py::list getRoaringMatrix( void ) const;

    // colwise functions
    py::dict bitCooccurrence_all( double threshold, int threads = 1 ) const;
    py::dict bitCooccurrence_subset( const std::vector<size_t>& col_indices, double threshold, int threads ) const;


private:
    // roaring matrix
    std::vector<roaring::Roaring> roaring_bitmap_;

    // attributes
    std::shared_ptr<const std::vector<int>> row_masses_;
    std::shared_ptr<const std::vector<int>> col_masses_;
    size_t n_rows_;
    size_t n_cols_;
    size_t _n_packed_cols;

    // distance functions
    uint32_t hammingDistance( size_t i, size_t j ) const;
    double jaccardIndex( size_t i, size_t j ) const;
    uint32_t innerProductRowwise( size_t i, size_t j ) const;

    // utility functions
    std::vector<roaring::Roaring> colwiseRoaringFromRowwise( void ) const;
};

} // namespace ardal