// RoaringMatrix.hpp
/*
Copyright 2025 Arthur V. Morris
*/
#pragma once

#include "roaring/roaring.hh"
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


namespace py = pybind11;
namespace _ardal {

class RoaringMatrix {
public:
    // constructor: takes a NumPy matrix
    RoaringMatrix( py::array_t<uint8_t> input_matrix, const std::vector<int>& row_masses );

    // distance functions
    py::array_t<uint32_t> hamming( int threads = 1 ) const;
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
    py::dict bitCooccurrence( double threshold, int threads = 1 ) const;


private:
    // roaring matrix
    std::vector<roaring::Roaring> _roaring_matrix;

    // attributes
    const std::vector<int>& _row_masses;
    size_t _n_rows;
    size_t _n_cols;

    // distance functions
    uint32_t hammingDistance( size_t i, size_t j ) const;
    uint32_t innerProductRowwise( size_t i, size_t j ) const;

    // utility functions
    std::vector<roaring::Roaring> colwiseRoaringFromRowwise() const;

};

} // namespace _ardal