// HybridMatrix.hpp
/*
Copyright 2025 Arthur V. Morris
*/
#pragma once

#include <memory>
#include "BitMatrix.hpp"
#include "RoaringMatrix.hpp"

namespace py = pybind11;
namespace _ardal {

class HybridMatrix {
public:
    HybridMatrix(py::array_t<uint8_t> matrix,
                 bool use_roaring_if_sparse = true,
                 double density_threshold = 0.02);

    // distance functions: BitMatrix and RoaringMatrix
    py::array_t<int> hamming( bool use_simd = true,
                              int threads = 1,
                              bool force_bit_backend = false ) const;

    py::array_t<int> innerProduct( bool use_simd = true,
                                   int threads = 1,
                                   bool force_bit_backend = false ) const;
 

    // neighbourhood functions: BitMatrix and RoaringMatrix
    py::array_t<int64_t> neighbourhood( size_t row_idx,
                                        int epsilon,
                                        bool use_simd = true,
                                        int threads = 1,
                                        bool force_bit_backend = false )  const;

    py::list innerProductNeighbourhood( size_t row_idx,
                                        int ip_epsilon,
                                        bool use_simd = true,
                                        bool force_bit_backend = false ) const;


    // set operation functions: BitMatrix
    std::vector<size_t> uniqueSharedBits( const std::vector<size_t>& row_indices,
                                          bool use_simd = true,
                                          bool force_bit_backend = false ) const;


    // statistical functions: BitMatrix
    py::array_t<double> colFrequency( std::vector<size_t>& row_indices ) const;

    py::array_t<double> columnEntropy( void ) const;

    py::dict bitCooccurrence( double threshold = 0.95, bool use_simd = true, int cache_bytes = 1e+9 ) const;

    py::array_t<double> klDivergence( const std::vector<size_t>& ingroup_indices ) const;


    // get functions: BitMatrix and RoaringMatrix
    py::array_t<size_t> getSetBitIndices( size_t row_idx, bool force_bit_backend = false ) const;
    

    // get functions: BitMatrix
    py::array_t<int> getRowMasses( void );

    py::array_t<int> getColumnMasses( void );

    double getDensity( void ) const;

    py::array_t<uint8_t> getBitMatrix( void ) const;

    py::list getRoaringMatrix( void ) const;

    bool roaringEnabled( void ) const;


private:
    std::unique_ptr<BitMatrix> bit_backend;
    std::unique_ptr<RoaringMatrix> roaring_backend;
    
    // attributes
    bool roaring_enabled;
    double density;
};

} // namespace _ardal