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
    HybridMatrix( py::array_t<uint8_t> matrix,
                  bool use_roaring_if_sparse = true,
                  double density_threshold = 0.02 );

    // distance functions: BitMatrix and RoaringMatrix
    py::array_t<uint32_t> hamming( bool use_simd = true,
                                   int threads = 1,
                                   const std::string& backend = "auto" ) const;
    
    py::array_t<uint32_t> hamming_subset( const std::vector<size_t> row_indices,
                                          const std::vector<size_t> col_indices,
                                          bool use_simd = true,
                                          int threads = 1,
                                          const std::string& backend = "auto" ) const;
    
    py::array_t<double> jaccard( bool use_simd = true,
                                 int threads = 1,
                                 const std::string& backend = "auto" ) const;
                                   
    py::array_t<int> innerProduct( bool use_simd = true,
                                   int threads = 1,
                                   const std::string& backend = "auto" ) const;
 

    // neighbourhood functions: BitMatrix and RoaringMatrix
    py::array_t<int64_t> neighbourhood( size_t row_idx,
                                        int epsilon,
                                        bool use_simd = true,
                                        int threads = 1,
                                        const std::string& backend = "auto" )  const;

    py::list innerProductNeighbourhood( size_t row_idx,
                                        int ip_epsilon,
                                        bool use_simd = true,
                                        const std::string& backend = "auto" ) const;


    // set operation functions: BitMatrix
    std::vector<size_t> uniqueSharedBits( const std::vector<size_t>& row_indices,
                                          bool use_simd = true,
                                          const std::string& backend = "auto" ) const;


    // statistical functions: BitMatrix
    py::array_t<double> colFrequency( std::vector<size_t>& row_indices ) const;

    py::array_t<double> columnEntropy( void ) const;

    py::dict bitCooccurrence_all( double threshold = 0.95,
                                  int threads = 1 ) const;

    py::dict bitCooccurrence_subset( const std::vector<size_t>& col_indices,
                                     double threshold,
                                     int threads ) const;

    py::array_t<double> klDivergence( const std::vector<size_t>& ingroup_indices ) const;

    py::array_t<double> jsDivergence( const std::vector<size_t>& ingroup_indices ) const;

    py::array_t<double> informationGain( const std::vector<size_t>& ingroup_indices ) const;


    // get functions: BitMatrix and RoaringMatrix
    py::array_t<size_t> getSetBitIndices( size_t row_idx, 
                                          const std::string& backend = "auto" ) const;
    

    // get functions: BitMatrix
    py::array_t<int> getRowMasses( void );

    py::array_t<int> getColumnMasses( void );

    double getDensity( void ) const;

    py::array_t<uint8_t> getBitMatrix( void ) const;

    py::list getRoaringMatrix( void ) const;

    bool roaringEnabled( void ) const;

    py::array_t<size_t> getSubsetPackedMatrix( const std::vector<size_t>& row_indices, 
                                               const std::vector<size_t>& col_indices,
                                               const int threads ) const;


private:
    std::unique_ptr<BitMatrix> bit_backend;
    std::unique_ptr<RoaringMatrix> roaring_backend;
    
    // attributes
    bool roaring_enabled;
    double density;
    size_t _n_cols;
    
    // backend selector
    enum class BackendType { AUTO, BIT, ROARING };
    BackendType selectBackend(size_t n_cols, double density) const;

    inline BackendType parse_backend(const std::string& backend) const {
        if (backend == "auto" || backend == "AUTO") return BackendType::AUTO;
        if (backend == "bit"  || backend == "BIT")  return BackendType::BIT;
        if (backend == "roaring" || backend == "ROARING") return BackendType::ROARING;
        throw std::invalid_argument("Unknown backend: " + backend);
    }
};

} // namespace _ardal