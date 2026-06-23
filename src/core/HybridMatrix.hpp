// HybridMatrix.hpp
/*
Copyright 2025 Arthur V. Morris
*/
#pragma once

#include <memory>
#include <iostream>
#include <cstdint>
#include <vector>
#include "BitMatrix.hpp"
#include "RoaringMatrix.hpp"
#include "MissingRanges.hpp"
#include "detail/flat_matrix.hpp"


namespace py = pybind11;
namespace ardal {

class SnvFallbackBackend;

class HybridMatrix {
public:
    HybridMatrix( py::array packed_or_dense_matrix,
                  bool is_bitpacked,
                  size_t n_cols_bits,
                  bool use_roaring_if_sparse = true,
                  double density_threshold = 0.02,
                  py::object missing_mask = py::object() );
    ~HybridMatrix();

    // distance functions: BitMatrix and RoaringMatrix
    py::array_t<uint32_t> hamming( bool use_simd = true,
                                   int threads = 1,
                                   const std::string& backend = "auto",
                                   bool mask_missing = false ) const;
    
    py::array_t<uint32_t> hamming_subset( const std::vector<size_t> row_indices,
                                          const std::vector<size_t> col_indices,
                                          bool use_simd = true,
                                          int threads = 1,
                                          const std::string& backend = "auto",
                                          bool mask_missing = false ) const;
    
    py::array_t<double> jaccard( bool use_simd = true,
                                 int threads = 1,
                                 const std::string& backend = "auto" ) const;
                                   
    py::array_t<int> innerProduct( bool use_simd = true,
                                   int threads = 1,
                                   const std::string& backend = "auto",
                                   bool mask_missing = false ) const;
    
    py::array_t<double> cosineDistance( bool use_simd = true,
                                        int threads = 1,
                                        const std::string& backend = "auto",
                                        bool mask_missing = false ) const;
    py::array_t<double> cosineDistance_subset( const std::vector<size_t> row_indices,
                                               const std::vector<size_t> col_indices,
                                               bool use_simd = true,
                                               int threads = 1,
                                               const std::string& backend = "auto",
                                               bool mask_missing = false ) const;
 

    // neighbourhood functions: BitMatrix and RoaringMatrix
    py::array_t<int64_t> neighbourhood( size_t row_idx,
                                        uint32_t epsilon,
                                        bool use_simd = true,
                                        int threads = 1,
                                        const std::string& backend = "auto" )  const;

    py::list innerProductNeighbourhood( size_t row_idx,
                                        uint32_t ip_epsilon,
                                        bool use_simd = true,
                                        const std::string& backend = "auto" ) const;

    py::list knn_hamming( size_t row,
                          uint32_t k,
                          bool use_simd = true,
                          int threads = 1,
                          const std::string& backend = "auto" ) const;

    py::list knn_inner_product( size_t row,
                                uint32_t k,
                                bool use_simd = true,
                                int threads = 1,
                                const std::string& backend = "auto" ) const;

    py::list knn_jaccard( size_t row,
                          uint32_t k,
                          bool use_simd = true,
                          int threads = 1,
                          const std::string& backend = "auto" ) const;

    py::list knn_cosine( size_t row,
                         uint32_t k,
                         bool use_simd = true,
                         int threads = 1,
                         const std::string& backend = "auto" ) const;

    py::list knn( size_t row,
                  uint32_t k,
                  const std::string& metric,
                  bool use_simd = true,
                  int threads = 1,
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
    py::array_t<int> getRowMasses( void ) const;

    py::array_t<int> getColumnMasses( void ) const;

    double getDensity( void ) const;

    py::array_t<uint8_t> getBitMatrix( void ) const;

    py::array_t<uint64_t> getPackedMatrix( void ) const;

    py::list getRoaringMatrix( void ) const;

    bool roaringEnabled( void ) const;

    py::array_t<size_t> getSubsetPackedMatrix( const std::vector<size_t>& row_indices, 
                                               const std::vector<size_t>& col_indices,
                                               const int threads ) const;
    
    py::array_t<uint64_t> getSubsetPackedMatrix_rows( const std::vector<size_t>& row_indices, 
                                                      const int threads,
                                                      bool silence_log = false ) const;

    void prepareSnvView( py::array_t<uint32_t> allele_to_locus,
                         py::array_t<uint8_t> allele_to_base );

    py::array_t<uint32_t> snvHamming( int threads = 1,
                                      bool mask_missing = false ) const;
    py::array_t<uint32_t> snvHamming_subset( const std::vector<size_t>& row_indices,
                                            const std::vector<size_t>& col_indices,
                                            int threads = 1,
                                            bool mask_missing = false ) const;
    py::array_t<int64_t> snvNeighbourhood( size_t row,
                                           uint32_t epsilon,
                                           int threads = 1,
                                           bool mask_missing = false ) const;
    py::list knnSnv( size_t row,
                     uint32_t k,
                     int threads = 1,
                     bool mask_missing = false ) const;

private:  
    // attributes
    size_t n_rows_ = 0;
    size_t n_cols_bits_ = 0;
    size_t words_per_row_ = 0;
    double density_;
    py::array owner_;
    ardal::detail::FlatMatrix flat_;
    std::vector<uint64_t> storage_;
    bool has_missing_mask_ = false;
    std::vector<std::vector<uint32_t>> mask_columns_;
    MissingRanges missing_ranges_;
    bool has_missing_ranges_ = false;
    
    // backends
    std::unique_ptr<BitMatrix> bit_backend;
    std::unique_ptr<RoaringMatrix> roaring_backend;
    bool roaring_enabled = false;

    // backend selector
    enum class BackendType { AUTO, BIT, ROARING };
    BackendType selectBackend(size_t n_cols, double density) const;

    mutable std::unique_ptr<SnvFallbackBackend> snv_fallback_backend_;

    inline BackendType parse_backend(const std::string& backend) const {
        if (backend == "auto" || backend == "AUTO") return BackendType::AUTO;
        if (backend == "bit"  || backend == "BIT")  return BackendType::BIT;
        if (backend == "roaring" || backend == "ROARING") return BackendType::ROARING;
        throw std::invalid_argument("Unknown backend: " + backend);
    }
};

} // namespace ardal
