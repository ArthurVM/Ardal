// BitMatrix.hpp
/*
Copyright 2025 Arthur V. Morris
*/
#pragma once

#include <algorithm>
#include <atomic>
#include <vector>
#include <unordered_map>
#include <queue>
#include <limits>
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
#include "detail/flat_matrix.hpp"


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
    BitMatrix( ardal::detail::FlatMatrix matrix_,
               std::shared_ptr<const std::vector<int>> row_masses,
               std::shared_ptr<const std::vector<int>> col_masses,
               const std::vector<std::vector<uint32_t>>* missing_mask = nullptr );

    // distance functions
    py::array_t<uint32_t> hamming( bool fill_cache = false,
                                   bool use_simd = true,
                                   int threads = 1,
                                   bool mask_missing = false ) const;
    py::array_t<uint32_t> hamming_subset( const std::vector<size_t>& row_indices,
                                          const std::vector<size_t>& col_indices,
                                          bool use_simd,
                                          int threads,
                                          bool mask_missing = false ) const;
    py::array_t<int> innerProduct( bool fill_cache = false,
                                   bool use_simd = true,
                                   int threads = 1,
                                   bool mask_missing = false ) const;
    py::array_t<double> cosineDistanceAll( bool use_simd = true,
                                           int threads = 1,
                                           bool mask_missing = false ) const;
    py::array_t<double> cosineDistance_subset( const std::vector<size_t>& row_indices,
                                               const std::vector<size_t>& col_indices,
                                               bool use_simd = true,
                                               int threads = 1,
                                               bool mask_missing = false ) const;
 
    // neighbourhood functions
    py::array_t<int64_t> neighbourhood( size_t row_idx,
                                        int epsilon,
                                        bool use_simd = true,
                                        int threads = 1 ) const;
    py::list innerProductNeighbourhood( size_t row_idx,
                                        int ip_epsilon,
                                        bool use_simd = true ) const;
    py::list knn_hamming( size_t row_idx,
                          int k,
                          bool use_simd = true,
                          int threads = 1 ) const;
    py::list knn_inner_product( size_t row_idx,
                                int k,
                                bool use_simd = true,
                                int threads = 1 ) const;
    py::list knn_jaccard( size_t row_idx,
                          int k,
                          bool use_simd = true,
                          int threads = 1 ) const;
    py::list knn_cosine( size_t row_idx,
                         int k,
                         bool use_simd = true,
                         int threads = 1 ) const;


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
    py::array_t<uint64_t> getSubsetMatrix_rows( const std::vector<size_t>& row_indices, 
                                                const int threads,
                                                bool silence_log = false ) const;
    py::array_t<uint64_t> getPackedMatrix( void ) const;

private:
    // bit-packed matrix
    ardal::detail::FlatMatrix matrix_;                     // packed matrix
    std::vector<const std::uint64_t*> row_ptrs_;           // row pointers
    using SparseMask = std::vector<std::pair<size_t, uint64_t>>;
    std::vector<SparseMask> row_masks_sparse_;
    std::vector<uint32_t> row_missing_ones_;
    bool has_missing_mask_ = false;

    // attributes
    const uint64_t* base_;                                 // pointer to the matrix base
    const size_t n_rows_;                                  // the number of rows (guids)
    const size_t n_cols_bits_;                             // the number of columns (alleles)
    const size_t wpr_;                                     // the number of packed columns
    std::shared_ptr<const std::vector<int>> row_masses_;   // the mass of each row
    std::shared_ptr<const std::vector<int>> col_masses_;   // the mass of each column
    double density_;                                       // the density of the matrix
    uint64_t tail_mask_ = ~0ULL;

    // access functions
    std::vector<size_t> getRowSetBitIndices( size_t row_idx ) const;
    // std::vector<uint64_t> getColumnVector(size_t col_idx) const;
    std::vector<size_t> getColSetBitIndices( size_t col_idx ) const;
    std::vector<std::vector<uint64_t>> subsetPackedMatrix( const std::vector<size_t>& row_indices,
                                                           const std::vector<size_t>& col_indices,
                                                           const int threads = 1 ) const;
    std::vector<uint64_t> buildColumnMask( const std::vector<size_t>& col_indices ) const;
    bool isFullColumnSelection( const std::vector<size_t>& col_indices ) const;
    uint32_t maskedHamming( const uint64_t* lhs,
                            const uint64_t* rhs,
                            const std::vector<uint64_t>& mask ) const;
    uint32_t maskedInnerProduct( const uint64_t* lhs,
                                 const uint64_t* rhs,
                                 const std::vector<uint64_t>& mask ) const;
    uint32_t maskedInnerProductPair( const uint64_t* lhs,
                                     const uint64_t* rhs,
                                     const SparseMask& lhs_mask,
                                     const SparseMask& rhs_mask,
                                     const std::vector<uint64_t>* column_mask = nullptr ) const;
    uint32_t maskedInnerProductPenaltySparse( const uint64_t* lhs,
                                              const uint64_t* rhs,
                                              const SparseMask& lhs_mask,
                                              const SparseMask& rhs_mask,
                                              const std::vector<uint64_t>* column_mask = nullptr ) const;
    uint32_t maskedRowMass( const uint64_t* row_ptr,
                            const std::vector<uint64_t>& mask ) const;
    uint32_t maskedRowMass( const uint64_t* row_ptr,
                            const SparseMask& row_mask,
                            const std::vector<uint64_t>* column_mask = nullptr ) const;
    uint32_t maskedHammingPair( const uint64_t* lhs,
                                const uint64_t* rhs,
                                const SparseMask& lhs_mask,
                                const SparseMask& rhs_mask,
                                const std::vector<uint64_t>* column_mask = nullptr ) const;
    uint32_t maskedHammingPenaltySparse( const uint64_t* lhs,
                                         const uint64_t* rhs,
                                         const SparseMask& lhs_mask,
                                         const SparseMask& rhs_mask,
                                         const std::vector<uint64_t>* column_mask = nullptr ) const;

    // statistics helper functions
    double density( void ) const;

    // hamming distance cache
    std::unordered_map<std::pair<size_t, size_t>, int, pair_hash> _hamming_cache;

    // internal helpers
    // bit unpacking
    inline bool getBit(std::size_t row, std::size_t col) const {
        if (row >= n_rows_ || col >= n_cols_bits_) return false;
        const std::size_t w = col >> 6;
        const unsigned b = static_cast<unsigned>(col & 63);
        uint64_t word = row_ptrs_[row][w];
        if (w + 1 == wpr_) word &= tail_mask_;
        return (word >> b) & 1ULL;
    }

    inline uint64_t getWord(std::size_t row, std::size_t w) const {
        uint64_t x = row_ptrs_[row][w];
        if (w + 1 == wpr_) x &= tail_mask_;
        return x;
    }
};

} // namespace ardal
