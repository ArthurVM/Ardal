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
#include <queue>
#include <limits>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <immintrin.h>
#include <omp.h>
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <sstream>
#include <algorithm>
#include "roaring/roaring.hh"
#include "detail/flat_matrix.hpp"
#include "utils/bitops.hpp"


namespace py = pybind11;
namespace ardal {

struct SnvEntry {
    uint32_t col;
    uint32_t locus;
    uint8_t base;
};

class RoaringMatrix {
public:
    // constructor: takes a NumPy matrix
    RoaringMatrix( ardal::detail::FlatMatrix flat_matrix,
                   std::shared_ptr<const std::vector<int>> row_masses,
                   std::shared_ptr<const std::vector<int>> col_masses,
                   const std::vector<std::vector<uint32_t>>* missing_mask = nullptr );

    // distance functions
    py::array_t<uint32_t> hamming( int threads = 1, bool mask_missing = false ) const;
    py::array_t<uint32_t> hamming_subset( const std::vector<size_t>& row_indices,
                                          const std::vector<size_t>& col_indices,
                                          int threads = 1,
                                          bool mask_missing = false ) const;
    py::array_t<double> jaccard( int threads = 1 ) const;
    py::array_t<int> innerProduct( int threads = 1 ) const;
    py::array_t<double> cosineDistance( int threads = 1 ) const;
    py::array_t<double> cosineDistance_subset( const std::vector<size_t>& row_indices,
                                               const std::vector<size_t>& col_indices,
                                               int threads = 1 ) const;
 
    // neighbourhood functions
    py::array_t<int64_t> neighbourhood( size_t row_idx, int epsilon, int threads = 1 ) const;
    py::list innerProductNeighbourhood( size_t row_idx, int ip_epsilon, int threads = 1 ) const;
    py::list knn_hamming( size_t row_idx, int k, int threads = 1 ) const;
    py::list knn_inner_product( size_t row_idx, int k, int threads = 1 ) const;
    py::list knn_jaccard( size_t row_idx, int k, int threads = 1 ) const;
    py::list knn_cosine( size_t row_idx, int k, int threads = 1 ) const;

    // set operation functions
    std::vector<size_t> uniqueSharedBits( const std::vector<size_t>& row_indices ) const;

    // get functions
    py::array_t<size_t> getSetBitIndices( size_t row_idx ) const;
    py::list getRoaringMatrix( void ) const;

    // SNV functions
    void prepareSnvView( py::array_t<uint32_t> allele_to_locus,
                         py::array_t<uint8_t> allele_to_base );
    py::array_t<uint32_t> snvHamming( int threads = 1,
                                      bool mask_missing = false ) const;
    py::array_t<uint32_t> snvHamming_subset( const std::vector<size_t>& row_indices,
                                             const std::vector<size_t>& col_indices,
                                             int threads = 1,
                                             bool mask_missing = false ) const;
    py::array_t<int64_t> snvNeighbourhood( size_t row_idx,
                                           uint32_t epsilon,
                                           int threads = 1,
                                           bool mask_missing = false ) const;
    py::list knnSnv( size_t row_idx,
                     int k,
                     int threads = 1,
                     bool mask_missing = false ) const;

    // colwise functions
    py::dict bitCooccurrence_all( double threshold, int threads = 1 ) const;
    py::dict bitCooccurrence_subset( const std::vector<size_t>& col_indices, double threshold, int threads ) const;


private:
    // roaring matrix
    std::vector<roaring::Roaring> roaring_bitmap_;

    // attributes
    const uint64_t* base_;                                 // pointer to the matrix base
    const size_t n_rows_;                                  // the number of rows (guids)
    const size_t n_cols_bits_;                             // the number of columns (alleles)
    const size_t wpr_;                                     // the number of packed columns
    std::shared_ptr<const std::vector<int>> row_masses_;   // the mass of each row
    std::shared_ptr<const std::vector<int>> col_masses_;   // the mass of each column
    std::vector<roaring::Roaring> missing_masks_;
    bool has_missing_mask_ = false;

    // distance functions
    uint32_t hammingDistance( size_t i, size_t j ) const;
    uint32_t maskedHammingDistance( size_t i, size_t j ) const;
    double jaccardIndex( size_t i, size_t j ) const;
    uint32_t innerProductRowwise( size_t i, size_t j ) const;

    // utility functions
    std::vector<roaring::Roaring> colwiseRoaringFromRowwise( void ) const;

    // SNV helpers
    void ensure_snv_vectors() const;
    uint32_t snvDistance(size_t i, size_t j) const;
    uint32_t snvDistanceMasked(size_t i, size_t j, const std::vector<uint8_t>* locus_mask) const;
    uint32_t snvDistanceRaw(const std::vector<uint64_t>& lhs,
                            const std::vector<uint64_t>& rhs,
                            const std::vector<uint8_t>* locus_mask) const;
    uint32_t snvDistanceMaskedEntries(const std::vector<SnvEntry>& lhs,
                                      const std::vector<SnvEntry>& rhs,
                                      const roaring::Roaring* mask_l,
                                      const roaring::Roaring* mask_r,
                                      const std::vector<uint8_t>* locus_mask) const;
    bool isFullColumnSelection( const std::vector<size_t>& col_indices ) const;
    roaring::Roaring buildColumnMask( const std::vector<size_t>& col_indices ) const;
    std::vector<uint8_t> buildLocusMask( const std::vector<size_t>& col_indices ) const;

    std::vector<uint32_t> allele_to_locus_;
    std::vector<uint8_t> allele_to_base_;
    mutable std::vector<std::vector<uint64_t>> snv_vectors_;
    mutable std::vector<std::vector<SnvEntry>> snv_entries_;
    std::vector<bool> missing_mask_empty_;
    mutable bool snv_lookup_loaded_ = false;
    mutable bool snv_vectors_ready_ = false;
    mutable bool snv_entries_ready_ = false;
};

} // namespace ardal
