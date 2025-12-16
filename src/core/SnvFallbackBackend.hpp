// SnvFallbackBackend.hpp
/*
Copyright 2025 Arthur V. Morris
*/
#pragma once

#include <vector>
#include <cstdint>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "detail/flat_matrix.hpp"

namespace ardal {

class SnvFallbackBackend {
public:
    SnvFallbackBackend( const ardal::detail::FlatMatrix& flat,
                        std::size_t n_cols_bits,
                        const std::vector<std::vector<uint32_t>>* missing_mask = nullptr );

    void prepare( pybind11::array_t<uint32_t> allele_to_locus,
                  pybind11::array_t<uint8_t> allele_to_base );

    pybind11::array_t<uint32_t> snvHamming( int threads, bool mask_missing = false ) const;
    pybind11::array_t<uint32_t> snvHamming_subset( const std::vector<size_t>& row_indices,
                                                   const std::vector<size_t>& col_indices,
                                                   int threads,
                                                   bool mask_missing = false ) const;
    pybind11::array_t<int64_t> snvNeighbourhood( std::size_t row_idx,
                                                 uint32_t epsilon,
                                                 int threads,
                                                 bool mask_missing = false ) const;
    pybind11::list knnSnv( std::size_t row_idx,
                           uint32_t k,
                           int threads,
                           bool mask_missing = false ) const;

private:
    void ensure_vectors() const;
    void collect_encoded_masked( std::size_t row_idx,
                                 const std::vector<uint32_t>* other_mask,
                                 const std::vector<uint8_t>* locus_mask,
                                 std::vector<uint64_t>& out ) const;
    uint32_t snvDistance( std::size_t i, std::size_t j ) const;
    uint32_t snvDistanceMasked( std::size_t i,
                                std::size_t j,
                                const std::vector<uint8_t>* locus_mask ) const;
    bool isFullColumnSelection( const std::vector<size_t>& col_indices ) const;
    std::vector<uint8_t> buildLocusMask( const std::vector<size_t>& col_indices ) const;
    bool is_missing( const std::vector<uint32_t>* mask, std::size_t col ) const;
    bool locus_allowed( uint32_t locus, const std::vector<uint8_t>* locus_mask ) const;

    const ardal::detail::FlatMatrix& flat_;
    const std::size_t n_rows_;
    const std::size_t n_cols_bits_;
    const std::size_t words_per_row_;
    const std::vector<std::vector<uint32_t>>* missing_mask_;
    std::vector<uint32_t> allele_to_locus_;
    std::vector<uint8_t> allele_to_base_;
    mutable std::vector<std::vector<uint64_t>> snv_vectors_;
    struct SnvEntry {
        uint32_t col;
        uint32_t locus;
        uint8_t base;
    };
    mutable std::vector<std::vector<SnvEntry>> snv_entries_;
    bool lookup_loaded_ = false;
    mutable bool vectors_ready_ = false;
    mutable bool entries_ready_ = false;
    bool has_missing_mask_ = false;
};

} // namespace ardal
