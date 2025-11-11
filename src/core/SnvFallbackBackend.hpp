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
                        std::size_t n_cols_bits );

    void prepare( pybind11::array_t<uint32_t> allele_to_locus,
                  pybind11::array_t<uint8_t> allele_to_base );

    pybind11::array_t<uint32_t> snvHamming( int threads ) const;
    pybind11::array_t<uint32_t> snvHamming_subset( const std::vector<size_t>& row_indices,
                                                  const std::vector<size_t>& col_indices,
                                                  int threads ) const;
    pybind11::array_t<int64_t> snvNeighbourhood( std::size_t row_idx,
                                                 uint32_t epsilon,
                                                 int threads ) const;
    pybind11::list knnSnv( std::size_t row_idx,
                           uint32_t k,
                           int threads ) const;

private:
    void ensure_vectors() const;
    uint32_t snvDistance( std::size_t i, std::size_t j ) const;

    const ardal::detail::FlatMatrix& flat_;
    const std::size_t n_rows_;
    const std::size_t n_cols_bits_;
    const std::size_t words_per_row_;
    std::vector<uint32_t> allele_to_locus_;
    std::vector<uint8_t> allele_to_base_;
    mutable std::vector<std::vector<uint64_t>> snv_vectors_;
    bool lookup_loaded_ = false;
    mutable bool vectors_ready_ = false;
};

} // namespace ardal
