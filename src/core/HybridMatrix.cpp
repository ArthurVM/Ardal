/*
Copyright 2025 Arthur V. Morris
*/
#include "HybridMatrix.hpp"
#include "utils/bitops.hpp"
#include "utils/PythonLogger.hpp"



namespace ardal {


/****************************************************************************************************
 * ardal::HybridMatrix::HybridMatrix
 *
 * Constructor for the HybridMatrix class.
 *
 * This constructor initializes the HybridMatrix object, which can use either a `BitMatrix`
 * or a `RoaringMatrix` as its backend, depending on the density of the input matrix.
 *
 * INPUT:
 *  packed_matrix (py::array_t<uint8_t>): A 2D 64bit packed binary matrix.
 *  use_roaring_if_sparse (bool): If true, the `RoaringMatrix` backend will be used if the
 *                                matrix density is below `density_threshold`.
 *  density_threshold (double): The density threshold below which `RoaringMatrix` is preferred.
 *
 * OUTPUT:
 *  None (constructor)
 *
 * EXCEPTIONS:
 *  Propagates exceptions from `BitMatrix` and `RoaringMatrix` constructors.
 ****************************************************************************************************/  
HybridMatrix::HybridMatrix( py::array packed_or_dense_matrix,
                            bool is_bitpacked,
                            size_t n_cols_bits,
                            bool use_roaring_if_sparse,
                            double density_threshold )
  : n_cols_bits_(n_cols_bits)
{
    const bool packed = is_bitpacked;

    if (packed) {
        // packed -> enforce u64 -> expose as a flat memmap like array with base pointer
        auto arr = pybind11::cast<pybind11::array_t<uint64_t, pybind11::array::c_style>>(packed_or_dense_matrix);
        ardal::detail::build_flat_from_packed( arr,
                                               n_cols_bits_,
                                               flat_,
                                               owner_ );
    } else {
        // dense -> pack into u64 -> expose as a flat memmap like array with base pointer
        ardal::detail::pack_dense_into_storage( packed_or_dense_matrix,
                                                n_cols_bits_,
                                                storage_,
                                                flat_ );
    }

    // metadata
    std::vector<int> row_masses, col_masses;
    int64_t total_mass = 0;
    compute_masses(flat_, row_masses, col_masses, total_mass);

    // std::stringstream ss;
    // ss << "HybridMatrix initialised: flat_.rows=" << flat_.n_rows << " flat_.wpr=" << flat_.wpr << " flat_.n_cols_bits=" << flat_.n_cols_bits;
    // ardal::utils::log_debug(static_cast<string>(ss.str()));

    density_ = static_cast<double>(total_mass) / ((static_cast<double>(flat_.n_rows) * n_cols_bits));

    // share masses
    auto rows_sp = std::make_shared<std::vector<int>>(std::move(row_masses));
    auto cols_sp = std::make_shared<std::vector<int>>(std::move(col_masses));

    // ----------- INIT BACKENDS -----------
    // handle roaring backend
    if (use_roaring_if_sparse && density_ < density_threshold) {
        roaring_enabled = true;
        ardal::utils::log_info("Initialising RoaringMatrix backend");
        roaring_backend = std::make_unique<RoaringMatrix>( flat_,
                                                           rows_sp,
                                                           cols_sp );
    } else {
        roaring_enabled = false;
        std::stringstream ss;
        ss << "Roaring not initialised: density=" << density_ << "; density_threshold=" << density_threshold;
        ardal::utils::log_debug(static_cast<string>(ss.str()));
    }

    // handle bit backend
    bit_backend = std::make_unique<BitMatrix>( flat_,
                                               rows_sp,
                                               cols_sp );
}


HybridMatrix::BackendType HybridMatrix::selectBackend( size_t n_cols,
                                                       double density ) const {
    // heuristic backend selection
    if (n_cols > 20000 && density < 0.05 && roaring_enabled) {
        return BackendType::ROARING;
    } else {
        return BackendType::BIT;
    }
}


double HybridMatrix::getDensity( void ) const {
    return density_;
}


py::array_t<uint32_t> HybridMatrix::hamming( bool use_simd,
                                             int threads,
                                             const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        ardal::utils::log_info("Computing Hamming distance matrix with RoaringMatrix backend.");
        return roaring_backend->hamming(threads);
    } else {
        ardal::utils::log_info("Computing Hamming distance matrix with BitMatrix backend.");
        return bit_backend->hamming(false, use_simd, threads);
    }
}


py::array_t<uint32_t> HybridMatrix::hamming_subset( const std::vector<size_t> row_indices,
                                                    const std::vector<size_t> col_indices,
                                                    bool use_simd,
                                                    int threads,
                                                    const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        // note: this assumes local density isnt significantly different from global
        // not a terrible assumption but may not always be the case
        chosen_backend = selectBackend(row_indices.size(), density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->hamming_subset(row_indices, col_indices, threads);
    } else {
        return bit_backend->hamming_subset(row_indices, col_indices, use_simd, threads);
    }
}


py::array_t<double> HybridMatrix::jaccard( bool use_simd,
                                           int threads,
                                           const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        // chosen_backend = selectBackend(n_cols_bits_, density_);
        chosen_backend = BackendType::ROARING;
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->jaccard(threads);
    } else {
        throw std::runtime_error("Jaccard distance is not currently supported in BitMatrix backend.");
    }
}


py::array_t<int> HybridMatrix::innerProduct( bool use_simd,
                                             int threads,
                                             const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->innerProduct(threads);
    } else {
        return bit_backend->innerProduct(false, use_simd, threads);
    }
}


py::array_t<int64_t> HybridMatrix::neighbourhood( size_t row,
                                                  int epsilon,
                                                  bool use_simd,
                                                  int threads,
                                                  const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->neighbourhood(row, epsilon, threads);
    } else {
        return bit_backend->neighbourhood(row, epsilon, use_simd, threads);
    }
}


py::list HybridMatrix::innerProductNeighbourhood( size_t row,
                                                  int ip_epsilon,
                                                  bool use_simd,
                                                  const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->innerProductNeighbourhood(row, ip_epsilon);
    } else {
        return bit_backend->innerProductNeighbourhood(row, ip_epsilon, use_simd);
    }
}


std::vector<size_t> HybridMatrix::uniqueSharedBits( const std::vector<size_t>& row_indices,
                                                    bool use_simd,
                                                    const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    // if (!force_bit_backend) {
    //     return roaring_backend->uniqueSharedBits(row_indices);
    // } else {
    //     return bit_backend->uniqueSharedBits(row_indices, use_simd);
    // }
    return bit_backend->uniqueSharedBits(row_indices, use_simd);
}


py::array_t<double> HybridMatrix::colFrequency( std::vector<size_t>& row_indices ) const {
    return bit_backend->colFrequency(row_indices);
}


py::array_t<double> HybridMatrix::columnEntropy( void ) const {
    return bit_backend->columnEntropy();
}


py::dict HybridMatrix::bitCooccurrence_all( double threshold,
                                            int threads ) const {
    
    if (roaring_enabled) {
        return roaring_backend->bitCooccurrence_all(threshold, threads);
    } else { 
        throw std::runtime_error("Bit co-occurrence is not supported in BitMatrix backend.");
    }
}


py::dict HybridMatrix::bitCooccurrence_subset( const std::vector<size_t>& col_indices,
                                               double threshold,
                                               int threads ) const {
    
    if (roaring_enabled) {
        return roaring_backend->bitCooccurrence_subset(col_indices, threshold, threads);
    } else {
        throw std::runtime_error("Bit co-occurrence is not supported in BitMatrix backend.");
    }
}


py::array_t<double> HybridMatrix::klDivergence( const std::vector<size_t>& input_guids ) const {
    return bit_backend->klDivergence(input_guids);
}


py::array_t<double> HybridMatrix::jsDivergence( const std::vector<size_t>& input_guids ) const {
    return bit_backend->jsDivergence(input_guids);
}


py::array_t<double> HybridMatrix::informationGain( const std::vector<size_t>& input_guids ) const {
    return bit_backend->informationGain(input_guids);
}


py::array_t<size_t> HybridMatrix::getSetBitIndices( size_t row_idx,
                                                    const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->getSetBitIndices(row_idx);
    } else {
        return bit_backend->getSetBitIndices(row_idx);
    }
}


py::array_t<int> HybridMatrix::getRowMasses( void ) const {
    return bit_backend->getRowMasses();
}


py::array_t<int> HybridMatrix::getColumnMasses( void ) const {
    return bit_backend->getColumnMasses();
}


py::array_t<uint8_t> HybridMatrix::getBitMatrix( void ) const {
    return bit_backend->getBitMatrix();
}


py::list HybridMatrix::getRoaringMatrix( void ) const {
    return roaring_backend->getRoaringMatrix();
}


bool HybridMatrix::roaringEnabled( void ) const {
    return roaring_enabled;
}


py::array_t<uint64_t> HybridMatrix::getPackedMatrix( void ) const {
    return bit_backend->getPackedMatrix();
}


py::array_t<size_t> HybridMatrix::getSubsetPackedMatrix( const std::vector<size_t>& row_indices, 
                                                         const std::vector<size_t>& col_indices,
                                                         const int threads ) const {
    return bit_backend->getSubsetMatrix(row_indices, col_indices, threads);
}

} // namespace ardal