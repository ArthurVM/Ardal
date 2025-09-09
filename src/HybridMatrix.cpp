/*
Copyright 2025 Arthur V. Morris
*/
#include "HybridMatrix.hpp"
#include "BitMatrix.hpp"
#include "RoaringMatrix.hpp"
#include <iostream>

namespace _ardal {


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
HybridMatrix::HybridMatrix( py::array_t<uint8_t> packed_npy_matrix,
                            bool use_roaring_if_sparse,
                            double density_threshold ) {
    auto buf = packed_npy_matrix.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Input matrix must be 2-dimensional");
    }

    // repack into expected structure for C++
    _packed_matrix.resize(n_rows, std::vector<uint64_t>(_n_packed_cols, 0));
    for (size_t i = 0; i < n_rows; ++i) {
        for (size_t j = 0; j < _n_packed_cols; ++j) {
            _packed_matrix[i][j] = data[i * _n_packed_cols + j];
        }
    }

    // get the matrix metadata once to pass to each backend
    // shape
    size_t n_rows = buf.shape[0];
    size_t n_packed_cols = buf.shape[1];
    size_t n_cols = n_packed_cols * 64;

    // col/row masses
    std::vector<int> row_masses(n_rows, 0);
    std::vector<int> col_masses(n_cols, 0);
    int64_t total_mass = 0;
    
    // matrix density
    const uint64_t* data = static_cast<const uint64_t*>(buf.ptr);
    for (size_t i = 0; i < n_rows; ++i) {
        for (size_t j = 0; j < n_packed_cols; ++j) {
            uint64_t word = data[i * n_packed_cols + j];
            for (size_t b = 0; b < 64; ++b) {
                size_t col = j * 64 + b;
                if (col >= n_cols) break;
                if ((word >> b) & 1) {
                    row_masses[i]++;
                    col_masses[col]++;
                    total_mass++;
                }
            }
        }
    }
    double density = static_cast<double>(total_mass) / (n_rows * n_cols);

    // initialise bit backend
    bit_backend = std::make_unique<_ardal::BitMatrix>(packed_matrix, row_masses, col_masses, n_rows, n_cols);

    // initialise roaring backend if required and density thresholds are met
    if (use_roaring_if_sparse && density < density_threshold) {
        roaring_enabled = true;
        roaring_backend = std::make_unique<_ardal::RoaringMatrix>(packed_matrix, row_masses, col_masses, n_rows, n_cols);
    } else {
        roaring_enabled = false;
    }
}


HybridMatrix::BackendType HybridMatrix::selectBackend(size_t n_cols, double density) const {
    // heuristic backend selection
    if (n_cols > 20000 && density < 0.02) {
        return BackendType::ROARING;
    } else {
        return BackendType::BIT;
    }
}


double HybridMatrix::getDensity( void ) const {
    return density;
}


py::array_t<uint32_t> HybridMatrix::hamming( bool use_simd,
                                             int threads,
                                             const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(_n_cols, density);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->hamming(threads);
    } else {
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
        chosen_backend = selectBackend(row_indices.size(), density);
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
        chosen_backend = selectBackend(_n_cols, density);
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
        chosen_backend = selectBackend(_n_cols, density);
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
        chosen_backend = selectBackend(_n_cols, density);
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
        chosen_backend = selectBackend(_n_cols, density);
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
        chosen_backend = selectBackend(_n_cols, density);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->getSetBitIndices(row_idx);
    } else {
        return bit_backend->getSetBitIndices(row_idx);
    }
}


py::array_t<int> HybridMatrix::getRowMasses( void ) {
    return bit_backend->getRowMasses();
}


py::array_t<int> HybridMatrix::getColumnMasses( void ) {
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


py::array_t<size_t> HybridMatrix::getSubsetPackedMatrix( const std::vector<size_t>& row_indices, 
                                                         const std::vector<size_t>& col_indices,
                                                         const int threads ) const {
    return bit_backend->getSubsetMatrix(row_indices, col_indices, threads);
}

} // namespace _ardal