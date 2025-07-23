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
 *  matrix (py::array_t<uint8_t>): A 2D NumPy array containing only binary values (0 or 1).
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
HybridMatrix::HybridMatrix( py::array_t<uint8_t> matrix,
                            bool use_roaring_if_sparse,
                            double density_threshold ) {
    bit_backend = std::make_unique<_ardal::BitMatrix>(matrix);
    density = bit_backend->getDensity();

    if (use_roaring_if_sparse && density < density_threshold) {
        roaring_enabled = true;
        const auto& row_masses = bit_backend->getRowMassesVector();
        roaring_backend = std::make_unique<_ardal::RoaringMatrix>(matrix, row_masses);
    } else {
        roaring_enabled = false;
    }
}


double HybridMatrix::getDensity() const {
    return density;
}


py::array_t<int> HybridMatrix::hamming( bool use_simd, int threads, bool force_bit_backend ) const {
    if (roaring_enabled && !force_bit_backend) {
        return roaring_backend->hamming(threads);
    } else {
        return bit_backend->hamming(false, use_simd, threads);
    }
}


py::array_t<int> HybridMatrix::innerProduct( bool use_simd, int threads, bool force_bit_backend ) const {
    if (roaring_enabled && !force_bit_backend) {
        return roaring_backend->innerProduct(threads);
    } else {
        return bit_backend->innerProduct(false, use_simd, threads);
    }
}


py::array_t<int64_t> HybridMatrix::neighbourhood( size_t row, int epsilon, bool use_simd, int threads, bool force_bit_backend )  const {
    if (roaring_enabled && !force_bit_backend) {
        return roaring_backend->neighbourhood(row, epsilon, threads);
    } else {
        return bit_backend->neighbourhood(row, epsilon, use_simd, threads);
    }
}


py::list HybridMatrix::innerProductNeighbourhood( size_t row, int ip_epsilon, bool use_simd, bool force_bit_backend ) const {
    if (roaring_enabled && !force_bit_backend) {
        return roaring_backend->innerProductNeighbourhood(row, ip_epsilon);
    } else {
        return bit_backend->innerProductNeighbourhood(row, ip_epsilon, use_simd);
    }
}


std::vector<size_t> HybridMatrix::uniqueSharedBits( const std::vector<size_t>& row_indices, bool use_simd, bool force_bit_backend ) const {
    if (!force_bit_backend) {
        return roaring_backend->uniqueSharedBits(row_indices);
    } else {
        return bit_backend->uniqueSharedBits(row_indices, use_simd);
    }
}


py::array_t<double> HybridMatrix::colFrequency( std::vector<size_t>& row_indices ) const {
    return bit_backend->colFrequency(row_indices);
}


py::array_t<double> HybridMatrix::columnEntropy( void ) const {
    return bit_backend->columnEntropy();
}


py::dict HybridMatrix::bitCooccurrence( double threshold, int threads ) const {
    std::cout << "HM: GOT HERE" << std::endl;
    return roaring_backend->bitCooccurrence(threshold, threads);
}


py::array_t<double> HybridMatrix::klDivergence( const std::vector<size_t>& input_guids ) const {
    return bit_backend->klDivergence(input_guids);
}


py::array_t<size_t> HybridMatrix::getSetBitIndices( size_t row_idx, bool force_bit_backend ) const {
    if (roaring_enabled && !force_bit_backend) {
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


} // namespace _ardal