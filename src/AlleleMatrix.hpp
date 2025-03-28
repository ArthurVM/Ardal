#ifndef ALLELE_MATRIX_HPP
#define ALLELE_MATRIX_HPP

#include <iostream>
#include <stdexcept>
#include <limits>
#include <map>
#include <immintrin.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace _ardal {

class AlleleMatrix {
 public:

    // Constructor
    AlleleMatrix(pybind11::array_t<uint8_t> matrix) : _matrix(matrix) {

        // dimensionality Check
        if (_matrix.ndim() != 2) {  // Access the member _matrix now.
            throw std::runtime_error("Matrix must be 2D.");
        }

        // capture dimensions
        _n = _matrix.shape(0);
        _m = _matrix.shape(1);

        // size check
        if (_n > std::numeric_limits<size_t>::max() / _m) {
            throw std::runtime_error("Matrix dimensions are too large, potential overflow.");
        }

        // data type check
        if (_matrix.dtype().kind() != 'u' || _matrix.dtype().itemsize() != 1) {
            throw std::runtime_error("Matrix must have uint8 data type.");
        }

        // get the mass of each row, referring to the number of alleles each guid exhibits
        _rmass = _mass();

        // memory contiguity check
        // if (!_matrix.check()) {
        //     throw std::runtime_error("Matrix must be contiguous in memory. np.ascontiguousarray(matrix) is recommended.");
        // }
    }

    ~AlleleMatrix() { /* AlleleMatrix class destructor */ }

    // access methods
    pybind11::array_t<uint8_t> access( pybind11::array_t<size_t> coords );
    std::set<int> accessGUID( int guid_index ) const;
    pybind11::array_t<uint8_t> getMatrix( void ) const;
    size_t getNumRows( void ) const;
    size_t getNumCols( void ) const;
    int getMass( int ) const;
    std::vector<int> gatherSNPs( const pybind11::array_t<int> guid_indices ) const;


 private:
    pybind11::array_t<uint8_t> _matrix;  // the matrix data
    size_t _n;                           // number of rows
    size_t _m;                           // number of columns
    std::vector<int> _rmass;             // the mass of each row

    // hamming distance cache
    mutable std::map<std::pair<int, int>, int> _hamming_cache;

    // row mass function
    std::vector<int> _mass( void ) const;

};  // class AlleleMatrix

}  // namespace _ardal

#endif  // SRC_ALLELEMATRIX_HPP_

