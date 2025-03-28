/*
Copyright 2025 Arthur V. Morris
*/

#include "src/AlleleMatrix.hpp"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace py = pybind11;
namespace _ardal {



/****************************************************************************************************
 * _ardal::AlleleMatrix::_mass
 *
 * Calculate the row masses of the matrix.
 *
 * This private helper function calculates the "mass" of each row in the matrix. In the context
 * of a binary allele matrix, the mass of a row represents the number of alleles (1s) present
 * in that row.
 *
 * INPUT: None (operates on the private member _matrix)
 *
 * OUTPUT:
 *  std::vector<int> : A vector containing the mass of each row in the matrix.
 ****************************************************************************************************/
std::vector<int> AlleleMatrix::_mass( void ) const {
    std::vector<int> row_mass;
    row_mass.reserve(_n);

    auto matrix_acc = _matrix.unchecked<2>();

    for (size_t i = 0; i < _n; ++i) {
        int mass = 0;
        for (size_t j = 0; j < _m; ++j) {
            mass += matrix_acc(i, j);
        }
        row_mass.push_back(mass);
    }
    return row_mass;
}



/****************************************************************************************************
 * _ardal::AlleleMatrix::access
 *
 * Access elements in a 2D matrix at specified coordinates.
 *
 * This function efficiently retrieves elements from a 2D matrix (represented as a NumPy array)
 * at given coordinates. It performs bounds checking to ensure all coordinates are within
 * the matrix dimensions.
 *
 * INPUT:
 * coords (py::array_t<size_t>) : A 2D array of coordinates (shape (k, 2)), where each row
 *                             represents a (row, col) coordinate pair.
 *
 * OUTPUT:
 * result (py::array_t<uint8_t>) : A 1D NumPy array containing the elements at the specified
 *                                coordinates. The length of the array is equal to the number
 *                                of coordinate pairs provided.
 *
 * EXCEPTIONS:
 * std::runtime_error : If the input matrix is not 2D, if the coordinates array is not 2D
 *                       or does not have a shape of (k, 2), or if any coordinate is out of bounds.
 ****************************************************************************************************/
py::array_t<uint8_t> AlleleMatrix::access( py::array_t<size_t> coords ) {
    // check coordinate dimensions
    if (coords.ndim() != 2 || coords.shape(1) != 2) {
        throw std::runtime_error("Coordinates must be 2D with shape (k, 2).");
    }

    size_t k = coords.shape(0);

    auto matrix_acc = _matrix.unchecked<2>();
    auto coords_acc = coords.unchecked<2>();

    // calculate min and max row/col
    size_t min_row = coords_acc(0, 0);
    size_t max_row = min_row;
    size_t min_col = coords_acc(0, 1);
    size_t max_col = min_col;

    for (size_t i = 1; i < k; ++i) {
        min_row = std::min(min_row, static_cast<size_t>(coords_acc(i, 0)));
        max_row = std::max(max_row, static_cast<size_t>(coords_acc(i, 0)));
        min_col = std::min(min_col, static_cast<size_t>(coords_acc(i, 1)));
        max_col = std::max(max_col, static_cast<size_t>(coords_acc(i, 1)));
    }

    // bounds check
    if (min_row >= _n || max_row >= _n || min_col >= _m || max_col >= _m) {
        throw std::runtime_error("Coordinates out of range.");
    }

    auto result = py::array_t<uint8_t>(k);
    auto result_acc = result.mutable_unchecked<1>();

    for (size_t i = 0; i < k; ++i) {
        result_acc(i) = matrix_acc(coords_acc(i, 0), coords_acc(i, 1));
    }

    return result;
}



/****************************************************************************************************
 * _ardal::AlleleMatrix::getMatrix
 * 
 * Return the allele matrix.
 *
 * OUTPUT:
 *  py::array_t<uint8_t> : A 2D numpy array representing a binary allele matrix.
 ****************************************************************************************************/
py::array_t<uint8_t> AlleleMatrix::getMatrix( void ) const {
    return _matrix;
}



/****************************************************************************************************
 * _ardal::AlleleMatrix::getNumRows
 * 
 * Return the number of rows in the allele matrix.
 *
 * OUTPUT:
 *  size_t : The number of rows in the allele matrix
 ****************************************************************************************************/
size_t AlleleMatrix::getNumRows( void ) const {
    return _n;
}



/****************************************************************************************************
 * _ardal::AlleleMatrix::getNumCols
 * 
 * Return the number of columns in the allele matrix.
 *
 * OUTPUT:
 *  size_t : The number of columns in the allele matrix
 ****************************************************************************************************/
size_t AlleleMatrix::getNumCols( void ) const {
    return _m;
}



/****************************************************************************************************
 * _ardal::AlleleMatrix::accessGUID
 *
 * Access the set of SNPs for a given GUID.
 *
 * This function retrieves the set of SNPs (Single Nucleotide Polymorphisms) that are present
 * in a specified GUID (row) of the allele matrix.
 *
 * INPUT:
 *   guid_index (int) : The index of the GUID (row) in the matrix.
 *
 * OUTPUT:
 *   std::set<int> : A set containing the indices of the SNPs present in the specified GUID.
 *
 * EXCEPTIONS:
 *   None (returns an empty set if the guid_index is out of bounds)
 ****************************************************************************************************/
std::set<int> AlleleMatrix::accessGUID( int guid_index ) const {
    std::set<int> guid_snps;
    if (guid_index >= 0 && guid_index < _n) {
        for (size_t j = 0; j < _m; ++j) {
            if (_matrix.at(guid_index, j)) {
                guid_snps.insert(j);
            }
        }
    }
    return guid_snps;
}



/****************************************************************************************************
 * _ardal::AlleleMatrix::gatherSNPs
 *
 * Gather SNPs from multiple GUIDs.
 *
 * This function collects the SNPs from a set of specified GUIDs (rows) in the allele matrix.
 * It iterates through the provided GUID indices, retrieves the SNPs for each GUID, and
 * aggregates them into a single vector.
 *
 * INPUT:
 *   guid_indices (const py::array_t<int>) : A NumPy array containing the indices of the GUIDs.
 *
 * OUTPUT:
 *   std::vector<int> : A vector containing the indices of all SNPs present in the specified GUIDs.
 *
 * EXCEPTIONS:
 *   None (returns an empty vector if no GUID indices are provided)
 ****************************************************************************************************/
std::vector<int> AlleleMatrix::gatherSNPs( const py::array_t<int> guid_indices ) const {
    std::vector<int> snp_vector;
    if (guid_indices.size() > 0) {
        py::buffer_info guid_buf = guid_indices.request();
        int* guids = (int *) guid_buf.ptr;

        // iterate through the rest of the guids
        for (size_t i = 0; i < guid_indices.size(); ++i) {
            // get the SNPs for the current guid
            std::set<int> guid_snps = accessGUID(guids[i]);

            // append them to the vector
            for (int snp : guid_snps) {
                snp_vector.push_back(snp);
            }
        }
    }
    return snp_vector;
}



/****************************************************************************************************
 * _ardal::AlleleMatrix::getMass
 *
 * Get the mass of a row in the matrix.
 *
 * This function returns the pre-calculated mass of a row in the matrix. The mass of a row
 * represents the number of alleles (1s) present in that row.
 *
 * INPUT: 
 *   guid_index (int) : The index for the query guid.
 *
 * OUTPUT:
 *   int : The mass of the row defined by the guid index in the matrix.
 * 
 * EXCEPTIONS:
 *   std::runtime_error : If the guid index is out of range.
 ****************************************************************************************************/
int AlleleMatrix::getMass( int guid_index ) const {
    if (guid_index >= 0 && guid_index < _n) {
        return _rmass[guid_index];
    }
    else {
        throw std::runtime_error("guid_index out of range.");
    }
}

} // namespace _ardal