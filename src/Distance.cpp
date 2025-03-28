/*
Copyright 2025 Arthur V. Morris
*/

#include "src/Distance.hpp"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace py = pybind11;
namespace _ardal {



/****************************************************************************************************
 * _ardal::Distance::hamming
 * 
 * Calculate the Hamming distances between all pairs of rows.
 *
 * This function calculates the pairwise Hamming distances between all rows of the allele matrix.
 * The Hamming distance between two rows is the number of positions at which the corresponding
 * elements differ.  The results are returned as a condensed distance matrix.
 *
 * INPUT: None (operates on the private member AlleleMatrix instance _allele_matrix)
 *
 * OUTPUT:
 *  py::array_t<int> : A 1D NumPy array representing the condensed distance matrix containing
 *                      the pairwise Hamming distances.  The length of the array is n*(n-1)/2,
 *                      where 'n' is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<int> Distance::hamming( void ) const {
    // access matrix (read only)
    auto matrix_acc = _allele_matrix.getMatrix().unchecked<2>();

    // get dist matrix size
    size_t n = _allele_matrix.getNumRows();
    size_t m = _allele_matrix.getNumCols();
    size_t dk = (n * (n - 1)) / 2;

    // initialise with mutable access
    py::array_t<int> dist_matrix(dk);
    auto dist_matrix_acc = dist_matrix.mutable_unchecked<1>();

    size_t k = 0;  // index for dist_matrix

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {  // iterate over rows
            int distance;

            // _cache lookup
            auto cached_dist = _cache.get(i, j);

            if (cached_dist != -1) {
                distance = cached_dist;  // dist found and assigned
            } else {
                distance = 0;  // dist not found so initialised as 0

                // hamming distance calculation
                for (size_t l = 0; l < m; ++l) {
                    if (matrix_acc(i, l) != matrix_acc(j, l)) {
                        distance++;
                    }
                }
                // store in cache
                _cache.put(i, j, distance);
            }
            dist_matrix_acc(k++) = distance;
        }
    }
    return dist_matrix;
}



/****************************************************************************************************
 * _ardal::Distance::jaccard
 * 
 * Calculate the Jaccard distances between all pairs of rows.
 *
 * This function calculates the pairwise Jaccard distances between all rows of the allele matrix.
 * The Jaccard distance between two rows is defined as 1 - (Intersection / Union), where
 * Intersection is the number of alleles present in both rows, and Union is the number of alleles
 * present in either row.  The result is returned as a condensed distance matrix.
 *
 * INPUT: None (operates on the private member AlleleMatrix instance _allele_matrix)
 *
 * OUTPUT:
 *  py::array_t<double> : A 1D NumPy array representing the condensed distance matrix containing the 
 *                        pairwise Jaccard distances. The length of the array is n*(n-1)/2, where 'n' 
 *                        is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<double> Distance::jaccard( void ) const {
    // access matrix (read only)
    auto matrix_acc = _allele_matrix.getMatrix().unchecked<2>();

    // get dist matrix size
    size_t n = _allele_matrix.getNumRows();
    size_t m = _allele_matrix.getNumCols();
    size_t dk = (n * (n - 1)) / 2;

    // initialise with mutable access
    py::array_t<double> dist_matrix(dk);  // Use double for Jaccard distance
    auto dist_matrix_acc = dist_matrix.mutable_unchecked<1>();

    size_t k = 0;  // index for dist_matrix

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {  // iterate over rows
            // initialise variables
            int union_n = 0;
            int intersection_n = 0;

            for (size_t l = 0; l < m; ++l) {  // iterate over columns
                // calculate union and intersection
                if (matrix_acc(i, l) == 1 || matrix_acc(j, l) == 1) {
                    union_n++;
                    if (matrix_acc(i, l) == 1 && matrix_acc(j, l) == 1) {
                        intersection_n++;
                    }
                }
            }

            // calculate the Jaccard distance and update the distance matrix
            double dist = (union_n == 0) ? 0.0 : 1 - (static_cast<double>(intersection_n) / union_n);
            dist_matrix_acc(k++) = dist;  // store and increment k
        }
    }

    return dist_matrix;
}

}  // namespace _ardal
