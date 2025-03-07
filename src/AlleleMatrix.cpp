/*
Copyright 2025 Arthur V. Morris
*/

#include "AlleleMatrix.hpp"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace py = pybind11;
namespace _ardal {



/****************************************************************************************************
 * ardal::AlleleMatrix::_mass
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
 * ardal::AlleleMatrix::access
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
 * ardal::AlleleMatrix::getMatrix
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
 * ardal::AlleleMatrix::hamming
 * 
 * Calculate the Hamming distances between all pairs of rows.
 *
 * This function calculates the pairwise Hamming distances between all rows of the allele matrix.
 * The Hamming distance between two rows is the number of positions at which the corresponding
 * elements differ.  The results are returned as a condensed distance matrix.
 *
 * INPUT: None (operates on the private member _matrix)
 *
 * OUTPUT:
 *  py::array_t<int> : A 1D NumPy array representing the condensed distance matrix containing
 *                      the pairwise Hamming distances.  The length of the array is n*(n-1)/2,
 *                      where 'n' is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<int> AlleleMatrix::hamming( void ) const {
    // access matrix (read only)
    auto matrix_acc = _matrix.unchecked<2>();

    // calculate dist matrix size
    size_t n = _matrix.shape(0);
    size_t dk = (n * (n - 1)) / 2;

    // initialise with mutable access
    py::array_t<int> dist_matrix(dk);
    auto dist_matrix_acc = dist_matrix.mutable_unchecked<1>();

    size_t k = 0;  // index for dist_matrix

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {  // iterate over rows
            int dist;

            // _hamming_cache lookup
            auto key = std::make_pair(std::min(i, j), std::max(i, j));
            auto cached_dist = _hamming_cache.find(key);

            if (cached_dist != _hamming_cache.end()) {
                dist = cached_dist->second;  // dist found and assigned
            } else {
                dist = 0;  // dist not found so initialised as 0

                // hamming distance calculation
                for (size_t l = 0; l < _m; ++l) {
                    if (matrix_acc(i, l) != matrix_acc(j, l)) {
                        dist++;
                    }
                }
                // store in cache
                _hamming_cache[key] = dist;
            }
            dist_matrix_acc(k++) = dist;
        }
    }
    return dist_matrix;
}



/****************************************************************************************************
 * ardal::AlleleMatrix::jaccard
 * 
 * Calculate the Jaccard distances between all pairs of rows.
 *
 * This function calculates the pairwise Jaccard distances between all rows of the allele matrix.
 * The Jaccard distance between two rows is defined as 1 - (Intersection / Union), where
 * Intersection is the number of alleles present in both rows, and Union is the number of alleles
 * present in either row.  The result is returned as a condensed distance matrix.
 *
 * INPUT: None (operates on the private member _matrix)
 *
 * OUTPUT:
 *  py::array_t<double> : A 1D NumPy array representing the condensed distance matrix containing the 
 *                        pairwise Jaccard distances. The length of the array is n*(n-1)/2, where 'n' 
 *                        is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<double> AlleleMatrix::jaccard( void ) const {
    // access matrix (read only)
    auto matrix_acc = _matrix.unchecked<2>();

    // calculate dist matrix size
    size_t n = _matrix.shape(0);
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

            for (size_t l = 0; l < _m; ++l) {  // iterate over columns
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



/****************************************************************************************************
 * ardal::AlleleMatrix::neighbourhood
 * 
 * Find the epsilon-neighborhood of a row using Hamming distance.
 *
 * This function identifies the rows in the matrix that are within a specified Hamming distance
 * (epsilon) of a given row.
 *
 * INPUT:
 *   row_coord (size_t) : The index of the target row.
 *   epsilon (int)   : The maximum Hamming distance threshold.
 *
 * OUTPUT:
 *   py::array_t<int> : A 1D NumPy array containing the indices of the rows that are within the 
 *                       epsilon-neighborhood of the target row.
 *
 * EXCEPTIONS:
 *   std::runtime_error : If row_coord is out of range.
 ****************************************************************************************************/
py::array_t<int> AlleleMatrix::neighbourhood( size_t row_coord, int epsilon ) const {
    // do some data cleanliness
    if (epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_coord < 0) {
        throw std::runtime_error("row_coord must be non-negative.");
        }
    if (row_coord >= _n) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }

    std::vector<int> ep_n;
    auto matrix_acc = _matrix.unchecked<2>();

    int q_mass = _rmass[row_coord];  // access pre-calculated row mass

    for (size_t i = 0; i < _n; ++i) {
        if (i != row_coord) { 
            // row mass filter
            int mass_d = std::abs(q_mass - _rmass[i]);
            if (mass_d > epsilon) {
                continue;  // skip hamming dist calculation
            }

            // hamming dist calculation
            int distance = 0;
            for (size_t j = 0; j < _m; ++j) {
                if (matrix_acc(row_coord, j) != matrix_acc(i, j)) {
                    distance++;
                    if (distance > epsilon) {
                        break;  // early break where epsilon is exceeded
                    }
                }
            }
            
            // DEBUGGING
            // std::cout << "qmass: " << q_mass << " " << i << " (" << _rmass[i] << ") : " << distance << std::endl;

            if (distance <= epsilon) {
                ep_n.push_back(i);
            }
        }
    }

    // check the neighbourhood isnt empty
    if (ep_n.empty()) {
        // return an np array of size 0 if no neighbors are found
        return py::array_t<int>(0);
    } else {
        return py::array(py::cast(ep_n));
    }
}



/****************************************************************************************************
 * ardal::AlleleMatrix::neighbourhoodSIMD
 * 
 * Find the epsilon-neighborhood of a row using Hamming distance (SIMD optimized).
 *
 * This function identifies the rows in the matrix that are within a specified Hamming distance
 * (epsilon) of a given row, using SIMD (AVX2) intrinsics for optimized performance.
 *
 * INPUT:
 *   row_coord (size_t) : The index of the target row.
 *   epsilon (int)   : The maximum Hamming distance threshold.
 *
 * OUTPUT:
 *   py::array_t<int> : A 1D NumPy array containing the indices of the rows that are within the
 *                       epsilon-neighborhood of the target row.
 *
 * EXCEPTIONS:
 *   std::runtime_error : If row_coord is out of range.
 ****************************************************************************************************/
py::list AlleleMatrix::neighbourhoodSIMD( size_t row_coord, int epsilon ) const {
    // do some data cleanliness
    if (epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_coord < 0) {
        throw std::runtime_error("row_coord must be non-negative.");
        }
    if (row_coord >= _n) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }

    py::list ep_n;
    auto matrix_acc = _matrix.unchecked<2>();

    int q_mass = _rmass[row_coord];  // access pre-calculated row mass

    for (size_t i = 0; i < _n; ++i) {
        if (i != row_coord) {
            // row mass filter
            int mass_d = std::abs(q_mass - _rmass[i]);
            if (mass_d > epsilon) {
                continue;  // skip hamming dist calculation
            }

            int distance;
            bool cached = false;

            // check _hamming_cache
            // since (row, col) will be identical to (col, row), use min and max to construct a single key for cache recovery
            auto key = std::make_pair(std::min(row_coord, i), std::max(row_coord, i));
            auto cached_dist = _hamming_cache.find(key);  // access cache

            if (cached_dist != _hamming_cache.end()) {
                distance = cached_dist->second;
                cached = true;
            }

            bool simd_complete = false;
            if (!cached) {
                // hamming distance using SIMD (AVX2)
                distance = 0;
                size_t j = 0;
                for (; j + 31 < _m; j += 32) {  // process in chunks of 32
                    // NOTE: assumes uint8_t. This assumption *should* be correct, but it will break spectacularly if it isnt
                    __m256i a = _mm256_loadu_si256((__m256i*)&matrix_acc(row_coord, j));
                    __m256i b = _mm256_loadu_si256((__m256i*)&matrix_acc(i, j));
                    __m256i xor_result = _mm256_xor_si256(a, b);

                    // set bit count
                    distance += _mm_popcnt_u32(_mm256_extract_epi32(xor_result, 0));
                    distance += _mm_popcnt_u32(_mm256_extract_epi32(xor_result, 1));
                    distance += _mm_popcnt_u32(_mm256_extract_epi32(xor_result, 2));
                    distance += _mm_popcnt_u32(_mm256_extract_epi32(xor_result, 3));
                    distance += _mm_popcnt_u32(_mm256_extract_epi32(xor_result, 4));
                    distance += _mm_popcnt_u32(_mm256_extract_epi32(xor_result, 5));
                    distance += _mm_popcnt_u32(_mm256_extract_epi32(xor_result, 6));
                    distance += _mm_popcnt_u32(_mm256_extract_epi32(xor_result, 7));

                    // early break if epsilon is exceeded
                    if (distance > epsilon) {
                        _hamming_cache[key] = distance;  // <--- CACHE PARTIAL RESULTS
                        break;
                    }
                }

                // ensure the final distance is cached
                if (j >= _m - (_m % 32)) { simd_complete = true; }
                if (!simd_complete) { _hamming_cache[key] = distance; }
            }

            // calculate hamming distance for any remaining elements:
            if (!simd_complete) {
                for (size_t j = j; j < _m; ++j) {
                    if (matrix_acc(row_coord, j) != matrix_acc(i, j)) {
                        distance++;
                    }
                    // early break if epsilon is exceeded
                    if (distance > epsilon) {
                        break;
                    }
                }
            }
            if (!cached) { _hamming_cache[key] = distance; }

            // DEBUGGING
            // std::cout << "qmass: " << q_mass << " " << i << " (" << _rmass[i] << ") : " << distance << std::endl;

            if (distance <= epsilon) {
                ep_n.append(py::make_tuple(i, distance)); // Append tuple
            }
        }
    }

    // check the neighbourhood isnt empty
    // if (ep_n.empty()) {
    //     // return an np array of size 0 if no neighbors are found
    //     return py::array_t<int>(0);
    // } else {
    //     return py::array(py::cast(ep_n));
    // }

    return ep_n;
}

/****************************************************************************************************
 * ardal::AlleleMatrix::accessGUID
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
 * ardal::AlleleMatrix::gatherSNPs
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
 * ardal::AlleleMatrix::getMass
 *
 * Get the mass of each row in the matrix.
 *
 * This function returns the pre-calculated mass of each row in the matrix. The mass of a row
 * represents the number of alleles (1s) present in that row.
 *
 * INPUT: None (operates on the private member _rmass)
 *
 * OUTPUT:
 *   std::vector<int> : A vector containing the mass of each row in the matrix.
 ****************************************************************************************************/
std::vector<int> AlleleMatrix::getMass( void ) {
    return _rmass;
}

} // namespace _ardal


// Pybind methods
PYBIND11_MODULE(_ardal, m) {  // _ardal module and method bindings
    py::class_<_ardal::AlleleMatrix>(m, "AlleleMatrix")
        .def(py::init<py::array_t<uint8_t>&>())
        .def("access", &_ardal::AlleleMatrix::access, py::arg("coords"))
        .def("accessGUID", &_ardal::AlleleMatrix::accessGUID, py::arg("guid_index"))
        .def("getMatrix", &_ardal::AlleleMatrix::getMatrix)
        .def("getMass", &_ardal::AlleleMatrix::getMass)
        .def("hamming", &_ardal::AlleleMatrix::hamming)
        .def("jaccard", &_ardal::AlleleMatrix::jaccard)
        .def("neighbourhood", &_ardal::AlleleMatrix::neighbourhood, py::arg("row_coord"), py::arg("epsilon"))
        .def("neighbourhoodSIMD", &_ardal::AlleleMatrix::neighbourhoodSIMD, py::arg("row_coord"), py::arg("epsilon"))
        .def("gatherSNPs", &_ardal::AlleleMatrix::gatherSNPs, py::arg("guid_indices"));
}