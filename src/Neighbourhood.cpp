/*
Copyright 2025 Arthur V. Morris
*/

#include "src/Neighbourhood.hpp"

namespace py = pybind11;
namespace _ardal {



/****************************************************************************************************
 * _ardal::Neighbourhood::neighbourhood
 * 
 * Find the epsilon-Neighbourhood of a row using Hamming distance.
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
 *                       epsilon-Neighbourhood of the target row.
 *
 * EXCEPTIONS:
 *   std::runtime_error : If row_coord is out of range.
 ****************************************************************************************************/
py::array_t<int> Neighbourhood::neighbourhood( size_t row_coord, int epsilon ) const {
    // get matrix dimensions
    size_t n = _allele_matrix.getNumRows();
    size_t m = _allele_matrix.getNumCols();

    // do some data cleanliness
    if (epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_coord < 0) {
        throw std::runtime_error("row_coord must be non-negative.");
        }
    if (row_coord >= n) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }

    py::list ep_n;
    
    // access matrix (read only)
    auto matrix_acc = _allele_matrix.getMatrix().unchecked<2>();
    
    // access pre-calculated row mass
    int q_mass = _allele_matrix.getMass(row_coord);

    for (size_t i = 0; i < n; ++i) {
        if (i != row_coord) {
            // row mass filter
            int i_mass = _allele_matrix.getMass(i);
            int mass_d = std::abs(q_mass - i_mass);
            if (mass_d > epsilon) {
                continue;  // skip hamming dist calculation
            }

            // hamming dist calculation
            int distance;                 // initialise distance
            bool ep_exceeded = false;     // flag which stores whether the epsilon neighbourhood was exceeded for early exit


            // check _cache
            auto cached_dist = _cache.get(row_coord, i);  // access cache

            // check if cached_dist exists
            // if it does then dont bother with distance calculation
            if (cached_dist != -1) {
                distance = cached_dist;
            } else {
                distance = 0;
                for (size_t j = 0; j < m; ++j) {
                    if (matrix_acc(row_coord, j) != matrix_acc(i, j)) {
                        distance++;
                        if (distance > epsilon) {
                            ep_exceeded = true;
                            break;  // early break where epsilon is exceeded
                        }
                    }
                }
            }

            if (ep_exceeded) {
                continue;  // don't check against epsilon again
            }

            // DEBUGGING
            // std::cout << "qmass: " << q_mass << " " << i << " (" << i_mass << ") : " << distance << std::endl;

            if (distance <= epsilon) {
                _cache.put(row_coord, i, distance);  // cache results
                ep_n.append(py::make_tuple(row_coord, i, distance));   // Append tuple
            }
        }
    }

    return ep_n;
}



/****************************************************************************************************
 * _ardal::Neighbourhood::neighbourhoodSIMD
 * 
 * Find the epsilon-Neighbourhood of a row using Hamming distance (SIMD optimized).
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
 *                       epsilon-Neighbourhood of the target row.
 *
 * EXCEPTIONS:
 *   std::runtime_error : If row_coord is out of range.
 ****************************************************************************************************/
py::list Neighbourhood::neighbourhoodSIMD( size_t row_coord, int epsilon ) const {
    // get matrix dimensions
    size_t n = _allele_matrix.getNumRows();
    size_t m = _allele_matrix.getNumCols();

    // do some data cleanliness
    if (epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_coord < 0) {
        throw std::runtime_error("row_coord must be non-negative.");
        }
    if (row_coord >= n) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }

    // access matrix (read only)
    auto matrix_acc = _allele_matrix.getMatrix().unchecked<2>();

    py::list ep_n;

    // access pre-calculated row mass
    int q_mass = _allele_matrix.getMass(row_coord);

    for (size_t i = 0; i < n; ++i) {
        if (i != row_coord) {
            // row mass filter
            int i_mass = _allele_matrix.getMass(i);
            int mass_d = std::abs(q_mass - i_mass);
            if (mass_d > epsilon) {
                continue;  // skip hamming dist calculation
            }

            int distance;                 // initialise distance
            bool ep_exceeded = false;     // flag which stores whether the epsilon neighbourhood was exceeded for early exit
            bool cached = false;          // flag which stores whether a cached distance is being used

            // check _cache
            auto cached_dist = _cache.get(row_coord, i);  // access cache

            if (cached_dist != -1) {
                distance = cached_dist;
                cached = true;
            } else {
                distance = 0;
                // hamming distance using SIMD (AVX2)
                size_t j = 0;
                for (; j + 31 < m; j += 32) {  // process in chunks of 32
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
                        ep_exceeded = true;
                        break;  // early break where epsilon is exceeded
                    }
                }

                // calculate hamming distance for any remaining elements:
                if(j < m && !ep_exceeded) {
                    for (; j < m; ++j) {
                        if (matrix_acc(row_coord, j) != matrix_acc(i, j)) {
                            distance++;
                        }
                        // early break if epsilon is exceeded
                        if (distance > epsilon) {
                            ep_exceeded = true;
                            break;
                        }
                    }
                }
            }

            if (ep_exceeded) {
                continue;  // don't check against epsilon again
            }

            // DEBUGGING
            // std::cout << "qmass: " << q_mass << " " << i << " (" << _rmass[i] << ") : " << distance << std::endl;

            if (distance <= epsilon) {
                if (!cached) { _cache.put(row_coord, i, distance); }  // cache result
                ep_n.append(py::make_tuple(i, distance));         // append results tuple
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

} // namespace _ardal