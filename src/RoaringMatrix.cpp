/*
Copyright 2025 Arthur V. Morris
*/
#include "RoaringMatrix.hpp"
#include <stdexcept>
#include <cmath>
#include <cstring>

namespace py = pybind11;
namespace _ardal {


/****************************************************************************************************
 * ardal::RoaringMatrix::RoaringMatrix
 *
 * Constructor for the RoaringMatrix class.
 *
 * This constructor initializes the RoaringMatrix object by taking a 2D NumPy array of uint8_t
 * (representing a binary matrix) and converting each row into a Roaring bitmap. This
 * representation is particularly efficient for sparse binary data.
 *
 * INPUT:
 *  input_matrix (py::array_t<uint8_t>): A 2D NumPy array containing only binary values (0 or 1).
 *
 * OUTPUT:
 *  None (constructor)
 *
 * EXCEPTIONS:
 *  std::runtime_error:
 *   - If the input matrix is not 2-dimensional.
 *   - If the input matrix contains values other than 0 or 1.
 ****************************************************************************************************/
RoaringMatrix::RoaringMatrix( py::array_t<uint8_t> input_matrix,
                              const std::vector<int>& row_masses )
    : _row_masses(row_masses) {
    auto buf = input_matrix.request();

    std::cout << "Building roaring" << std::endl;

    if (buf.ndim != 2) {
        throw std::runtime_error("Input matrix must be 2-dimensional");
    }

    // capture matrix dimensions
    _n_rows = buf.shape[0];
    _n_cols = buf.shape[1];

    // allocate memory for roaring bitmaps
    _roaring_matrix.resize(_n_rows);
    
    // populate roaring bitmaps
    auto ptr = static_cast<uint8_t*>(buf.ptr);
    for (size_t i = 0; i < _n_rows; ++i) {
        roaring::Roaring& bitmap = _roaring_matrix[i];
        for (size_t j = 0; j < _n_cols; ++j) {
            uint8_t val = ptr[i * _n_cols + j];
            if (val != 0 && val != 1) {
                throw std::runtime_error("Input matrix must only contain binary values (0 or 1)");
            }
            if (val) {
                bitmap.add(j);
            }
        }
    }
}


/****************************************************************************************************
 * ardal::RoaringMatrix::hamming
 *
 * Calculates the Hamming distances between all pairs of rows using Roaring bitmaps.
 *
 * This function calculates the pairwise Hamming distances between all rows of the matrix.
 * The Hamming distance between two rows (Roaring bitmaps) is computed as the cardinality of
 * their symmetric difference (A XOR B).
 *
 * INPUT:
 *  None (operates on the private member _roaring_matrix)
 *
 * PARAMETERS:
 *  threads (int) : The number of threads to use for parallel computation.
 *
 * OUTPUT:
 *  py::array_t<int> : A 1D NumPy array representing the condensed distance matrix containing
 *                     the pairwise Hamming distances. The length of the array is n*(n-1)/2,
 *                     where 'n' is the number of rows in the matrix.
 ****************************************************************************************************/
py::array_t<uint32_t> RoaringMatrix::hamming( int threads ) const {

    std::cout << "Calculating hamming" << std::endl;

    const size_t total_pairs = _n_rows * (_n_rows - 1) / 2;
    py::array_t<uint32_t> dist_matrix(total_pairs);     
    {
        py::gil_scoped_release release;
        omp_set_num_threads(threads);

        auto dist_ptr = dist_matrix.mutable_data();

        #pragma omp parallel for schedule(guided)
        for (size_t pair_idx = 0; pair_idx < total_pairs; ++pair_idx) {
            size_t i = floor((2.0 * _n_rows - 1 - sqrt(pow(2.0 * _n_rows - 1, 2) - 8.0 * pair_idx)) / 2.0);
            size_t j = pair_idx + i + 1 - _n_rows * i + i * i / 2;

            dist_ptr[pair_idx] = hammingDistance(i, j);
        }
    }     
    return dist_matrix;
}


uint32_t RoaringMatrix::hammingDistance( size_t i, size_t j ) const {
    return (_roaring_matrix[i] ^ _roaring_matrix[j]).cardinality();
}


py::array_t<int64_t> RoaringMatrix::neighbourhood( size_t row_idx, int epsilon, int threads ) const {
    // do some data cleanliness
    if (epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_idx < 0) {
        throw std::runtime_error("Row index must be non-negative.");
        }
    if (row_idx >= _n_rows) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }
    if (threads <= 0)
        throw std::runtime_error("Number of threads must be positive.");
    
    int q_mass = _row_masses[row_idx];    // access pre-calculated row mass
    std::vector<std::vector<std::pair<size_t, int>>> thread_results;

    {
        py::gil_scoped_release release;   // release GIL for full parallel region
        omp_set_num_threads(threads);     // Explicitly control number of threads

        #pragma omp parallel
        {
            const int thread_id = omp_get_thread_num();

            #pragma omp single
            thread_results.resize(omp_get_num_threads());

            auto& local = thread_results[thread_id];

            #pragma omp for schedule(static)
            for (size_t i = 0; i < _n_rows; ++i) {
                if (i == row_idx) continue;

                int mass_d = std::abs(q_mass - _row_masses[i]);
                if (mass_d > epsilon) continue;
                
                // due to the restrictions on using a roaring bitmap, early exit is not trivial :(
                int distance = hammingDistance(row_idx, i);

                if (distance <= epsilon)
                    local.emplace_back(i, distance);
            }
        }
    }   // GIL reestablished here

    // count total neighbours to preallocate np array
    size_t total_neighbours = 0;
    for (const auto& vec : thread_results) {
        total_neighbours += vec.size();
    }

    // create result np array of shape (total_neighbours, 2)
    py::array_t<int64_t> ep_n({total_neighbours, (size_t)2});
    auto result_ptr = ep_n.mutable_data();

    // populate array with (idex, dist) pairs
    size_t curr_idx = 0;
    for (const auto& vec : thread_results) {
        for (const auto& [idx, dist] : vec) {
            result_ptr[curr_idx * 2 + 0] = idx;
            result_ptr[curr_idx * 2 + 1] = dist;
            curr_idx++;
        }
    }       
    return ep_n;
}


py::array_t<int> RoaringMatrix::innerProduct( int threads ) const {
    py::array_t<int> inner_product_matrix;
    return inner_product_matrix;
}

py::list RoaringMatrix::innerProductNeighbourhood( size_t row_idx, int ip_epsilon, int threads ) const {
    // do some data cleanliness
    if (ip_epsilon < 0) {
        throw std::runtime_error("epsilon must be non-negative.");
        }
    if (row_idx < 0) {
        throw std::runtime_error("Row index must be non-negative.");
        }
    if (row_idx >= _n_rows) {
        throw std::runtime_error("Coordinate dimensions exceed the number of rows.");
        }
    if (threads <= 0)
        throw std::runtime_error("Number of threads must be positive.");
    
    std::vector<std::vector<std::pair<size_t, int>>> thread_results;

    {
        py::gil_scoped_release release;   // release GIL for full parallel region
        omp_set_num_threads(threads);     // Explicitly control number of threads

        #pragma omp parallel
        {
            const int thread_id = omp_get_thread_num();

            #pragma omp single
            thread_results.resize(omp_get_num_threads());

            auto& local = thread_results[thread_id];

            #pragma omp for schedule(static)
            for (size_t i = 0; i < _n_rows; ++i) {
                if (i == row_idx) continue;
                
                int distance = innerProductRowwise(row_idx, i);

                if (distance >= ip_epsilon)
                    local.emplace_back(i, distance);
            }
        }
    }   // GIL reestablished here

    // count total neighbours to preallocate np array
    size_t total_neighbours = 0;
    for (const auto& vec : thread_results) {
        total_neighbours += vec.size();
    }

    // create result np array of shape (total_neighbours, 2)
    py::array_t<int64_t> ipe_n({total_neighbours, (size_t)2});
    auto result_ptr = ipe_n.mutable_data();

    // populate array with (idex, dist) pairs
    size_t curr_idx = 0;
    for (const auto& vec : thread_results) {
        for (const auto& [idx, dist] : vec) {
            result_ptr[curr_idx * 2 + 0] = idx;
            result_ptr[curr_idx * 2 + 1] = dist;
            curr_idx++;
        }
    }       
    return ipe_n;
}


uint32_t RoaringMatrix::innerProductRowwise( size_t i, size_t j ) const {
    return (_roaring_matrix[i] & _roaring_matrix[j]).cardinality();
}


std::vector<size_t> RoaringMatrix::uniqueSharedBits( const std::vector<size_t>& row_indices ) const {
    std::vector<size_t> shared_bits;
    return shared_bits;
}


py::array_t<size_t> RoaringMatrix::getSetBitIndices( size_t row_idx ) const {
    const auto& bitmap = _roaring_matrix.at(row_idx);

    // Create a numpy array of uint32_t to hold the values.
    py::array_t<uint32_t> values_u32(bitmap.cardinality());

    // Fill the numpy array with values from the roaring bitmap.
    bitmap.toUint32Array(values_u32.mutable_data());

    // Cast the numpy array to size_t, which is what the python side expects.
    return values_u32.cast<py::array_t<size_t>>();
}


py::list RoaringMatrix::getRoaringMatrix( void ) const {
    py::list roaring_matrix;
    for (size_t i = 0; i < _n_rows; ++i) {
        const auto& bitmap = _roaring_matrix[i];
        roaring_matrix.append(getSetBitIndices(i));
    }
    return roaring_matrix;
}

} // namespace _ardal