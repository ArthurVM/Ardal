#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>
#include <iostream>

namespace py = pybind11;

/****************************************************************************************************
 * ACCESSALLELEMATRIX
 * 
 * Access elements in a 2D matrix at specified coordinates.
 *
 * This function efficiently retrieves elements from a 2D matrix (represented as a NumPy array)
 * at given coordinates. It performs bounds checking to ensure all coordinates are within
 * the matrix dimensions.
 *
 * INPUT:
 *  matrix (py::array_t<uint8_t>) : The 2D matrix (NumPy array) to access.
 *  coords (py::array_t<int>)     : A 2D array of coordinates (shape (k, 2)), where each row
 *                                  represents a (row, col) coordinate pair.
 *
 * OUTPUT:
 *  result (py::array_t<uint8_t>) : A 1D NumPy array containing the elements at the specified
 *                                 coordinates. The length of the array is equal to the number
 *                                 of coordinate pairs provided.
 *
 * EXCEPTIONS:
 *  std::runtime_error : If the input matrix is not 2D, if the coordinates array is not 2D
 *                        or does not have a shape of (k, 2), or if any coordinate is out of bounds.
 ****************************************************************************************************/
py::array_t<uint8_t> accessAlleleMatrix( py::array_t<uint8_t> matrix, py::array_t<int> coords ) {

    auto func_start = std::chrono::high_resolution_clock::now();

    // Check input dimensions
    if (matrix.ndim() != 2) {
        throw std::runtime_error("Matrix must be 2D.");
    }
    if (coords.ndim() != 2 || coords.shape(1) != 2) {
        throw std::runtime_error("Coordinates must be 2D with shape (k, 2).");
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    // Get matrix dimensions
    size_t n = matrix.shape(0);
    size_t m = matrix.shape(1);
    size_t k = coords.shape(0); // Number of queries

    auto t2 = std::chrono::high_resolution_clock::now();

    // Access data (read-only)
    auto matrix_acc = matrix.unchecked<2>();
    auto coords_acc = coords.unchecked<2>();

    auto t3 = std::chrono::high_resolution_clock::now();

    // Pre-calculate min and max row/col
    int min_row = coords_acc(0, 0);
    int max_row = min_row;
    int min_col = coords_acc(0, 1);
    int max_col = min_col;

    auto t4 = std::chrono::high_resolution_clock::now();

    for (size_t i = 1; i < k; ++i) {
        min_row = std::min(min_row, coords_acc(i, 0));
        max_row = std::max(max_row, coords_acc(i, 0));
        min_col = std::min(min_col, coords_acc(i, 1));
        max_col = std::max(max_col, coords_acc(i, 1));
    }

    auto t5 = std::chrono::high_resolution_clock::now();

    // Perform bounds check *once*
    if (min_row < 0 || max_row >= n || min_col < 0 || max_col >= m) {
        throw std::runtime_error("Coordinates out of range."); // More concise error message
    }

    auto t6 = std::chrono::high_resolution_clock::now();

    // Create output array
    auto result = py::array_t<uint8_t>(k);
    auto result_acc = result.mutable_unchecked<1>();

    auto t7 = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < k; ++i) {
        int row = coords_acc(i, 0);
        int col = coords_acc(i, 1);
        result_acc(i) = matrix_acc(row, col);   

    }

    auto func_end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(func_end - func_start);
    
    std::cout << "TIMINGS:" << std::endl;
    std::cout << "T1-T2: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << " nanoseconds" << std::endl
              << "T2-T3: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() << " nanoseconds" << std::endl
              << "T3-T4: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count() << " nanoseconds" << std::endl
              << "T4-T5: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t4).count() << " nanoseconds" << std::endl
              << "T5-T6: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t6 - t5).count() << " nanoseconds" << std::endl
              << "T6-T7: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t7 - t6).count() << " nanoseconds" << std::endl
              << "T7-END: " << std::chrono::duration_cast<std::chrono::nanoseconds>(func_end - t7).count() << " nanoseconds" << std::endl
              << "Total: " << elapsed.count() << " nanoseconds" << std::endl;

    return result;
}


PYBIND11_MODULE(ardal, m) {
    m.def("accessAlleleMatrix", &accessAlleleMatrix, "Efficient allele lookups");
}