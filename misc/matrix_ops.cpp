#include <vector>
#include <cstddef>

extern "C" void access_matrix(
    const double* matrix,          // Matrix data pointer
    size_t rows, size_t cols,      // Dimensions of the matrix
    const size_t* coordinates,     // Coordinate vector
    size_t num_coords,             // Number of coordinate pairs
    double* output                 // Output vector
) {
    for (size_t i = 0; i < num_coords; ++i) {
        size_t x = coordinates[2 * i];       // Row index
        size_t y = coordinates[2 * i + 1];   // Column index
        output[i] = matrix[x * cols + y];    // Access the matrix
    }
}
