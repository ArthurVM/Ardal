#include "gtest/gtest.h"
#include "AlleleMatrix.hpp" // Include your header file
#include <pybind11/numpy.h>


TEST(AlleleMatrixTest, Initialization) {
    py::array_t<uint8_t> input_matrix({2, 3}, {3, 2}, "input data description"); // Example input
    _ardal::AlleleMatrix matrix(input_matrix);  // Assuming _ardal is your namespace
    // Add assertions to check initialization, e.g., dimensions
    ASSERT_EQ(matrix.get_matrix().shape(0), 2); // Assumes existence of a shape method.
    ASSERT_EQ(matrix.get_matrix().shape(1), 3);
}


TEST(AlleleMatrixTest, HammingZeroDistance) {
    // Test case where Hamming distance should be zero (identical arrays)
    py::array_t<uint8_t> data({2, 3}, {6}, {0, 1, 0, 0, 1, 0});
    _ardal::AlleleMatrix matrix1(data);
    _ardal::AlleleMatrix matrix2(data);

    py::array_t<int> result = matrix1.hamming(matrix2); // or similar call depending on how you implemented it
    ASSERT_EQ(result.at(0), 0); // Check first distance
    ASSERT_EQ(result.at(1), 0);

}


TEST(AlleleMatrixTest, HammingNonZeroDistance) {
    py::array_t<uint8_t> data1({2, 3}, {6}, {0, 1, 0, 0, 1, 0});
    _ardal::AlleleMatrix matrix1(data1);
    py::array_t<uint8_t> data2({2, 3}, {6}, {1, 1, 0, 0, 0, 0}); // Different data
    _ardal::AlleleMatrix matrix2(data2);


    py::array_t<int> result = matrix1.hamming(matrix2); 
    ASSERT_EQ(result.at(0), 1); // Expect a non-zero distance
    ASSERT_EQ(result.at(1), 1);

}


// Add more test cases for different scenarios and edge cases:
// - Empty input matrix
// - Matrices of different dimensions (if allowed by your implementation)
// - Different data types (if applicable)
// - Test other methods like jaccard, etc.


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
