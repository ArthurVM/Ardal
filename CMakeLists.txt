cmake_minimum_required(VERSION 3.10)
project(ArdalTests)

find_package(GTest REQUIRED)
find_package(pybind11 REQUIRED)

add_executable(test_allele_matrix src/AlleleMatrix.cpp src/AlleleMatrix_test.cpp)
target_link_libraries(test_allele_matrix GTest::gtest_main pybind11::pybind11)