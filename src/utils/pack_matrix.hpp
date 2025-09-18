// pack_matrix.hpp
#pragma once
#include <cstdint>
#include <vector>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace py = pybind11;

inline void pack_dense_row_to_words(const uint8_t* row, size_t n_cols, uint64_t* out_words) {
    const size_t words = (n_cols + 63) / 64;
    for (size_t w = 0; w < words; ++w) {
        uint64_t word = 0;
        const size_t base = w * 64;
        const size_t limit = (base + 64 <= n_cols) ? 64 : (n_cols - base);
        for (size_t b = 0; b < limit; ++b) {
            // non-zero -> set bit
            word |= (static_cast<uint64_t>(row[base + b] != 0) << b);
        }
        out_words[w] = word;
    }
}

inline py::array pack_dense_to_words(py::array dense) {
    // accept bool/uint8; ensure C-contiguous
    py::buffer_info info = dense.request();
    if (info.ndim != 2) throw std::runtime_error("dense must be 2D");
    if (!(info.format == py::format_descriptor<uint8_t>::format() ||
          info.format == py::format_descriptor<bool>::format())) {
        throw std::runtime_error("dense dtype must be uint8 or bool");
    }
    const size_t n_rows = info.shape[0];
    const size_t n_cols = info.shape[1];
    const size_t words  = (n_cols + 63) / 64;

    // output <u8 row-major
    py::array out(py::dtype("<u8"), { (py::ssize_t)n_rows, (py::ssize_t)words });
    auto out_info = out.request();

    const uint8_t* in_ptr = static_cast<const uint8_t*>(info.ptr);
    uint64_t* out_ptr = static_cast<uint64_t*>(out_info.ptr);
    const size_t in_row_stride  = (size_t)info.strides[0];
    const size_t out_row_stride = (size_t)out_info.strides[0];

    #pragma omp parallel for if(n_rows > 64)
    for (py::ssize_t r = 0; r < (py::ssize_t)n_rows; ++r) {
        const uint8_t* row = reinterpret_cast<const uint8_t*>(
            reinterpret_cast<const char*>(in_ptr) + r * in_row_stride);
        uint64_t* out_row = reinterpret_cast<uint64_t*>(
            reinterpret_cast<char*>(out_ptr) + r * out_row_stride);
        pack_dense_row_to_words(row, n_cols, out_row);
    }
    return out;
}