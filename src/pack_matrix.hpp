// pack.hpp
#pragma once
#include <cstdint>
#include <vector>
#include <stdexcept>


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