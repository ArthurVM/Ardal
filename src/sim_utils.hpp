#pragma once

#include <cstdint>
#include <vector>
#include <immintrin.h>


namespace _ardal {
namespace simd_utils {

// HAMMING

inline uint32_t hamming_scalar( const uint64_t* row1,
                                const uint64_t* row2,
                                size_t packed_cols ) {
    uint32_t dist = 0;
    for (size_t k = 0; k < packed_cols; ++k) {
        dist += __builtin_popcountll(row1[k] ^ row2[k]);
    }
    return dist;
}


inline uint32_t hamming_avx2(const uint64_t* row1, const uint64_t* row2, size_t packed_cols) {
    uint32_t dist = 0;
    size_t k = 0;
    for (; k + 4 <= packed_cols; k += 4) {
        __m256i a = _mm256_loadu_si256((__m256i const*)(row1 + k));
        __m256i b = _mm256_loadu_si256((__m256i const*)(row2 + k));
        __m256i xor_result = _mm256_xor_si256(a, b);
        uint64_t chunks[4];
        _mm256_storeu_si256((__m256i*)chunks, xor_result);
        dist += __builtin_popcountll(chunks[0]) + __builtin_popcountll(chunks[1]) +
                __builtin_popcountll(chunks[2]) + __builtin_popcountll(chunks[3]);
    }
    for (; k < packed_cols; ++k) {
        dist += __builtin_popcountll(row1[k] ^ row2[k]);
    }
    return dist;
}



// EPSILON NEIGHBOURHOOD HAMMING

inline int epsilon_hamming_scalar( const uint64_t* row1,
                                   const uint64_t* row2,
                                   size_t packed_cols,
                                   int epsilon ) {
    int dist = 0;
    for (size_t k = 0; k < packed_cols; ++k) {
        dist += __builtin_popcountll(row1[k] ^ row2[k]);
        if (dist > epsilon) return dist;
    }
    return dist;
}


inline int epsilon_hamming_avx2( const uint64_t* row1,
                                 const uint64_t* row2,
                                 size_t packed_cols,
                                 int epsilon ) {
    int dist = 0;
    size_t k = 0;
    for (; k + 4 <= packed_cols; k += 4) {
        __m256i a = _mm256_loadu_si256((__m256i const*)(row1 + k));
        __m256i b = _mm256_loadu_si256((__m256i const*)(row2 + k));
        __m256i xor_result = _mm256_xor_si256(a, b);
        uint64_t chunks[4];
        _mm256_storeu_si256((__m256i*)chunks, xor_result);
        dist += __builtin_popcountll(chunks[0]) + __builtin_popcountll(chunks[1]) +
                __builtin_popcountll(chunks[2]) + __builtin_popcountll(chunks[3]);
        if (dist > epsilon) return dist;
    }
    for (; k < packed_cols; ++k) {
        dist += __builtin_popcountll(row1[k] ^ row2[k]);
        if (dist > epsilon) return dist;
    }
    return dist;
}



// INNER PRODUCT

inline int inner_product_scalar(const uint64_t* row1, const uint64_t* row2, size_t packed_cols) {
    int ip = 0;
    for (size_t k = 0; k < packed_cols; ++k) {
        ip += __builtin_popcountll(row1[k] & row2[k]);
    }
    return ip;
}


inline int inner_product_avx2(const uint64_t* row1, const uint64_t* row2, size_t packed_cols) {
    int ip = 0;
    size_t k = 0;
    for (; k + 4 <= packed_cols; k += 4) {
        __m256i a = _mm256_loadu_si256((__m256i const*)(row1 + k));
        __m256i b = _mm256_loadu_si256((__m256i const*)(row2 + k));
        __m256i and_result = _mm256_and_si256(a, b);
        uint64_t chunks[4];
        _mm256_storeu_si256((__m256i*)chunks, and_result);
        ip += __builtin_popcountll(chunks[0]) + __builtin_popcountll(chunks[1]) +
              __builtin_popcountll(chunks[2]) + __builtin_popcountll(chunks[3]);
    }
    for (; k < packed_cols; ++k) {
        ip += __builtin_popcountll(row1[k] & row2[k]);
    }
    return ip;
}



// JACCARD

inline void jaccard_popcount(const std::vector<uint64_t>& A,
                             const std::vector<uint64_t>& B,
                             size_t valid_tail_bits,
                             size_t& intersection_size,
                             size_t& union_size) {
    size_t n = A.size();
    intersection_size = 0;
    union_size = 0;

    for (size_t i = 0; i < n; ++i) {
        uint64_t a = A[i];
        uint64_t b = B[i];
        if (i == n - 1 && valid_tail_bits > 0) {
            uint64_t mask = (1ULL << valid_tail_bits) - 1;
            a &= mask;
            b &= mask;
        }
        intersection_size += __builtin_popcountll(a & b);
        union_size        += __builtin_popcountll(a | b);
    }
}

} // namespace simd_utils
} // namespace _ardal

