// simd_utils.hpp
/*
Utilities to support bit popcounts: SIMD and scalar implementations.

Copyright 2025 Arthur V. Morris
*/
#pragma once

#include <cstdint>
#include <vector>
#include <iostream>
#include "utils/bitops.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifdef __AVX512VPOPCNTDQ__
#include <immintrin.h>
#endif

namespace ardal {
namespace simd_utils {



// HAMMING

inline uint32_t hamming_scalar( const uint64_t* __restrict r1,
                                const uint64_t* __restrict r2,
                                size_t words_per_row,
                                size_t n_cols_bits ) noexcept {
    if (!words_per_row) return 0;
    uint32_t d = 0;
    const size_t last = words_per_row - 1;

    for (size_t k = 0; k < last; ++k) {
        d += popcnt64(r1[k] ^ r2[k]);
    }
    
    const uint64_t m = ardal::tail_mask(n_cols_bits);
    d += popcnt64((r1[last] ^ r2[last]) & m);
    return d;
}


#ifdef __AVX2__
// popcount 256 bits using PSHUFB nibble table; returns sum of 32 bytes
static inline __m256i popcnt256_bytewise( __m256i v ) {
    // 4-bit popcount LUT in 16 bytes
    const __m128i lut4 = _mm_setr_epi8(
        0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4
    );
    const __m256i lut  = _mm256_broadcastsi128_si256(lut4);

    const __m256i lo_nibble = _mm256_and_si256(v, _mm256_set1_epi8(0x0f));
    const __m256i hi_nibble = _mm256_and_si256(_mm256_srli_epi16(v, 4), _mm256_set1_epi8(0x0f));

    const __m256i cnt_lo = _mm256_shuffle_epi8(lut, lo_nibble);
    const __m256i cnt_hi = _mm256_shuffle_epi8(lut, hi_nibble);

    return _mm256_add_epi8(cnt_lo, cnt_hi); // 32 x u8 partial popcounts
}

// horizontally sum 32 unsigned bytes in a 256-bit register to a 32-bit integer
static inline uint32_t hsum256_epu8( __m256i v8 ) {
    // widen to u16 and sum
    const __m256i zero = _mm256_setzero_si256();
    __m256i v16 = _mm256_sad_epu8(v8, zero); // sums per 16 bytes into 4x u64 lanes
    // add the two 128-bit halves
    __m128i lo = _mm256_castsi256_si128(v16);
    __m128i hi = _mm256_extracti128_si256(v16, 1);
    __m128i sum128 = _mm_add_epi64(lo, hi);
    uint64_t s = (uint64_t)_mm_cvtsi128_si64(sum128) + (uint64_t)_mm_extract_epi64(sum128, 1);
    return (uint32_t)s;
}

inline uint32_t hamming_avx2( const uint64_t* __restrict row1,
                              const uint64_t* __restrict row2,
                              size_t words_per_row,
                              size_t n_cols_bits ) {
    uint32_t dist = 0;
    size_t k = 0;

    // process 4x u64 = 256 bits per iteration
    const size_t vec_step = 4;
    const size_t last_vec = (words_per_row / vec_step) * vec_step;

    for (; k < last_vec; k += vec_step) {
        __m256i a = _mm256_loadu_si256((const __m256i*)(row1 + k));
        __m256i b = _mm256_loadu_si256((const __m256i*)(row2 + k));
        __m256i x = _mm256_xor_si256(a, b);
        __m256i c = popcnt256_bytewise(x);
        dist += hsum256_epu8(c);
    }

    // tail words (and apply bit tail mask on the final word)
    for (; k < words_per_row; ++k) {
        uint64_t x = row1[k] ^ row2[k];
        if (k + 1 == words_per_row) x &= tail_mask(n_cols_bits);
        dist += __builtin_popcountll(x);
    }
    return dist;
}
#endif


#ifdef __AVX512VPOPCNTDQ__
// mad 8×64-bit popcounts per instruction for AVX512
// NOTE: CURRENTLY UNTESTED
inline uint32_t hamming_avx512( const uint64_t* __restrict row1,
                                const uint64_t* __restrict row2,
                                size_t words_per_row,
                                size_t n_cols_bits ) {
    uint32_t dist = 0;
    size_t k = 0;
    const size_t vec_step = 8; // 8x u64 per 512-bit register
    const size_t last_vec = (words_per_row / vec_step) * vec_step;

    for (; k < last_vec; k += vec_step) {
        __m512i a = _mm512_loadu_si512((const void*)(row1 + k));
        __m512i b = _mm512_loadu_si512((const void*)(row2 + k));
        __m512i x = _mm512_xor_si512(a, b);
        __m512i pc = _mm512_popcnt_epi64(x);               // 8 lanes of 64-bit popcounts
        // horizontal sum of 8 u64 lanes
        dist += (uint32_t)_mm512_reduce_add_epi64(pc);
    }

    for (; k < words_per_row; ++k) {
        uint64_t x = row1[k] ^ row2[k];
        if (k + 1 == words_per_row) x &= tail_mask(n_cols_bits);
        dist += __builtin_popcountll(x);
    }
    return dist;
}
#endif



// EPSILON NEIGHBOURHOOD HAMMING

inline int epsilon_hamming_scalar( const uint64_t* __restrict row1,
                                   const uint64_t* __restrict row2,
                                   size_t words_per_row,
                                   size_t n_cols_bits,
                                   int epsilon ) noexcept {
    if (words_per_row == 0) return 0;
    int dist = 0;
    const size_t last = words_per_row - 1;
    for (size_t k = 0; k < last; ++k) {
        dist += popcnt64(row1[k] ^ row2[k]);
        if (dist > epsilon) return dist;
    }
    const uint64_t m = tail_mask(n_cols_bits);
    dist += popcnt64((row1[last] ^ row2[last]) & m);
    return dist;
}


#ifdef __AVX2__
// 256-bit bytewise popcount via PSHUFB (returns 32 u8 lane counts)
static inline __m256i _ardal_popcnt256_u8( __m256i v ) noexcept {
    const __m128i lut4 = _mm_setr_epi8(0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);
    const __m256i lut  = _mm256_broadcastsi128_si256(lut4);
    const __m256i lo   = _mm256_and_si256(v, _mm256_set1_epi8(0x0f));
    const __m256i hi   = _mm256_and_si256(_mm256_srli_epi16(v, 4), _mm256_set1_epi8(0x0f));
    return _mm256_add_epi8(_mm256_shuffle_epi8(lut, lo),
                           _mm256_shuffle_epi8(lut, hi));
}

// horizontal sum of 32 u8 lane counts to uint32
static inline uint32_t _ardal_hsum256_u8( __m256i v8 ) noexcept {
    const __m256i zero = _mm256_setzero_si256();
    __m256i v16 = _mm256_sad_epu8(v8, zero); // 4x u64 partial sums
    __m128i lo = _mm256_castsi256_si128(v16);
    __m128i hi = _mm256_extracti128_si256(v16, 1);
    __m128i sum128 = _mm_add_epi64(lo, hi);
    uint64_t s = (uint64_t)_mm_cvtsi128_si64(sum128) + (uint64_t)_mm_extract_epi64(sum128, 1);
    return (uint32_t)s;
}

inline int epsilon_hamming_avx2( const uint64_t* __restrict row1,
                                 const uint64_t* __restrict row2,
                                 size_t words_per_row,
                                 size_t n_cols_bits,
                                 int epsilon) noexcept {
    int dist = 0;
    size_t k = 0;
    const size_t vec_step = 4; // 4x u64 per 256-bit
    const size_t last_vec = (words_per_row / vec_step) * vec_step;

    for (; k < last_vec; k += vec_step) {
        __m256i a = _mm256_loadu_si256((const __m256i*)(row1 + k));
        __m256i b = _mm256_loadu_si256((const __m256i*)(row2 + k));
        __m256i x = _mm256_xor_si256(a, b);
        dist += (int)_ardal_hsum256_u8(_ardal_popcnt256_u8(x));
        if (dist > epsilon) return dist;
    }
    // scalar tail, with tail-bit mask on the very last word
    for (; k < words_per_row; ++k) {
        uint64_t x = row1[k] ^ row2[k];
        if (k + 1 == words_per_row) x &= tail_mask(n_cols_bits);
        dist += popcnt64(x);
        if (dist > epsilon) return dist;
    }
    return dist;
}
#endif



// INNER PRODUCT

inline int inner_product_scalar( const uint64_t* __restrict row1,
                                 const uint64_t* __restrict row2,
                                 size_t words_per_row,
                                 size_t n_cols_bits) noexcept {
    if (words_per_row == 0) return 0;
    int ip = 0;
    const size_t last = words_per_row - 1;
    for (size_t k = 0; k < last; ++k)
        ip += popcnt64(row1[k] & row2[k]);
    const uint64_t m = tail_mask(n_cols_bits);
    ip += popcnt64((row1[last] & row2[last]) & m);
    return ip;
}


#ifdef __AVX2__
inline int inner_product_avx2( const uint64_t* __restrict row1,
                               const uint64_t* __restrict row2,
                               size_t words_per_row,
                               size_t n_cols_bits) noexcept {
    int ip = 0;
    size_t k = 0;
    const size_t vec_step = 4;
    const size_t last_vec = (words_per_row / vec_step) * vec_step;

    for (; k < last_vec; k += vec_step) {
        __m256i a = _mm256_loadu_si256((const __m256i*)(row1 + k));
        __m256i b = _mm256_loadu_si256((const __m256i*)(row2 + k));
        __m256i x = _mm256_and_si256(a, b);
        ip += (int)_ardal_hsum256_u8(_ardal_popcnt256_u8(x));
    }
    for (; k < words_per_row; ++k) {
        uint64_t x = row1[k] & row2[k];
        if (k + 1 == words_per_row) x &= tail_mask(n_cols_bits);
        ip += popcnt64(x);
    }
    return ip;
}
#endif



// JACCARD

inline void jaccard_popcount( const std::vector<uint64_t>& A,
                              const std::vector<uint64_t>& B,
                              size_t n_cols_bits,
                              size_t& intersection_size,
                              size_t& union_size) noexcept {
    const size_t n = A.size();
    intersection_size = 0;
    union_size = 0;
    if (n == 0) return;

    const uint64_t m = tail_mask(n_cols_bits);
    const size_t last = n - 1;

    for (size_t i = 0; i < last; ++i) {
        uint64_t a = A[i], b = B[i];
        intersection_size += popcnt64(a & b);
        union_size        += popcnt64(a | b);
    }
    // last word masked
    {
        uint64_t a = A[last] & m;
        uint64_t b = B[last] & m;
        intersection_size += popcnt64(a & b);
        union_size        += popcnt64(a | b);
    }
}

inline void jaccard_popcount_ptr( const uint64_t* __restrict A,
                                  const uint64_t* __restrict B,
                                  size_t words_per_row,
                                  size_t n_cols_bits,
                                  size_t& intersection_size,
                                  size_t& union_size ) noexcept {
    intersection_size = 0;
    union_size = 0;
    if (words_per_row == 0) return;

    const size_t last = words_per_row - 1;
    for (size_t i = 0; i < last; ++i) {
        intersection_size += popcnt64(A[i] & B[i]);
        union_size        += popcnt64(A[i] | B[i]);
    }
    const uint64_t m = tail_mask(n_cols_bits);
    intersection_size += popcnt64((A[last] & B[last]) & m);
    union_size        += popcnt64((A[last] | B[last]) & m);
}


} // namespace simd_utils
} // namespace ardal
