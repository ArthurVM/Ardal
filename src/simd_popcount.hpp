// for the the BitMatrix alleleCooccurrence function, which is not currently in use
#pragma once

#include <cstdint>
#include <immintrin.h>

// basic scalar fallback
inline int popcount_scalar(uint64_t x) {
    return __builtin_popcountll(x);
}

#if defined(__AVX512VL__) && defined(__AVX512VPOPCNTDQ__)
inline int popcount_avx512(const uint64_t* a, const uint64_t* b, size_t len,
                           uint64_t* and_out, uint64_t* or_out) {
    int inter = 0, uni = 0;
    for (size_t i = 0; i < len; i += 8) {
        __m512i va = _mm512_loadu_si512((__m512i*)(a + i));
        __m512i vb = _mm512_loadu_si512((__m512i*)(b + i));

        __m512i vand = _mm512_and_si512(va, vb);
        __m512i vor  = _mm512_or_si512(va, vb);

        inter += _mm512_reduce_add_epi64(_mm512_popcnt_epi64(vand));
        uni   += _mm512_reduce_add_epi64(_mm512_popcnt_epi64(vor));
    }
    return {inter, uni};
}
#endif

inline void popcount_pairwise(const std::vector<uint64_t>& A,
                              const std::vector<uint64_t>& B,
                              size_t valid_tail_bits,
                              bool use_simd,
                              size_t& intersection_size,
                              size_t& union_size) {
    size_t n = A.size();
    size_t i = 0;

    // SIMD unrolled
    if (use_simd) {
        for (; i + 3 < n; i += 4) {
            uint64_t a0 = A[i],     b0 = B[i];
            uint64_t a1 = A[i+1],   b1 = B[i+1];
            uint64_t a2 = A[i+2],   b2 = B[i+2];
            uint64_t a3 = A[i+3],   b3 = B[i+3];

            intersection_size += __builtin_popcountll(a0 & b0);
            union_size        += __builtin_popcountll(a0 | b0);

            intersection_size += __builtin_popcountll(a1 & b1);
            union_size        += __builtin_popcountll(a1 | b1);

            intersection_size += __builtin_popcountll(a2 & b2);
            union_size        += __builtin_popcountll(a2 | b2);

            intersection_size += __builtin_popcountll(a3 & b3);
            union_size        += __builtin_popcountll(a3 | b3);
        }
    }

    // scalar tail + final word masking
    for (; i < n; ++i) {
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
