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
#include "utils/PythonLogger.hpp"


namespace ardal {

// AVX2 kernels
#if defined(ARDAL_HAVE_AVX2_KERNELS) && ARDAL_HAVE_AVX2_KERNELS
    #include <immintrin.h>
    // [[gnu::target("avx2")]]
    inline uint32_t hamming_avx2(const uint64_t* __restrict, const uint64_t* __restrict, size_t, size_t) noexcept;
    inline uint32_t epsilon_hamming_avx2(const uint64_t* __restrict, const uint64_t* __restrict, size_t, size_t, uint32_t) noexcept;
    inline uint32_t inner_product_avx2(const uint64_t* __restrict, const uint64_t* __restrict, size_t, size_t) noexcept;
#endif

// AVX-512 VPOPCNTDQ kernels
#if defined(ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS) && ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS
    #include <immintrin.h>
    // [[gnu::target("avx512vpopcntdq")]]
    inline uint32_t hamming_avx512(const uint64_t* __restrict, const uint64_t* __restrict, size_t, size_t) noexcept;
    inline uint32_t epsilon_hamming_avx512(const uint64_t* __restrict, const uint64_t* __restrict, size_t, size_t, uint32_t) noexcept;
    inline uint32_t inner_product_avx512(const uint64_t* __restrict, const uint64_t* __restrict, size_t, size_t) noexcept;
#endif

// HELPER FUNCTIONS
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
    const __m256i zero = _mm256_setzero_si256();
    __m256i v16 = _mm256_sad_epu8(v8, zero); // sums per 16 bytes into 4x u64 lanes
    // add the two 128-bit halves
    __m128i lo = _mm256_castsi256_si128(v16);
    __m128i hi = _mm256_extracti128_si256(v16, 1);
    __m128i sum128 = _mm_add_epi64(lo, hi);
    uint64_t s = (uint64_t)_mm_cvtsi128_si64(sum128) + (uint64_t)_mm_extract_epi64(sum128, 1);
    return (uint32_t)s;
}



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


inline uint32_t hamming_avx2( const uint64_t* __restrict row1,
                              const uint64_t* __restrict row2,
                              size_t words_per_row,
                              size_t n_cols_bits ) noexcept {
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
        dist += popcnt64(x);
    }
    return dist;
}


// mad 8×64-bit popcounts per instruction for AVX512
// NOTE: CURRENTLY UNTESTED
inline uint32_t hamming_avx512( const uint64_t* __restrict row1,
                                const uint64_t* __restrict row2,
                                size_t words_per_row,
                                size_t n_cols_bits ) noexcept {
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
        dist += popcnt64(x);
    }
    return dist;
}



// EPSILON NEIGHBOURHOOD HAMMING
inline uint32_t epsilon_hamming_scalar( const uint64_t* __restrict row1,
                                   const uint64_t* __restrict row2,
                                   size_t words_per_row,
                                   size_t n_cols_bits,
                                   uint32_t epsilon ) noexcept {
    if (words_per_row == 0) return 0;
    uint32_t dist = 0;
    const size_t last = words_per_row - 1;
    for (size_t k = 0; k < last; ++k) {
        dist += popcnt64(row1[k] ^ row2[k]);
        if (dist > epsilon) return dist;
    }
    const uint64_t m = tail_mask(n_cols_bits);
    dist += popcnt64((row1[last] ^ row2[last]) & m);
    return dist;
}


inline uint32_t epsilon_hamming_avx2( const uint64_t* __restrict row1,
                                      const uint64_t* __restrict row2,
                                      size_t words_per_row,
                                      size_t n_cols_bits,
                                      uint32_t epsilon ) noexcept {
    uint32_t dist = 0;
    size_t k = 0;
    const size_t vec_step = 4; // 4x u64 per 256-bit
    const size_t last_vec = (words_per_row / vec_step) * vec_step;

    for (; k < last_vec; k += vec_step) {
        __m256i a = _mm256_loadu_si256((const __m256i*)(row1 + k));
        __m256i b = _mm256_loadu_si256((const __m256i*)(row2 + k));
        __m256i x = _mm256_xor_si256(a, b);
        dist += (uint32_t)hsum256_epu8(popcnt256_bytewise(x));
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


// INNER PRODUCT
inline uint32_t inner_product_scalar( const uint64_t* __restrict row1,
                                      const uint64_t* __restrict row2,
                                      size_t words_per_row,
                                      size_t n_cols_bits ) noexcept {
    if (words_per_row == 0) return 0;
    uint32_t ip = 0;
    const size_t last = words_per_row - 1;
    for (size_t k = 0; k < last; ++k)
        ip += popcnt64(row1[k] & row2[k]);
    const uint64_t m = tail_mask(n_cols_bits);
    ip += popcnt64((row1[last] & row2[last]) & m);
    return ip;
}


inline uint32_t inner_product_avx512( const uint64_t* __restrict row1,
                                      const uint64_t* __restrict row2,
                                      size_t words_per_row,
                                      size_t n_cols_bits ) noexcept {
    uint32_t ip = 0;
    size_t k = 0;
    const size_t vec_step = 8;
    const size_t last_vec = (words_per_row / vec_step) * vec_step;

    for (; k < last_vec; k += vec_step) {
        __m512i a = _mm512_loadu_si512((const void*)(row1 + k));
        __m512i b = _mm512_loadu_si512((const void*)(row2 + k));
        __m512i x = _mm512_and_si512(a, b);
        __m512i pc = _mm512_popcnt_epi64(x);   // 8 lanes
        ip += (uint32_t)_mm512_reduce_add_epi64(pc);
    }
    for (; k < words_per_row; ++k) {
        uint64_t x = row1[k] & row2[k];
        if (k + 1 == words_per_row) x &= tail_mask(n_cols_bits);
        ip += popcnt64(x);
    }
    return ip;
}


inline uint32_t inner_product_avx2( const uint64_t* __restrict row1,
                                    const uint64_t* __restrict row2,
                                    size_t words_per_row,
                                    size_t n_cols_bits ) noexcept {
    uint32_t ip = 0;
    size_t k = 0;
    const size_t vec_step = 4;
    const size_t last_vec = (words_per_row / vec_step) * vec_step;

    for (; k < last_vec; k += vec_step) {
        __m256i a = _mm256_loadu_si256((const __m256i*)(row1 + k));
        __m256i b = _mm256_loadu_si256((const __m256i*)(row2 + k));
        __m256i x = _mm256_and_si256(a, b);
        ip += (uint32_t)hsum256_epu8(popcnt256_bytewise(x));
    }
    for (; k < words_per_row; ++k) {
        uint64_t x = row1[k] & row2[k];
        if (k + 1 == words_per_row) x &= tail_mask(n_cols_bits);
        ip += popcnt64(x);
    }
    return ip;
}



// JACCARD
inline void jaccard_popcount( const std::vector<uint64_t>& A,
                              const std::vector<uint64_t>& B,
                              size_t n_cols_bits,
                              size_t& intersection_size,
                              size_t& union_size ) noexcept {
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

} // namespace ardal


// SIMD DISPATCHERS
namespace ardal::simd_dispatchers {

// --- function pointer types ---
using HammingFn      = uint32_t(*)(const uint64_t*, const uint64_t*, size_t, size_t) noexcept;
using HammingEpsFn   = uint32_t(*)(const uint64_t*, const uint64_t*, size_t, size_t, uint32_t) noexcept;
using InnerProductFn = uint32_t(*)(const uint64_t*, const uint64_t*, size_t, size_t) noexcept;

// --- kernel declarations ---
using ardal::hamming_scalar;
using ardal::epsilon_hamming_scalar;
using ardal::inner_product_scalar;

#if ARDAL_HAVE_AVX2_KERNELS
using ardal::hamming_avx2;
using ardal::epsilon_hamming_avx2;
using ardal::inner_product_avx2;
#endif

#if ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS
using ardal::hamming_avx512;
using ardal::epsilon_hamming_avx512;
using ardal::inner_product_avx512;
#endif

// --- runtime CPU feature probes ---
#if defined(__GNUC__) && !defined(__clang__)
inline void cpu_init_once() {
    static bool init = ([](){ __builtin_cpu_init(); return true; })();
    (void)init;
}
#else
inline void cpu_init_once() {}
#endif

inline bool cpu_has_avx2() {
#if defined(__GNUC__) || defined(__clang__)
    cpu_init_once();
    return __builtin_cpu_supports("avx2");
#else
    return false;
#endif
}

inline bool cpu_has_avx512_vpopcnt() {
#if defined(__GNUC__) || defined(__clang__)
    cpu_init_once();
    return __builtin_cpu_supports("avx512vpopcntdq");
#else
    return false;
#endif
}

// --- epsilon-hamming dispatcher ---
inline HammingEpsFn epsilon_hamming_dispatcher(bool use_simd = true) {
    static HammingEpsFn simd_fp = []() -> HammingEpsFn {
#if ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS
        if (cpu_has_avx512_vpopcnt()) {
            ardal::utils::log_info("epsilon_hamming: dispatching to AVX-512 VPOPCNTDQ.");
            return &ardal::epsilon_hamming_avx512;  // fully qualified
        }
#endif
#if ARDAL_HAVE_AVX2_KERNELS
        if (cpu_has_avx2()) {
            ardal::utils::log_info("epsilon_hamming: dispatching to AVX2");
            return &ardal::epsilon_hamming_avx2;    // fully qualified
        }
#endif
        return nullptr;
    }();

    if (use_simd && simd_fp) return simd_fp;
    static HammingEpsFn scalar_fp = &ardal::epsilon_hamming_scalar; // fully qualified
    if (use_simd && !simd_fp) ardal::utils::log_info("epsilon_hamming: dispatching to scalar.");
    return scalar_fp;
}

// --- hamming dispatcher ---
inline HammingFn hamming_dispatcher(bool use_simd = true) {
    static HammingFn simd_fp = []() -> HammingFn {
#if ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS
        if (cpu_has_avx512_vpopcnt()) {
            ardal::utils::log_info("hamming: dispatching to AVX-512 VPOPCNTDQ.");
            return &ardal::hamming_avx512;
        }
#endif
#if ARDAL_HAVE_AVX2_KERNELS
        if (cpu_has_avx2()) {
            ardal::utils::log_info("hamming: dispatching to AVX2");
            return &ardal::hamming_avx2;
        }
#endif
        return nullptr;
    }();

    if (use_simd && simd_fp) return simd_fp;
    static HammingFn scalar_fp = &ardal::hamming_scalar;
    if (use_simd && !simd_fp) ardal::utils::log_info("hamming: dispatching to scalar.");
    return scalar_fp;
}

// --- inner product dispatcher ---
inline InnerProductFn inner_product_dispatcher(bool use_simd = true) {
    static InnerProductFn simd_fp = []() -> InnerProductFn {
#if ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS
        if (cpu_has_avx512_vpopcnt()) {
            ardal::utils::log_info("inner_product: dispatching to AVX-512.");
            return &ardal::inner_product_avx512;
        }
#endif
#if ARDAL_HAVE_AVX2_KERNELS
        if (cpu_has_avx2()) {
            ardal::utils::log_info("inner_product: dispatching to AVX2.");
            return &ardal::inner_product_avx2;
        }
#endif
        return nullptr;
    }();

    if (use_simd && simd_fp) return simd_fp;
    static InnerProductFn scalar_fp = &ardal::inner_product_scalar;
    if (use_simd && !simd_fp) ardal::utils::log_info("inner_product: dispatching to scalar.");
    return scalar_fp;
}


// SIMD diagnostics
inline void simd_diag() {
    std::stringstream ss;
#if defined(ARDAL_HAVE_AVX2_KERNELS)
    ss << "ARDAL_HAVE_AVX2_KERNELS = " << int(ARDAL_HAVE_AVX2_KERNELS);
    ardal::utils::log_debug(static_cast<string>(ss.str()));
#else
    ardal::utils::log_debug("ARDAL_HAVE_AVX2_KERNELS is NOT defined");
#endif
#if defined(ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS)
    std::stringstream ss512;
    ss512 << "ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS = " << int(ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS);
    ardal::utils::log_debug(static_cast<string>(ss.str()));
#else
    ardal::utils::log_debug("ARDAL_HAVE_AVX512VPOPCNTDQ_KERNELS is NOT defined");
#endif
    std::stringstream ss2;
    ss2 << "cpu_has_avx2 = " << int(cpu_has_avx2()) << "; cpu_has_avx512_vpopcnt = " << int(cpu_has_avx512_vpopcnt());
    ardal::utils::log_debug(static_cast<string>(ss2.str()));
}

} // namespace ardal::simd_dispatchers