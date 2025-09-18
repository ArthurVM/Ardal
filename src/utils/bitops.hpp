// ardal/utils/bitops.hpp
#pragma once
#include <cstdint>
#include <cstddef>

namespace ardal {

// --- popcount ---
#if __cpp_lib_bitops >= 201907L
  #include <bit>
  inline uint32_t popcnt64(uint64_t x) noexcept { return std::popcount(x); }
#elif defined(_MSC_VER)
  #include <intrin.h>
  inline uint32_t popcnt64(uint64_t x) noexcept { return static_cast<uint32_t>(__popcnt64(x)); }
#else
  inline uint32_t popcnt64(uint64_t x) noexcept { return static_cast<uint32_t>(__builtin_popcountll(x)); }
#endif

// --- ctz with nonzero precondition (never call with x==0) ---
#if defined(_MSC_VER)
  #include <intrin.h>
  inline unsigned ctz64_nonzero(uint64_t x) noexcept {
    unsigned long idx;
    // BMI1 tzcnt if available; falls back to bsf semantics (well-defined for x!=0)
    return _BitScanForward64(&idx, x) ? static_cast<unsigned>(idx) : 0u;
  }
#else
  inline unsigned ctz64_nonzero(uint64_t x) noexcept {
    return static_cast<unsigned>(__builtin_ctzll(x)); // precondition: x!=0
  }
#endif

// --- iterate set bits (O(popcount)) ---
template <class F>
inline void for_each_set_bit(uint64_t x, F&& f) noexcept {
  while (x) {
    const unsigned b = ctz64_nonzero(x);
    f(b);
    x &= (x - 1); // clear lowest set bit
  }
}

// --- mask for final word (safe when n_cols_bits % 64 == 0) ---
constexpr inline uint64_t tail_mask(std::size_t n_cols_bits) noexcept {
  const std::size_t t = n_cols_bits & 63;
  return t ? ((1ULL << t) - 1ULL) : ~0ULL;
}

} // namespace ardal
