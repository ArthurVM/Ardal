// ardal/detail/flat_matrix.hpp
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <pybind11/numpy.h>
#include "utils/bitops.hpp"


namespace ardal::detail {


struct FlatMatrix {
    const uint64_t* base = nullptr;  // row-major words
    std::size_t n_rows = 0;
    std::size_t wpr = 0;             // words per row
    std::size_t n_cols_bits = 0;
};


/****************************************************************************************************
 * ardal::detail::build_flat_from_packed
 *
 * Construct a zero-copy FlatMatrix view over a pre-packed <u64> NumPy array (e.g., a memmap).
 * Validates that the input is a 2D, tightly row-major array of shape [n_rows, wpr], where
 * wpr == ceil(n_cols_bits / 64). Pins the lifetime of the underlying buffer by returning
 * a borrowed owner handle.
 *
 * INPUT:
 *   arr          (pybind11::array_t<uint64_t, pybind11::array::c_style>):
 *                2D C-contiguous array of uint64 words, shape [n_rows, wpr].
 *   n_cols_bits  (std::size_t): Logical number of bit-columns (before 64-bit packing).
 *
 * OUTPUT:
 *   out_view   (ardal::detail::FlatMatrix&): Filled with base pointer, n_rows, wpr, n_cols_bits.
 *   out_owner  (pybind11::array&): Receives a reference to `arr` to keep its buffer alive.
 *
 * EXCEPTIONS:
 *   std::runtime_error if:
 *     - `arr` is not 2D,
 *     - dtype/layout is not tightly row-major <u64>,
 *     - words-per-row does not equal ceil(n_cols_bits/64).
 *
 * NOTES:
 *   - Performs no copies and does not modify the input; tail bits are not zeroed here.
 *   - Callers must mask the final word per row with ardal::tail_mask(n_cols_bits) when needed.
 ****************************************************************************************************/
inline void build_flat_from_packed(pybind11::array_t<uint64_t, pybind11::array::c_style> arr,
                                   std::size_t n_cols_bits,
                                   FlatMatrix& out_view,
                                   pybind11::array& out_owner)
{
    namespace py = pybind11;
    auto buf = arr.request();
    if (buf.ndim != 2) throw std::runtime_error("packed array must be 2D <u64>");

    const std::size_t n_rows = static_cast<std::size_t>(buf.shape[0]);
    const std::size_t wpr    = static_cast<std::size_t>(buf.shape[1]);
    const std::size_t expect = (n_cols_bits + 63u) / 64u;
    if (wpr != expect) throw std::runtime_error("words_per_row mismatch for packed input");

    // Tight C row-major: strides [ row = wpr*8, col = 8 ]
    if (buf.strides[1] != static_cast<py::ssize_t>(sizeof(uint64_t)) ||
        buf.strides[0] != static_cast<py::ssize_t>(wpr * sizeof(uint64_t))) {
        throw std::runtime_error("packed array must be tightly row-major <u64>");
    }

    out_view.base        = static_cast<const uint64_t*>(buf.ptr);
    out_view.n_rows      = n_rows;
    out_view.wpr         = wpr;
    out_view.n_cols_bits = n_cols_bits;

    out_owner = std::move(arr); // pin lifetime (zero-copy)
}


/****************************************************************************************************
 * ardal::detail::pack_dense_into_storage
 *
 * Pack a dense (bool/int) 2D NumPy array into 64-bit words and expose a FlatMatrix view
 * over the contiguous storage (row-major, words-per-row = ceil(n_cols_bits / 64)).
 * Respects arbitrary NumPy strides; supports itemsize ∈ {1,2,4,8}. Masks the final word
 * of each row to clear bits beyond n_cols_bits.
 *
 * INPUT:
 *   dense        (pybind11::array): 2D array of bool or integer type (1/2/4/8-byte items).
 *   n_cols_bits  (std::size_t): Logical number of bit-columns to pack.
 *
 * OUTPUT:
 *   storage_out  (std::vector<uint64_t>&): Filled to size n_rows * wpr with packed words.
 *   out_view     (ardal::detail::FlatMatrix&): View into `storage_out` (base, n_rows, wpr, n_cols_bits).
 *
 * EXCEPTIONS:
 *   std::runtime_error if:
 *     - `dense` is not 2D,
 *     - `dense.shape[1] != n_cols_bits`,
 *     - dtype is not bool or an integer with itemsize ∈ {1,2,4,8}.
 *
 * NOTES:
 *   - Performs a single allocation/assign into `storage_out`; no aliasing with `dense`.
 *   - Does not assume C-contiguity; reads via strides and width-based loads.
 ****************************************************************************************************/
inline void pack_dense_into_storage(pybind11::array dense,
                                    std::size_t n_cols_bits,
                                    std::vector<uint64_t>& storage_out,
                                    FlatMatrix& out_view)
{
    namespace py = pybind11;
    auto buf = dense.request();
    if (buf.ndim != 2) throw std::runtime_error("dense array must be 2D");

    const std::size_t n_rows = static_cast<std::size_t>(buf.shape[0]);
    const std::size_t n_cols = static_cast<std::size_t>(buf.shape[1]);
    if (n_cols != n_cols_bits) throw std::runtime_error("n_cols_bits mismatch");

    const std::size_t wpr = (n_cols_bits + 63u) / 64u;
    storage_out.assign(n_rows * wpr, 0ULL);

    const ptrdiff_t rs = buf.strides[0];              // bytes between rows
    const ptrdiff_t cs = buf.strides[1];              // bytes between columns
    const char* base_in = static_cast<const char*>(buf.ptr);

    // Safe width-based reader (bool or integer 1/2/4/8 bytes)
    auto read_nz = [&](const char* p) -> bool {
        switch (buf.itemsize) {
            case 1: { uint8_t  v; std::memcpy(&v, p, 1); return v != 0; }
            case 2: { uint16_t v; std::memcpy(&v, p, 2); return v != 0; }
            case 4: { uint32_t v; std::memcpy(&v, p, 4); return v != 0; }
            case 8: { uint64_t v; std::memcpy(&v, p, 8); return v != 0; }
            default: throw std::runtime_error("dense dtype must be bool/int (1/2/4/8 bytes)");
        }
    };

    const bool has_tail = (n_cols_bits & 63u) != 0;
    const uint64_t tailmask = has_tail ? ardal::tail_mask(n_cols_bits) : ~0ULL;

    for (ptrdiff_t r = 0; r < static_cast<ptrdiff_t>(n_rows); ++r) {
        const char* row0 = base_in + r * rs;
        uint64_t* out_row = storage_out.data() + static_cast<std::size_t>(r) * wpr;

        for (std::size_t w = 0; w < wpr; ++w) {
            const std::size_t bit0 = w * 64;
            const std::size_t lim  = std::min<std::size_t>(64, n_cols_bits - bit0);
            uint64_t word = 0ULL;

            const char* cell0 = row0 + static_cast<ptrdiff_t>(bit0) * cs;
            for (std::size_t b = 0; b < lim; ++b)
                if (read_nz(cell0 + static_cast<ptrdiff_t>(b) * cs))
                    word |= (1ULL << b);

            if (has_tail && (w + 1 == wpr)) word &= tailmask;
            out_row[w] = word;
        }
    }

    out_view.base        = storage_out.data();
    out_view.n_rows      = n_rows;
    out_view.wpr         = wpr;
    out_view.n_cols_bits = n_cols_bits;
}


/****************************************************************************************************
 * ardal::detail::compute_masses
 *
 * Compute per-row masses (row sums), per-column masses (column sums), and the total mass
 * (number of set bits) from a FlatMatrix. Applies a tail mask to the final word in each row
 * so that bits beyond n_cols_bits are ignored. Does not modify the input.
 *
 * INPUT:
 *   fm          (const ardal::detail::FlatMatrix&): Flat, row-major bit-packed matrix view.
 *
 * OUTPUT:
 *   row_masses  (std::vector<uint32_t>&): Resized to fm.n_rows; row_masses[r] = #1-bits in row r.
 *   col_masses  (std::vector<uint32_t>&): Resized to fm.n_cols_bits; per-column #1-bits.
 *   total_mass  (uint64_t&): Total number of set bits across the matrix.
 *
 * EXCEPTIONS:
 *   None thrown on valid input; vectors are resized/assigned as needed.
 *
 * COMPLEXITY:
 *   O(fm.n_rows * fm.wpr + nnz), using popcount and x &= (x-1) to iterate set bits.
 ****************************************************************************************************/
inline void compute_masses(const FlatMatrix& fm,
                    std::vector<int>& row_masses,
                    std::vector<int>& col_masses,
                    int64_t& total_mass)
{
    const std::size_t n_rows      = fm.n_rows;
    const std::size_t wpr         = fm.wpr;
    const std::size_t n_cols_bits = fm.n_cols_bits;

    row_masses.assign(n_rows, 0);
    col_masses.assign(n_cols_bits, 0);
    total_mass = 0;

    if (n_rows == 0 || wpr == 0) return;

    const bool has_tail     = (n_cols_bits & 63u) != 0u;
    const uint64_t tailmask = has_tail ? ardal::tail_mask(n_cols_bits) : ~0ULL;

    for (std::size_t r = 0; r < n_rows; ++r) {
        const uint64_t* row = fm.base + r * wpr;
        int row_sum = 0;

        for (std::size_t w = 0; w < wpr; ++w) {
            uint64_t x = row[w];
            if (has_tail && (w + 1 == wpr)) x &= tailmask;

            row_sum += static_cast<int>(ardal::popcnt64(x));

            while (x) {
                const unsigned b = ardal::ctz64_nonzero(x);
                const std::size_t col = (w << 6) | b;
                if (col < n_cols_bits) ++col_masses[col];
                x &= (x - 1);  // clear lowest set bit
            }
        }

        row_masses[r] = row_sum;
        total_mass   += static_cast<int64_t>(row_sum);
    }
}

} // namespace ardal::detail