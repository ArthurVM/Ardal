/*
Copyright 2025 Arthur V. Morris
*/
#include "HybridMatrix.hpp"
#include "BitMatrix.hpp"
#include "RoaringMatrix.hpp"
#include "utils/bitops.hpp"
#include "utils/PythonLogger.hpp"


// helpers for the vector<vector<uint64_t>> logic
// NOTE: this will be worked out in future versions as I port Ardal to using pointer+stride logic


ardal::detail::PackedVV make_vv_from_packed_numpy( py::array_t<uint64_t, py::array::c_style | py::array::forcecast> packed,
                                    std::size_t n_cols_bits) {
    auto buf = packed.request();
    if (buf.ndim != 2) throw std::runtime_error("packed array must be 2D <u64>");

    const std::size_t n_rows = buf.shape[0];
    const std::size_t wpr    = buf.shape[1];
    const std::size_t expect = (n_cols_bits + 63) / 64;

    if (wpr != expect) throw std::runtime_error("words_per_row mismatch");

    if (buf.strides[0] != static_cast<py::ssize_t>(wpr * sizeof(uint64_t)))
        throw std::runtime_error("packed array must be tightly row-major");

    const uint64_t* base = static_cast<const uint64_t*>(buf.ptr);
    const uint64_t mask  = ardal::tail_mask(n_cols_bits);

    ardal::detail::PackedVV out;
    out.vv.assign(n_rows, std::vector<uint64_t>(wpr));
    out.n_rows = n_rows; out.wpr = wpr;

    for (std::size_t r = 0; r < n_rows; ++r) {
        const uint64_t* src = base + r * wpr;
        std::memcpy(out.vv[r].data(), src, wpr * sizeof(uint64_t));
        if ((n_cols_bits & 63) != 0) out.vv[r].back() &= mask;
    }
    return out;
}


std::vector<std::vector<uint64_t>> pack_dense_to_vv( py::array dense,
                                                     std::size_t n_cols_bits ) {
    auto buf = dense.request();
    if (buf.ndim != 2) throw std::runtime_error("dense array must be 2D");

    const std::size_t n_rows = buf.shape[0];
    const std::size_t n_cols = buf.shape[1];

    if (n_cols_bits != n_cols) throw std::runtime_error("n_cols_bits mismatch");

    const std::size_t wpr  = (n_cols + 63) / 64;
    const uint64_t mask    = ardal::tail_mask(n_cols_bits);
    const ptrdiff_t rs = buf.strides[0], cs = buf.strides[1];
    const char kind = buf.format[0];

    auto read_nz = [&](const char* p)->bool {
        switch (kind) {
          case '?': return *reinterpret_cast<const bool*>(p);
          case 'b': return *reinterpret_cast<const int8_t*>(p) != 0;
          case 'B': return *reinterpret_cast<const uint8_t*>(p) != 0;
          case 'h': return *reinterpret_cast<const int16_t*>(p) != 0;
          case 'H': return *reinterpret_cast<const uint16_t*>(p) != 0;
          case 'i': return *reinterpret_cast<const int32_t*>(p) != 0;
          case 'I': return *reinterpret_cast<const uint32_t*>(p) != 0;
          case 'l': return *reinterpret_cast<const long*>(p) != 0;
          case 'L': return *reinterpret_cast<const unsigned long*>(p) != 0;
          case 'q': return *reinterpret_cast<const long long*>(p) != 0;
          case 'Q': return *reinterpret_cast<const unsigned long long*>(p) != 0;
          default: throw std::runtime_error("dense dtype must be bool or integer");
        }
    };

    std::vector<std::vector<uint64_t>> vv(n_rows, std::vector<uint64_t>(wpr, 0ULL));

    for (ptrdiff_t r = 0; r < static_cast<ptrdiff_t>(n_rows); ++r) {
        const char* row0 = static_cast<const char*>(buf.ptr) + r * rs;
        uint64_t* out = vv[static_cast<std::size_t>(r)].data();
        for (std::size_t w = 0; w < wpr; ++w) {
            const std::size_t base = w * 64;
            const std::size_t lim  = std::min<std::size_t>(64, n_cols - base);
            uint64_t word = 0ULL;
            const char* cell = row0 + static_cast<ptrdiff_t>(base) * cs;
            for (std::size_t b = 0; b < lim; ++b)
                if (read_nz(cell + static_cast<ptrdiff_t>(b) * cs)) word |= (1ULL << b);
            out[w] = (w + 1 == wpr && (n_cols_bits & 63)) ? (word & mask) : word;
        }
    }
    return vv;
}

void compute_masses(const std::vector<std::vector<uint64_t>>& vv,
                    std::size_t n_cols_bits,
                    std::vector<int>& row_masses,
                    std::vector<int>& col_masses,
                    int64_t& total_mass)
{
    const std::size_t n_rows = vv.size();
    if (!n_rows) { total_mass = 0; return; }
    const std::size_t wpr = vv[0].size();
    const uint64_t mask   = ardal::tail_mask(n_cols_bits);

    row_masses.assign(n_rows, 0);
    col_masses.assign(n_cols_bits, 0);
    total_mass = 0;

    for (std::size_t r = 0; r < n_rows; ++r) {
        const auto& row = vv[r];
        uint32_t rm = 0;
        for (std::size_t w = 0; w < wpr; ++w) {
            uint64_t x = (w + 1 == wpr && (n_cols_bits & 63)) ? (row[w] & mask) : row[w];
            rm += ardal::popcnt64(x);
            while (x) {
                const uint64_t lsb = x & (~x + 1);
                const unsigned b   = static_cast<unsigned>(__builtin_ctzll(x));
                const std::size_t col = (w << 6) | b;
                if (col < n_cols_bits) ++col_masses[col];
                x ^= lsb;
            }
        }
        row_masses[r] = static_cast<int>(rm);
        total_mass += rm;
    }
}

} // unnamed namespace


namespace ardal {


/****************************************************************************************************
 * ardal::HybridMatrix::HybridMatrix
 *
 * Constructor for the HybridMatrix class.
 *
 * This constructor initializes the HybridMatrix object, which can use either a `BitMatrix`
 * or a `RoaringMatrix` as its backend, depending on the density of the input matrix.
 *
 * INPUT:
 *  packed_matrix (py::array_t<uint8_t>): A 2D 64bit packed binary matrix.
 *  use_roaring_if_sparse (bool): If true, the `RoaringMatrix` backend will be used if the
 *                                matrix density is below `density_threshold`.
 *  density_threshold (double): The density threshold below which `RoaringMatrix` is preferred.
 *
 * OUTPUT:
 *  None (constructor)
 *
 * EXCEPTIONS:
 *  Propagates exceptions from `BitMatrix` and `RoaringMatrix` constructors.
 ****************************************************************************************************/  
HybridMatrix::HybridMatrix( py::array packed_or_dense_matrix,
                            bool is_bitpacked,
                            size_t n_cols_bits,
                            bool use_roaring_if_sparse,
                            double density_threshold )
  : n_cols_bits_(n_cols_bits)
{
    // build the VV
    ardal::detail::WordsVV vv_matrix;
    size_t n_rows = 0, wpr = 0;

    if (is_bitpacked) {
        ardal::utils::log_info("Packed array detected.");
        auto P = make_vv_from_packed_numpy(packed_or_dense_matrix, n_cols_bits);
        vv_matrix = std::move(P.vv);
        n_rows = P.n_rows; wpr = P.wpr;
    } else {
        ardal::utils::log_info("Dense array detected.");
        vv_matrix = pack_dense_to_vv(packed_or_dense_matrix, n_cols_bits);
        n_rows = vv_matrix.size();
        wpr    = vv_matrix.empty() ? 0 : vv_matrix[0].size();
    }

    // metadata
    std::vector<int> row_masses, col_masses;
    int64_t total_mass = 0;
    compute_masses(vv_matrix, n_cols_bits, row_masses, col_masses, total_mass);

    const double density =
        (n_rows && n_cols_bits)
            ? static_cast<double>(total_mass) / (static_cast<double>(n_rows) * n_cols_bits)
            : 0.0;

    // promote to shared ownership (single copy move into shared_ptr)
    packed_matrix_ = std::make_shared<ardal::detail::WordsVV>(std::move(vv_matrix));

    // share masses
    auto rows_sp = std::make_shared<std::vector<int>>(std::move(row_masses));
    auto cols_sp = std::make_shared<std::vector<int>>(std::move(col_masses));

    // init backends without copying vv (pass shared_ptrs)
    if (use_roaring_if_sparse && density < density_threshold) {
        roaring_enabled = true;
        ardal::utils::log_info("Initialising RoaringMatrix backend");
        roaring_backend = std::make_unique<RoaringMatrix>(
            packed_matrix_, rows_sp, cols_sp, n_rows, n_cols_bits);
    } else {
        roaring_enabled = false;
        ardal::utils::log_info("Roaring not initialised.");
    }

    bit_backend = std::make_unique<BitMatrix>(
        vv_, rows_sp, cols_sp, n_rows, n_cols_bits);
}


HybridMatrix::BackendType HybridMatrix::selectBackend( size_t n_cols,
                                                       double density ) const {
    // heuristic backend selection
    if (n_cols > 20000 && density < 0.02) {
        return BackendType::ROARING;
    } else {
        return BackendType::BIT;
    }
}


double HybridMatrix::getDensity( void ) const {
    return density_;
}


py::array_t<uint32_t> HybridMatrix::hamming( bool use_simd,
                                             int threads,
                                             const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        ardal::utils::log_info("Computing Hamming distance matrix with RoaringMatrix backend.");
        return roaring_backend->hamming(threads);
    } else {
        ardal::utils::log_info("Computing Hamming distance matrix with BitMatrix backend.");
        return bit_backend->hamming(false, use_simd, threads);
    }
}


py::array_t<uint32_t> HybridMatrix::hamming_subset( const std::vector<size_t> row_indices,
                                                    const std::vector<size_t> col_indices,
                                                    bool use_simd,
                                                    int threads,
                                                    const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        // note: this assumes local density isnt significantly different from global
        // not a terrible assumption but may not always be the case
        chosen_backend = selectBackend(row_indices.size(), density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->hamming_subset(row_indices, col_indices, threads);
    } else {
        return bit_backend->hamming_subset(row_indices, col_indices, use_simd, threads);
    }
}


py::array_t<double> HybridMatrix::jaccard( bool use_simd,
                                           int threads,
                                           const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        // chosen_backend = selectBackend(n_cols_bits_, density_);
        chosen_backend = BackendType::ROARING;
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->jaccard(threads);
    } else {
        throw std::runtime_error("Jaccard distance is not currently supported in BitMatrix backend.");
    }
}


py::array_t<int> HybridMatrix::innerProduct( bool use_simd,
                                             int threads,
                                             const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->innerProduct(threads);
    } else {
        return bit_backend->innerProduct(false, use_simd, threads);
    }
}


py::array_t<int64_t> HybridMatrix::neighbourhood( size_t row,
                                                  int epsilon,
                                                  bool use_simd,
                                                  int threads,
                                                  const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->neighbourhood(row, epsilon, threads);
    } else {
        return bit_backend->neighbourhood(row, epsilon, use_simd, threads);
    }
}


py::list HybridMatrix::innerProductNeighbourhood( size_t row,
                                                  int ip_epsilon,
                                                  bool use_simd,
                                                  const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->innerProductNeighbourhood(row, ip_epsilon);
    } else {
        return bit_backend->innerProductNeighbourhood(row, ip_epsilon, use_simd);
    }
}


std::vector<size_t> HybridMatrix::uniqueSharedBits( const std::vector<size_t>& row_indices,
                                                    bool use_simd,
                                                    const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    // if (!force_bit_backend) {
    //     return roaring_backend->uniqueSharedBits(row_indices);
    // } else {
    //     return bit_backend->uniqueSharedBits(row_indices, use_simd);
    // }
    return bit_backend->uniqueSharedBits(row_indices, use_simd);
}


py::array_t<double> HybridMatrix::colFrequency( std::vector<size_t>& row_indices ) const {
    return bit_backend->colFrequency(row_indices);
}


py::array_t<double> HybridMatrix::columnEntropy( void ) const {
    return bit_backend->columnEntropy();
}


py::dict HybridMatrix::bitCooccurrence_all( double threshold,
                                            int threads ) const {
    
    if (roaring_enabled) {
        return roaring_backend->bitCooccurrence_all(threshold, threads);
    } else { 
        throw std::runtime_error("Bit co-occurrence is not supported in BitMatrix backend.");
    }
}


py::dict HybridMatrix::bitCooccurrence_subset( const std::vector<size_t>& col_indices,
                                               double threshold,
                                               int threads ) const {
    
    if (roaring_enabled) {
        return roaring_backend->bitCooccurrence_subset(col_indices, threshold, threads);
    } else {
        throw std::runtime_error("Bit co-occurrence is not supported in BitMatrix backend.");
    }
}


py::array_t<double> HybridMatrix::klDivergence( const std::vector<size_t>& input_guids ) const {
    return bit_backend->klDivergence(input_guids);
}


py::array_t<double> HybridMatrix::jsDivergence( const std::vector<size_t>& input_guids ) const {
    return bit_backend->jsDivergence(input_guids);
}


py::array_t<double> HybridMatrix::informationGain( const std::vector<size_t>& input_guids ) const {
    return bit_backend->informationGain(input_guids);
}


py::array_t<size_t> HybridMatrix::getSetBitIndices( size_t row_idx,
                                                    const std::string& backend ) const {
    BackendType chosen_backend = parse_backend(backend);
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->getSetBitIndices(row_idx);
    } else {
        return bit_backend->getSetBitIndices(row_idx);
    }
}


py::array_t<int> HybridMatrix::getRowMasses( void ) const {
    return bit_backend->getRowMasses();
}


py::array_t<int> HybridMatrix::getColumnMasses( void ) const {
    return bit_backend->getColumnMasses();
}


py::array_t<uint8_t> HybridMatrix::getBitMatrix( void ) const {
    return bit_backend->getBitMatrix();
}


py::list HybridMatrix::getRoaringMatrix( void ) const {
    return roaring_backend->getRoaringMatrix();
}


bool HybridMatrix::roaringEnabled( void ) const {
    return roaring_enabled;
}


py::array_t<uint64_t> HybridMatrix::getPackedMatrix( void ) const {
    return bit_backend->getPackedMatrix();
}


py::array_t<size_t> HybridMatrix::getSubsetPackedMatrix( const std::vector<size_t>& row_indices, 
                                                         const std::vector<size_t>& col_indices,
                                                         const int threads ) const {
    return bit_backend->getSubsetMatrix(row_indices, col_indices, threads);
}

} // namespace ardal