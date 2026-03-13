/*
Copyright 2025 Arthur V. Morris
*/
#include "HybridMatrix.hpp"
#include <algorithm>
#include <cctype>
#include <utility>
#include "SnvFallbackBackend.hpp"
#include "utils/bitops.hpp"
#include "utils/PythonLogger.hpp"
#include "utils/simd_utils.hpp"


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
                            double density_threshold,
                            py::object missing_mask )
  : n_cols_bits_(n_cols_bits)
{
    const bool packed = is_bitpacked;

    if (packed) {
        // packed -> enforce u64 -> expose as a flat memmap like array with base pointer
        auto arr = pybind11::cast<pybind11::array_t<uint64_t, pybind11::array::c_style>>(packed_or_dense_matrix);
        ardal::detail::build_flat_from_packed( arr,
                                               n_cols_bits_,
                                               flat_,
                                               owner_ );
    } else {
        // dense -> pack into u64 -> expose as a flat memmap like array with base pointer
        ardal::detail::pack_dense_into_storage( packed_or_dense_matrix,
                                                n_cols_bits_,
                                                storage_,
                                                flat_ );
    }

    // metadata
    std::vector<int> row_masses, col_masses;
    int64_t total_mass = 0;
    compute_masses(flat_, row_masses, col_masses, total_mass);

    // std::stringstream ss;
    // ss << "HybridMatrix initialised: flat_.rows=" << flat_.n_rows << " flat_.wpr=" << flat_.wpr << " flat_.n_cols_bits=" << flat_.n_cols_bits;
    // ardal::utils::log_debug(static_cast<string>(ss.str()));

    n_rows_ = flat_.n_rows;
    words_per_row_ = flat_.wpr;
    density_ = static_cast<double>(total_mass) / ((static_cast<double>(flat_.n_rows) * n_cols_bits));

    // share masses
    auto rows_sp = std::make_shared<std::vector<int>>(std::move(row_masses));
    auto cols_sp = std::make_shared<std::vector<int>>(std::move(col_masses));

    if (!missing_mask.is_none()) {
        try {
            py::sequence rows_seq = py::cast<py::sequence>(missing_mask);
            if (py::len(rows_seq) != static_cast<py::ssize_t>(n_rows_)) {
                throw std::runtime_error("missing_mask rows length mismatch");
            }
            mask_columns_.resize(n_rows_);
            bool any = false;
            for (size_t i = 0; i < n_rows_; ++i) {
                py::object row_obj = rows_seq[i];
                if (row_obj.is_none()) {
                    continue;
                }
                py::sequence cols_seq = py::cast<py::sequence>(row_obj);
                auto& out = mask_columns_[i];
                out.reserve(static_cast<size_t>(py::len(cols_seq)));
                for (py::handle item : cols_seq) {
                    int col = py::cast<int>(item);
                    if (col < 0 || static_cast<size_t>(col) >= n_cols_bits_) {
                        continue;
                    }
                    out.push_back(static_cast<uint32_t>(col));
                }
                if (!out.empty()) {
                    std::sort(out.begin(), out.end());
                    out.erase(std::unique(out.begin(), out.end()), out.end());
                    any = true;
                }
            }
            has_missing_mask_ = any;
        } catch (const std::exception& exc) {
            throw std::runtime_error(std::string("HybridMatrix: invalid missing_mask payload: ") + exc.what());
        }
    }

    const std::vector<std::vector<uint32_t>>* mask_ptr = has_missing_mask_ ? &mask_columns_ : nullptr;

    // ----------- INIT BACKENDS -----------
    // handle roaring backend
    if (use_roaring_if_sparse && density_ < density_threshold) {
        roaring_enabled = true;
        ardal::utils::log_info("Initialising RoaringMatrix backend");
        roaring_backend = std::make_unique<RoaringMatrix>( flat_,
                                                           rows_sp,
                                                           cols_sp,
                                                           mask_ptr );
    } else {
        roaring_enabled = false;
        std::stringstream ss;
        ss << "Roaring not initialised: density=" << density_ << "; density_threshold=" << density_threshold;
        ardal::utils::log_debug(static_cast<string>(ss.str()));
    }

    // handle bit backend
    bit_backend = std::make_unique<BitMatrix>( flat_,
                                               rows_sp,
                                               cols_sp,
                                               mask_ptr );
    
    // run SIMD diagnostics
    simd_dispatchers::simd_diag();
}


HybridMatrix::~HybridMatrix() = default;


HybridMatrix::BackendType HybridMatrix::selectBackend( size_t n_cols,
                                                       double density ) const {
    // heuristic backend selection
    if (n_cols > 20000 && density < 0.05 && roaring_enabled) {
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
                                             const std::string& backend,
                                             bool mask_missing ) const {
    BackendType chosen_backend = parse_backend(backend);
    const bool apply_mask = mask_missing && has_missing_mask_;
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
        if (apply_mask) {
            chosen_backend = BackendType::BIT;
        }
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        ardal::utils::log_info("Computing Hamming distance matrix with RoaringMatrix backend.");
        return roaring_backend->hamming(threads, apply_mask);
    } else {
        ardal::utils::log_info("Computing Hamming distance matrix with BitMatrix backend.");
        return bit_backend->hamming(false, use_simd, threads, apply_mask);
    }
}


py::array_t<uint32_t> HybridMatrix::hamming_subset( const std::vector<size_t> row_indices,
                                                    const std::vector<size_t> col_indices,
                                                    bool use_simd,
                                                    int threads,
                                                    const std::string& backend,
                                                    bool mask_missing ) const {
    BackendType chosen_backend = parse_backend(backend);
    const bool apply_mask = mask_missing && has_missing_mask_;
    if (chosen_backend == BackendType::AUTO) {
        // note: this assumes local density isnt significantly different from global
        // not a terrible assumption but may not always be the case
        chosen_backend = selectBackend(row_indices.size(), density_);
        if (apply_mask) {
            chosen_backend = BackendType::BIT;
        }
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->hamming_subset(row_indices, col_indices, threads, apply_mask);
    } else {
        return bit_backend->hamming_subset(row_indices, col_indices, use_simd, threads, apply_mask);
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
                                             const std::string& backend,
                                             bool mask_missing ) const {
    BackendType chosen_backend = parse_backend(backend);
    const bool apply_mask = mask_missing && has_missing_mask_;
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
        if (apply_mask) {
            chosen_backend = BackendType::BIT;
        }
    }
    if (apply_mask && chosen_backend == BackendType::ROARING) {
        throw std::runtime_error("mask_missing is not supported for Roaring innerProduct backend.");
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->innerProduct(threads);
    } else {
        return bit_backend->innerProduct(false, use_simd, threads, mask_missing);
    }
}

py::array_t<double> HybridMatrix::cosineDistance( bool use_simd,
                                                  int threads,
                                                  const std::string& backend,
                                                  bool mask_missing ) const {
    BackendType chosen_backend = parse_backend(backend);
    const bool apply_mask = mask_missing && has_missing_mask_;
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(n_cols_bits_, density_);
        if (apply_mask) {
            chosen_backend = BackendType::BIT;
        }
    }
    if (apply_mask && chosen_backend == BackendType::ROARING) {
        throw std::runtime_error("mask_missing is not supported for Roaring cosine backend.");
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->cosineDistance(threads);
    }
    return bit_backend->cosineDistanceAll(use_simd, threads, mask_missing);
}


py::array_t<double> HybridMatrix::cosineDistance_subset( const std::vector<size_t> row_indices,
                                                         const std::vector<size_t> col_indices,
                                                         bool use_simd,
                                                         int threads,
                                                         const std::string& backend,
                                                         bool mask_missing ) const {
    BackendType chosen_backend = parse_backend(backend);
    const bool apply_mask = mask_missing && has_missing_mask_;
    if (chosen_backend == BackendType::AUTO) {
        chosen_backend = selectBackend(col_indices.size(), density_);
        if (apply_mask) {
            chosen_backend = BackendType::BIT;
        }
    }
    if (apply_mask && chosen_backend == BackendType::ROARING) {
        throw std::runtime_error("mask_missing is not supported for Roaring cosine subset backend.");
    }
    if (chosen_backend == BackendType::ROARING) {
        if (!roaring_enabled)
            throw std::runtime_error("Roaring backend not available.");
        return roaring_backend->cosineDistance_subset(row_indices, col_indices, threads);
    }
    return bit_backend->cosineDistance_subset(row_indices, col_indices, use_simd, threads, mask_missing);
}


py::array_t<int64_t> HybridMatrix::neighbourhood( size_t row,
                                                  uint32_t epsilon,
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


py::array_t<int64_t> HybridMatrix::snvNeighbourhood( size_t row,
                                                     uint32_t epsilon,
                                                     int threads,
                                                     bool mask_missing ) const {
    if (mask_missing) {
        if (!snv_fallback_backend_) {
            throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
        }
        return snv_fallback_backend_->snvNeighbourhood(row, epsilon, threads, true);
    }
    if (roaring_enabled && roaring_backend) {
        return roaring_backend->snvNeighbourhood(row, epsilon, threads, false);
    }
    if (!snv_fallback_backend_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }
    return snv_fallback_backend_->snvNeighbourhood(row, epsilon, threads, false);
}


py::list HybridMatrix::innerProductNeighbourhood( size_t row,
                                                  uint32_t ip_epsilon,
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


py::list HybridMatrix::knn_hamming( size_t row,
                                    uint32_t k,
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
        return roaring_backend->knn_hamming(row, static_cast<int>(k), threads);
    } else {
        return bit_backend->knn_hamming(row, static_cast<int>(k), use_simd, threads);
    }
}


py::list HybridMatrix::knn_inner_product( size_t row,
                                          uint32_t k,
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
        return roaring_backend->knn_inner_product(row, static_cast<int>(k), threads);
    } else {
        return bit_backend->knn_inner_product(row, static_cast<int>(k), use_simd, threads);
    }
}


py::list HybridMatrix::knn_jaccard( size_t row,
                                    uint32_t k,
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
        return roaring_backend->knn_jaccard(row, static_cast<int>(k), threads);
    } else {
        return bit_backend->knn_jaccard(row, static_cast<int>(k), use_simd, threads);
    }
}


py::list HybridMatrix::knnSnv( size_t row,
                               uint32_t k,
                               int threads,
                               bool mask_missing ) const {
    if (mask_missing) {
        if (!snv_fallback_backend_) {
            throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
        }
        return snv_fallback_backend_->knnSnv(row, k, threads, true);
    }
    if (roaring_enabled && roaring_backend) {
        return roaring_backend->knnSnv(row, static_cast<int>(k), threads, false);
    }
    if (!snv_fallback_backend_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }
    return snv_fallback_backend_->knnSnv(row, k, threads, false);
}


py::list HybridMatrix::knn_cosine( size_t row,
                                   uint32_t k,
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
        return roaring_backend->knn_cosine(row, static_cast<int>(k), threads);
    } else {
        return bit_backend->knn_cosine(row, static_cast<int>(k), use_simd, threads);
    }
}


py::list HybridMatrix::knn( size_t row,
                            uint32_t k,
                            const std::string& metric,
                            bool use_simd,
                            int threads,
                            const std::string& backend ) const {
    std::string metric_lc = metric;
    std::transform(metric_lc.begin(), metric_lc.end(), metric_lc.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });

    if (metric_lc == "hamming") {
        return knn_hamming(row, k, use_simd, threads, backend);
    } else if (metric_lc == "inner_product" || metric_lc == "innerproduct") {
        return knn_inner_product(row, k, use_simd, threads, backend);
    } else if (metric_lc == "jaccard") {
        return knn_jaccard(row, k, use_simd, threads, backend);
    } else if (metric_lc == "cosine") {
        return knn_cosine(row, k, use_simd, threads, backend);
    } else if (metric_lc == "snv") {
        return knnSnv(row, k, threads);
    }
    throw std::invalid_argument("Unsupported metric for knn: " + metric);
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


void HybridMatrix::prepareSnvView( py::array_t<uint32_t> allele_to_locus,
                                   py::array_t<uint8_t> allele_to_base ) {
    py::array_t<uint32_t> locus_for_roaring = allele_to_locus;
    py::array_t<uint8_t> base_for_roaring = allele_to_base;
    if (roaring_enabled) {
        roaring_backend->prepareSnvView(std::move(locus_for_roaring), std::move(base_for_roaring));
    }
    if (!snv_fallback_backend_) {
        snv_fallback_backend_ = std::make_unique<SnvFallbackBackend>(flat_, n_cols_bits_, has_missing_mask_ ? &mask_columns_ : nullptr);
    }
    snv_fallback_backend_->prepare(std::move(allele_to_locus), std::move(allele_to_base));
}


py::array_t<uint32_t> HybridMatrix::snvHamming( int threads,
                                                bool mask_missing ) const {
    if (mask_missing) {
        if (!snv_fallback_backend_) {
            throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
        }
        return snv_fallback_backend_->snvHamming(threads, true);
    }
    if (roaring_enabled) {
        return roaring_backend->snvHamming(threads, false);
    }
    if (!snv_fallback_backend_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }
    return snv_fallback_backend_->snvHamming(threads, false);
}


py::array_t<uint32_t> HybridMatrix::snvHamming_subset( const std::vector<size_t>& row_indices,
                                                      const std::vector<size_t>& col_indices,
                                                      int threads,
                                                      bool mask_missing ) const {
    if (mask_missing) {
        if (!snv_fallback_backend_) {
            throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
        }
        return snv_fallback_backend_->snvHamming_subset(row_indices, col_indices, threads, true);
    }
    if (roaring_enabled) {
        return roaring_backend->snvHamming_subset(row_indices, col_indices, threads, false);
    }
    if (!snv_fallback_backend_) {
        throw std::runtime_error("SNV lookup not initialised. Call prepareSnvView first.");
    }
    return snv_fallback_backend_->snvHamming_subset(row_indices, col_indices, threads, false);
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


py::array_t<uint64_t> HybridMatrix::getSubsetPackedMatrix_rows( const std::vector<size_t>& row_indices, 
                                                                const int threads,
                                                                bool silence_log ) const {
    return bit_backend->getSubsetMatrix_rows(row_indices, threads, silence_log);
}

} // namespace ardal
