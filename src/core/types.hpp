#pragma once
#include <vector>
#include <cstdint>

namespace ardal::detail {
    using WordsVV = std::vector<std::vector<std::uint64_t>>;

    struct PackedVV {
        std::vector<std::vector<uint64_t>> vv;
        std::size_t n_rows = 0;
        std::size_t wpr    = 0;
    };
}