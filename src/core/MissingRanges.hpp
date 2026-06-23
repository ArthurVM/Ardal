#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace ardal {

struct MissingRanges {
    std::vector<uint64_t> offsets;
    std::vector<std::pair<uint32_t, uint32_t>> ranges;

    bool empty() const {
        return ranges.empty();
    }
};

} // namespace ardal
