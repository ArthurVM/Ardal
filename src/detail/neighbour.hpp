// ardal/detail/neighbour.hpp
#pragma once
#include <cstdint>

namespace ardal::knn {

struct Neighbor {
    using IdT = uint32_t;
    using DistT = uint32_t;

    IdT   id;
    DistT d;
};

// max heap by distance 
// keeps worst/top-k threshold on top
struct ByMaxDistance {
    bool operator()(const Neighbor& a, const Neighbor& b) const noexcept {
        // "true" means a has LOWER priority than b
        return a.d < b.d; // => largest d is top
    }
};

// min heap by distance
// probably not that useful but might be useful elsewhere
// NOTE: not currently in use in Ardal
struct ByMinDistance {
    bool operator()(const Neighbor& a, const Neighbor& b) const noexcept {
        return a.d > b.d; // => smallest d is top
    }
};

// final output ordering: ascending distance, tie-break by id
struct AscDistanceId {
    bool operator()(const Neighbor& x, const Neighbor& y) const noexcept {
        return (x.d < y.d) || (x.d == y.d && x.id < y.id);
    }
};

// helpers
inline bool equal_id_dist(const Neighbor& a, const Neighbor& b) noexcept {
    return a.id == b.id && a.d == b.d;
}

} // namespace ardal::knn
