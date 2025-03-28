#ifndef DISTANCE_CACHE_HPP
#define DISTANCE_CACHE_HPP

#include <map>
#include <utility>

namespace _ardal {

class DistanceCache {
 public:
    DistanceCache() {}
    ~DistanceCache() { /* DistanceCache class destructor */  }

    int get(int row1, int row2) const;
    void put(int row1, int row2, int distance);
    void clear();

 private:
    mutable std::map<std::pair<int, int>, int> _cache;

};  // class DistanceCache

} // namespace _ardal

#endif // DISTANCE_CACHE_HPP
