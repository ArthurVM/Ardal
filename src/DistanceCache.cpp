/*
Copyright 2025 Arthur V. Morris
*/

#include "DistanceCache.hpp"
#include <algorithm>

namespace _ardal {



/****************************************************************************************************
 * _ardal::DistanceCache::get
 *
 * Retrieve a Hamming distance from the cache.
 *
 * This function retrieves a previously calculated Hamming distance between two rows from the
 * _hamming_cache. If the distance is not found in the cache, it returns -1.
 *
 * INPUT:
 *   row1 (int) : The index of the first row.
 *   row2 (int) : The index of the second row.
 *
 * OUTPUT:
 *   int : The Hamming distance between the two rows if found in the cache, -1 otherwise.
 ****************************************************************************************************/
int DistanceCache::get(int row1, int row2) const {
    // since (row, col) will be identical to (col, row), use min and max to construct a single key for cache recovery
    auto key = std::make_pair(std::min(row1, row2), std::max(row1, row2));
    
    auto it = _cache.find(key);
    if (it != _cache.end()) {
        return it->second;
    }
    return -1;
}



/****************************************************************************************************
 * _ardal::DistanceCache::put
 *
 * Store a Hamming distance in the cache.
 *
 * This function stores a calculated Hamming distance between two rows in the _hamming_cache.
 * The cache uses a pair of row indices as the key and the Hamming distance as the value.
 *
 * INPUT:
 *   row1 (int) : The index of the first row.
 *   row2 (int) : The index of the second row.
 *   distance (int) : The Hamming distance between the two rows.
 *
 * OUTPUT: None (void)
 ****************************************************************************************************/
void DistanceCache::put(int row1, int row2, int distance) {
    auto key = std::make_pair(std::min(row1, row2), std::max(row1, row2));
    _cache[key] = distance;
}



/****************************************************************************************************
 * _ardal::DistanceCache::clear
 *
 * Clear the Hamming distance cache.
 *
 * This function removes all entries from the _hamming_cache, forcing the recomputation of
 * Hamming distances the next time they are requested.
 *
 * INPUT: None (operates on the private member _hamming_cache)
 *
 * OUTPUT: None (void)
 ****************************************************************************************************/
void DistanceCache::clear( void ) {
    _cache.clear();
}

} // namespace _ardal
