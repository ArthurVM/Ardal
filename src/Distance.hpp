#ifndef DISTANCE_HPP
#define DISTANCE_HPP

#include "src/AlleleMatrix.hpp"
#include "src/DistanceCache.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace _ardal {

class Distance {
 public:

    // Constructor
    Distance(const AlleleMatrix& allele_matrix, DistanceCache& cache)
    : _allele_matrix(allele_matrix), _cache(cache) {}
    ~Distance() { /* Distance class destructor */ }

    // distance methods
    pybind11::array_t<int> hamming( void ) const;
    pybind11::array_t<double> jaccard( void ) const;


 private:
    const AlleleMatrix& _allele_matrix;
    DistanceCache& _cache;

};  // class Distance

}  // namespace _ardal

#endif  // SRC_DISTANCE_HPP_