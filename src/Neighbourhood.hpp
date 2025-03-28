#ifndef Neighbourhood_HPP
#define Neighbourhood_HPP

#include "src/AlleleMatrix.hpp"
#include "src/DistanceCache.hpp" // Include DistanceCache
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace _ardal {

class Neighbourhood {
 public:
    // Constructor
    Neighbourhood(const AlleleMatrix& allele_matrix, DistanceCache& cache)
    : _allele_matrix(allele_matrix), _cache(cache) {}
    ~Neighbourhood() { /* Neighbourhood class destructor */ }

    pybind11::array_t<int> neighbourhood(size_t row_coord, int epsilon) const;
    pybind11::list neighbourhoodSIMD(size_t row_coord, int epsilon) const;

 private:
    const AlleleMatrix& _allele_matrix;
    DistanceCache& _cache;

};  // class Neighbourhood

}  // namespace _ardal

#endif  // SRC_Neighbourhood_HPP
