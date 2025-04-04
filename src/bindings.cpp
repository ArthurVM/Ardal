#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "AlleleMatrix.hpp"
#include "Distance.hpp"
#include "DistanceCache.hpp"
#include "Neighbourhood.hpp"

namespace py = pybind11;
namespace _ardal {

PYBIND11_MODULE(_ardal, m) {
    // AlleleMatrix
    py::class_<AlleleMatrix>(m, "AlleleMatrix")
        .def(py::init<py::array_t<uint8_t>&>())
        .def("access", &AlleleMatrix::access, py::arg("coords"))
        .def("accessGUID", &AlleleMatrix::accessGUID, py::arg("guid_index"))
        .def("getMatrix", &AlleleMatrix::getMatrix)
        .def("getNumRows", &AlleleMatrix::getNumRows)
        .def("getNumCols", &AlleleMatrix::getNumCols)
        .def("getMass", &AlleleMatrix::getMass)
        .def("gatherSNPs", &AlleleMatrix::gatherSNPs, py::arg("guid_indices"));

    // DistanceCache
    py::class_<DistanceCache>(m, "DistanceCache")
        .def(py::init<>())
        .def("get", &DistanceCache::get)
        .def("put", &DistanceCache::put)
        .def("clear", &DistanceCache::clear);

    // Distance
    py::class_<Distance>(m, "Distance")
        .def(py::init<const AlleleMatrix&, DistanceCache&>())
        .def("hamming", &Distance::hamming, py::arg("nocache"))
        .def("jaccard", &Distance::jaccard);

    // Neighbourhood
    py::class_<Neighbourhood>(m, "Neighbourhood")
        .def(py::init<const AlleleMatrix&, DistanceCache&>())
        .def("neighbourhood", &Neighbourhood::neighbourhood, py::arg("row_coord"), py::arg("epsilon"), py::arg("nocache"))
        .def("neighbourhoodSIMD", &Neighbourhood::neighbourhoodSIMD, py::arg("row_coord"), py::arg("epsilon"), py::arg("nocache"));
}

} // namespace _ardal
