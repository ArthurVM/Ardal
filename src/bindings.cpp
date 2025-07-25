/*
Copyright 2025 Arthur V. Morris
*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "BitMatrix.hpp"
#include "HybridMatrix.hpp"

namespace py = pybind11;


// pybind functions
PYBIND11_MODULE(_ardal, m) {  // _ardal module and function bindings
    py::class_<_ardal::HybridMatrix>(m, "HybridMatrix")
        .def(py::init<py::array_t<uint8_t>, bool, double>(), 
             py::arg("matrix"), 
             py::arg("use_roaring_if_sparse") = true, 
             py::arg("density_threshold") = 0.02)

        // BitMatrix and RoaringMatrix
        .def("hamming", &_ardal::HybridMatrix::hamming, py::arg("use_simd") = true, py::arg("threads") = 1, py::arg("force_bit_backend") = false)
        .def("jaccard", &_ardal::HybridMatrix::jaccard, py::arg("use_simd") = true, py::arg("threads") = 1, py::arg("force_bit_backend") = false)
        .def("innerProduct", &_ardal::HybridMatrix::innerProduct, py::arg("use_simd") = true, py::arg("threads") = 1, py::arg("force_bit_backend") = false)
        .def("neighbourhood", &_ardal::HybridMatrix::neighbourhood, py::arg("row_idx"), py::arg("epsilon"), py::arg("use_simd") = true, py::arg("threads") = 1, py::arg("force_bit_backend") = false)
        .def("innerProductNeighbourhood", &_ardal::HybridMatrix::innerProductNeighbourhood, py::arg("row_idx"), py::arg("ip_epsilon"), py::arg("use_simd") = true, py::arg("force_bit_backend") = false)
        .def("uniqueSharedBits", &_ardal::HybridMatrix::uniqueSharedBits, py::arg("row_indices"), py::arg("use_simd") = true, py::arg("force_bit_backend") = false)

        // BitMatrix only
        .def("colFrequency", &_ardal::HybridMatrix::colFrequency, py::arg("row_indices"))
        .def("columnEntropy", &_ardal::HybridMatrix::columnEntropy)
        .def("klDivergence", &_ardal::HybridMatrix::klDivergence, py::arg("ingroup_indices"))
        .def("jsDivergence", &_ardal::HybridMatrix::jsDivergence, py::arg("ingroup_indices"))
        .def("informationGain", &_ardal::HybridMatrix::informationGain, py::arg("ingroup_indices"))
        .def("getSetBitIndices", &_ardal::HybridMatrix::getSetBitIndices, py::arg("row_idx"), py::arg("force_bit_backend") = false)
        .def("getRowMasses", &_ardal::HybridMatrix::getRowMasses)
        .def("getColumnMasses", &_ardal::HybridMatrix::getColumnMasses)
        .def("getDensity", &_ardal::HybridMatrix::getDensity)
        .def("bitCooccurrence_all", &_ardal::HybridMatrix::bitCooccurrence_all, py::arg("threshold") = 0.95, py::arg("threads") = 1)
        .def("bitCooccurrence_subset", &_ardal::HybridMatrix::bitCooccurrence_subset, py::arg("col_indices"), py::arg("threshold") = 0.95, py::arg("threads") = 1)
        .def("getBitMatrix", &_ardal::HybridMatrix::getBitMatrix)
        .def("getRoaringMatrix", &_ardal::HybridMatrix::getRoaringMatrix)
        .def("roaringEnabled", &_ardal::HybridMatrix::roaringEnabled);
}