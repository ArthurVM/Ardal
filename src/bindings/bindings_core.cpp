/*
Copyright 2025 Arthur V. Morris
*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../core/BitMatrix.hpp"
#include "../core/HybridMatrix.hpp"
#include "../utils/Logger.hpp"

namespace py = pybind11;


// pybind functions
PYBIND11_MODULE(ardal, m) {  // ardal module and function bindings
    m.def("set_backend_verbosity", &ardal::utils::set_log_level_from_str, "Set backend verbosity");
    py::class_<ardal::HybridMatrix>(m, "HybridMatrix")
        .def(py::init<py::array, bool, size_t, bool, double>(), 
             py::arg("matrix"),
             py::arg("is_bitpacked"),
             py::arg("n_cols_bits"),
             py::arg("use_roaring_if_sparse") = true,
             py::arg("density_threshold") = 0.02)

        // BitMatrix and RoaringMatrix
        .def("hamming", &ardal::HybridMatrix::hamming,
            py::arg("use_simd") = true,
            py::arg("threads") = 1,
            py::arg("backend") = "auto")

        .def("hamming_subset", &ardal::HybridMatrix::hamming_subset,
            py::arg("row_indices"),
            py::arg("col_indices"),
            py::arg("use_simd") = true,
            py::arg("threads") = 1,
            py::arg("backend") = "auto")

        .def("jaccard", &ardal::HybridMatrix::jaccard,
            py::arg("use_simd") = true,
            py::arg("threads") = 1,
            py::arg("backend") = "auto")

        .def("innerProduct", &ardal::HybridMatrix::innerProduct,
            py::arg("use_simd") = true,
            py::arg("threads") = 1,
            py::arg("backend") = "auto")

        .def("neighbourhood", &ardal::HybridMatrix::neighbourhood,
            py::arg("row_idx"),
            py::arg("epsilon"),
            py::arg("use_simd") = true,
            py::arg("threads") = 1,
            py::arg("backend") = "auto")

        .def("innerProductNeighbourhood", &ardal::HybridMatrix::innerProductNeighbourhood,
            py::arg("row_idx"),
            py::arg("ip_epsilon"),
            py::arg("use_simd") = true,
            py::arg("backend") = "auto")

        .def("uniqueSharedBits", &ardal::HybridMatrix::uniqueSharedBits,
            py::arg("row_indices"),
            py::arg("use_simd") = true,
            py::arg("backend") = "auto")


        // BitMatrix only
        .def("colFrequency", &ardal::HybridMatrix::colFrequency,
            py::arg("row_indices"))

        .def("columnEntropy", &ardal::HybridMatrix::columnEntropy)

        .def("klDivergence", &ardal::HybridMatrix::klDivergence,
            py::arg("ingroup_indices"))

        .def("jsDivergence", &ardal::HybridMatrix::jsDivergence,
            py::arg("ingroup_indices"))

        .def("informationGain", &ardal::HybridMatrix::informationGain,
            py::arg("ingroup_indices"))

        .def("getSetBitIndices", &ardal::HybridMatrix::getSetBitIndices,
            py::arg("row_idx"),
            py::arg("backend") = "auto")

        .def("getRowMasses", &ardal::HybridMatrix::getRowMasses)

        .def("getColumnMasses", &ardal::HybridMatrix::getColumnMasses)

        .def("getDensity", &ardal::HybridMatrix::getDensity)

        .def("bitCooccurrence_all", &ardal::HybridMatrix::bitCooccurrence_all,
            py::arg("threshold") = 0.95,
            py::arg("threads") = 1)

        .def("bitCooccurrence_subset", &ardal::HybridMatrix::bitCooccurrence_subset,
            py::arg("col_indices"),
            py::arg("threshold") = 0.95,
            py::arg("threads") = 1)

        .def("getBitMatrix", &ardal::HybridMatrix::getBitMatrix)

        .def("getPackedMatrix", &ardal::HybridMatrix::getPackedMatrix)

        .def("getRoaringMatrix", &ardal::HybridMatrix::getRoaringMatrix)
        
        .def("roaringEnabled", &ardal::HybridMatrix::roaringEnabled)

        .def("getSubsetPackedMatrix", &ardal::HybridMatrix::getSubsetPackedMatrix,
            py::arg("row_indices"),
            py::arg("col_indices"),
            py::arg("threads") = 1);
}