/*
Copyright 2025 Arthur V. Morris
*/
#include <pybind11/pybind11.h>
#include "utils/pack_matrix.hpp"


PYBIND11_MODULE(_ardal_pack, m) {
    m.def("pack_dense_to_words", &pack_dense_to_words, "Pack dense uint8/bool -> <u8 words>");
}