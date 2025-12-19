#ifndef PYTHON_LOGGER_HPP
#define PYTHON_LOGGER_HPP

#include <pybind11/pybind11.h>
#include <string>
#include "Logger.hpp"

using namespace std;
namespace py = pybind11;


namespace ardal {
namespace utils {

// Function to log messages using Python's logging module
inline void log_info(const std::string& message) {
    if (current_log_level > INFO) return;
    try {
        py::object logger = py::module_::import("logging").attr("getLogger")("ardal");
        logger.attr("info")("[BACKEND] " + message);
    } catch (py::error_already_set& e) {
        // Handle cases where the logging module can't be imported or used
        py::print("Error in C++ logging:", e.what());
    }
}

inline void log_warning(const std::string& message) {
    if (current_log_level > WARNING) return;
    try {
        py::object logger = py::module_::import("logging").attr("getLogger")("ardal");
        logger.attr("warning")("[BACKEND] " + message);
    } catch (py::error_already_set& e) {
        py::print("Error in C++ logging:", e.what());
    }
}

inline void log_error(const std::string& message) {
    if (current_log_level > ERROR) return;
    try {
        py::object logger = py::module_::import("logging").attr("getLogger")("ardal");
        logger.attr("error")("[BACKEND] " + message);
    } catch (py::error_already_set& e) {
        py::print("Error in C++ logging:", e.what());
    }
}

inline void log_debug(const std::string& message) {
    if (current_log_level > DEBUG) return;
    try {
        py::object logger = py::module_::import("logging").attr("getLogger")("ardal");
        logger.attr("debug")("[BACKEND] " + message);
    } catch (py::error_already_set& e) {
        py::print("Error in C++ logging:", e.what());
    }
}

} // namespace utils
} // namespace ardal

#endif // PYTHON_LOGGER_HPP
