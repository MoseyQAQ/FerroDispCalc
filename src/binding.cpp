#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "get_averaged_structure.hpp"

namespace py = pybind11;

PYBIND11_MODULE(fdc, m) {
    m.doc() = "Averaging the structure of a trajectory";

    m.def("get_averaged_structure", &get_averaged_structure, "A function to average the structure of a trajectory");
}