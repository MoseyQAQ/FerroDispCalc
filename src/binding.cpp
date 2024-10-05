#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "get_averaged_structure.hpp"
#include "get_displacement.hpp"
#include "get_polarization.hpp"
namespace py = pybind11;

PYBIND11_MODULE(fdc, m) {
    m.doc() = "Averaging the structure of a trajectory";

    m.def("get_averaged_structure", &get_averaged_structure, "A function to average the structure of a trajectory");
    m.def("get_displacement", &get_displacement, "A function to calculate the displacement of a trajectory");
    m.def("get_polarization", &get_polarization, "A function to calculate the polarization of a trajectory");
}