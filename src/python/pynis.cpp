#include <string>
#include <sstream>
#include <fstream>
#include <functional>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "nis/nis.hpp"

namespace py = pybind11;

using namespace nis;

PYBIND11_MODULE(pynis, m) {

    py::module m_core = m.def_submodule("core");

    m_core.def("decay_mat_from_file", &endf::construct_decay_matrix_from_file);
    m_core.def("compute_inventory", &solver::compute_inventory);

}

