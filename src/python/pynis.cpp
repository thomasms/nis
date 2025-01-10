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

    m.doc() = R"pbdoc(
        NIS library using pybind11

        ToDo: write

        A list:
            - Entry 1
            - Entry 2

    )pbdoc"; // optional module docstring

    // version
    m.attr("__version__") = "0.0.1";

    py::module m_core = m.def_submodule("core");

    using BasicModelPyType = core::BasicModel<double>;
    py::class_<BasicModelPyType>(m_core, "BasicModel")
        .def(py::init<>())
        .def_readwrite("nuclides", &BasicModelPyType::nuclides)
        .def_readwrite("dt", &BasicModelPyType::dt);

    using ModelResultPyType = core::ModelResult<double>;
    py::class_<ModelResultPyType>(m_core, "ModelResult")
        .def(py::init<>())
        .def_readwrite("nuclides", &ModelResultPyType::nuclides);

    m_core.def("decay_mat_from_file", &endf::construct_decay_matrix_from_file);
    m_core.def("compute_inventory", py::overload_cast<const BasicModelPyType&, const Eigen::MatrixX<double>&>(&solver::compute_inventory<double>));

    py::module m_data = m.def_submodule("data");
    // would be nice to make these all read only attributes, if possible
    m_data.def("ALL_ZAIS",  []() { return data::ALL_ZAIS; });
    m_data.def("Z_ELEMENT_MAP",  []() { return data::Z_ELEMENT_MAP; });
    m_data.attr("ZAI_SIZE") = py::int_(data::ZAI_SIZE);
    m_data.attr("SINK_INDEX") = py::int_(data::SINK_INDEX);

    py::module m_zai = m.def_submodule("zai");
    m_zai.def("META_STABLE_MAP",  []() { return zai::META_STABLE_MAP; });
    m_zai.def("get_z", &zai::get_z);
    m_zai.def("get_a", &zai::get_a);
    m_zai.def("get_i", &zai::get_i);
    m_zai.def("make", py::overload_cast<int, int, int>(&zai::make));
    m_zai.def("decompose", &zai::decompose);
    m_zai.def("get_from_name", &zai::get_from_name);
    m_zai.def("get_name_from", &zai::get_name_from);
}
