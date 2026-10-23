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

// keep the cross sections in C++, they are big
PYBIND11_MAKE_OPAQUE(nis::endf::ReactionData);

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

    using SpectrumPyType = core::Spectrum<double>;
    py::class_<SpectrumPyType, std::shared_ptr<SpectrumPyType>>(m_core, "Spectrum")
        .def(py::init<const std::vector<double>&, const std::vector<double>&>(), py::arg("x"), py::arg("y"))
        .def_property_readonly("x", &SpectrumPyType::xcopy)
        .def_property_readonly("y", &SpectrumPyType::ycopy);

    m_core.def("convert_by_energy", &core::convert_by_energy<double>);
    m_core.def("convert_by_lethargy", &core::convert_by_lethargy<double>);

    using BasicModelPyType = core::BasicModel<double>;
    py::class_<BasicModelPyType>(m_core, "BasicModel")
        .def(py::init<>())
        .def_readwrite("nuclides", &BasicModelPyType::nuclides)
        .def_readwrite("dt", &BasicModelPyType::dt)
        .def_readwrite("spectrum", &BasicModelPyType::spectrum);

    using ModelResultPyType = core::ModelResult<double>;
    py::class_<ModelResultPyType>(m_core, "ModelResult")
        .def(py::init<>())
        .def_readwrite("nuclides", &ModelResultPyType::nuclides);

    m_core.def("filter_nuclides", &core::filter_nuclides<double>);

    m_core.def("decay_mat_from_file", &endf::construct_decay_matrix_from_file);

    using DecayEnergyDataPyType = core::DecayEnergyData<double>;
    py::class_<DecayEnergyDataPyType>(m_core, "DecayEnergyData")
        .def(py::init<>())
        .def_readwrite("decay_constants", &DecayEnergyDataPyType::decay_constants)
        .def_readwrite("light_particle", &DecayEnergyDataPyType::light_particle)
        .def_readwrite("electromagnetic", &DecayEnergyDataPyType::electromagnetic)
        .def_readwrite("heavy_particle", &DecayEnergyDataPyType::heavy_particle);

    using DecayHeatPyType = core::DecayHeat<double>;
    py::class_<DecayHeatPyType>(m_core, "DecayHeat")
        .def(py::init<>())
        .def_readwrite("light_particle", &DecayHeatPyType::light_particle)
        .def_readwrite("electromagnetic", &DecayHeatPyType::electromagnetic)
        .def_readwrite("heavy_particle", &DecayHeatPyType::heavy_particle)
        .def_property_readonly("total", &DecayHeatPyType::total);

    m_core.def("decay_energy_data_from_file", &endf::read_decay_energy_data_from_file);
    m_core.def("compute_activities", &core::compute_activities<double>);
    m_core.def("compute_activity", &core::compute_activity<double>);
    m_core.def("compute_decay_heat", &core::compute_decay_heat<double>);
    m_core.def("vector_from_inventory", &io::helpers::vector_from_inventory<double>);

    py::class_<endf::CrossSection>(m_core, "CrossSection")
        .def_readonly("parent_zai", &endf::CrossSection::parent_zai)
        .def_readonly("mt_number", &endf::CrossSection::mt_number)
        .def_readonly("energies", &endf::CrossSection::energies)
        .def_readonly("values", &endf::CrossSection::values);

    // ReactionData is opaque so it stays in C++ memory and is not copied to a python list
    py::class_<endf::ReactionData>(m_core, "ReactionData")
        .def("__len__", [](const endf::ReactionData& d) { return d.size(); })
        .def("__getitem__", [](const endf::ReactionData& d, size_t i) {
            if(i >= d.size()){
                throw py::index_error();
            }
            return d[i];
        });

    m_core.def("reaction_data_from_file", &endf::read_reaction_data_from_file);
    m_core.def("reaction_mat", &endf::construct_reaction_matrix);
    m_core.def("reaction_mat_from_file", &endf::construct_reaction_matrix_from_file);
    m_core.def("collapse_cross_section", &endf::collapse_cross_section);

    // solvers
    m_core.def("solver_cram48", &core::cram48::solve<double>);
    m_core.def("compute_inventory", py::overload_cast<const BasicModelPyType&, const Eigen::MatrixX<double>&, solver::solver_func<double>&>(&solver::compute_inventory<double>));

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
