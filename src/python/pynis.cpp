#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "nis/nis.hpp"

// The nuclear data stays in C++, so it is not copied, and a change in python changes the data:
// i.e. decay_data[zai].half_life = 10.0 or reaction_data[i].values = new_values.
// a list of CrossSection in python
PYBIND11_MAKE_OPAQUE(nis::data::ReactionData);
// a dict of FissionYields by parent ZAI in python
PYBIND11_MAKE_OPAQUE(nis::data::FissionYieldData);
// a dict of DecayNuclide by ZAI in python
PYBIND11_MAKE_OPAQUE(nis::data::DecayData);
// a list of DecayMode in python
PYBIND11_MAKE_OPAQUE(std::vector<nis::data::DecayMode>);

namespace py = pybind11;

using namespace nis;

namespace {

    std::string _repr_mode(const data::DecayMode& mode){
        std::ostringstream os;
        os << "DecayMode(rtyp=" << mode.rtyp << ", branching_ratio=" << mode.branching_ratio
           << ", final_state=" << mode.final_state << ")";
        return os.str();
    }

    // The nuclear data types, with no dependency on ENDF. Each can be read from a file, or made in python,
    // then changed, and then given to decay_mat, decay_energy_data or reaction_mat.
    void bind_nuclear_data(py::module& m_core){

        // ---- decay ----
        py::class_<data::DecayMode>(m_core, "DecayMode",
            "One decay mode of a nuclide. rtyp is the decay type with the ENDF numbers: 0 gamma, 1 beta-, "
            "2 beta+ (and EC), 3 IT, 4 alpha, 5 neutron, 6 spontaneous fission, 7 proton, 10 unknown (to the sink). "
            "A decay that is more than one of these, one after the other, has more digits: 1.5 is beta- then neutron.")
            .def(py::init<>())
            .def(py::init([](double rtyp, double branching_ratio, int final_state){
                    return data::DecayMode{rtyp, branching_ratio, final_state};
                 }), py::arg("rtyp"), py::arg("branching_ratio"), py::arg("final_state") = 0)
            .def_readwrite("rtyp", &data::DecayMode::rtyp)
            .def_readwrite("branching_ratio", &data::DecayMode::branching_ratio)
            .def_readwrite("final_state", &data::DecayMode::final_state, "the isomeric state of the daughter")
            .def("__repr__", &_repr_mode);

        // a change to an item of the list changes the nuclide, i.e. nuclide.modes[0].branching_ratio = 0.5
        py::bind_vector<std::vector<data::DecayMode>>(m_core, "DecayModeList");
        // so a python list of DecayMode can be given, i.e. nuclide.modes = [DecayMode(1, 1.0)]
        py::implicitly_convertible<py::list, std::vector<data::DecayMode>>();

        py::class_<data::DecayNuclide>(m_core, "DecayNuclide",
            "The decay data of one nuclide. half_life is in s, infinity is stable. "
            "The average energies per decay are in eV, they are only used for decay heat.")
            .def(py::init<>())
            .def(py::init([](double half_life, const std::vector<data::DecayMode>& modes,
                             double light_particle, double electromagnetic, double heavy_particle){
                    data::DecayNuclide nuclide;
                    nuclide.half_life = half_life;
                    nuclide.modes = modes;
                    nuclide.light_particle = light_particle;
                    nuclide.electromagnetic = electromagnetic;
                    nuclide.heavy_particle = heavy_particle;
                    return nuclide;
                 }), py::arg("half_life"), py::arg("modes") = std::vector<data::DecayMode>(),
                 py::arg("light_particle") = 0.0, py::arg("electromagnetic") = 0.0, py::arg("heavy_particle") = 0.0)
            .def_readwrite("half_life", &data::DecayNuclide::half_life)
            .def_readwrite("modes", &data::DecayNuclide::modes)
            .def_readwrite("light_particle", &data::DecayNuclide::light_particle)
            .def_readwrite("electromagnetic", &data::DecayNuclide::electromagnetic)
            .def_readwrite("heavy_particle", &data::DecayNuclide::heavy_particle)
            .def_property("decay_constant", &data::DecayNuclide::decay_constant, &data::DecayNuclide::set_decay_constant,
                          "1/s, 0 is stable. Setting it sets the half-life.")
            .def_property_readonly("is_stable", &data::DecayNuclide::is_stable)
            .def_property_readonly("total_branching_ratio", &data::DecayNuclide::total_branching_ratio,
                                   "the sum of the branching ratios, 1 for good data")
            .def("__repr__", [](const data::DecayNuclide& n){
                std::ostringstream os;
                os << "DecayNuclide(half_life=" << n.half_life << ", modes=[";
                for(std::size_t i = 0; i < n.modes.size(); ++i){
                    os << (i > 0 ? ", " : "") << _repr_mode(n.modes[i]);
                }
                os << "])";
                return os.str();
            });

        // a dict of DecayNuclide by ZAI. decay_data[zai] gives the nuclide in the data, not a copy,
        // so decay_data[zai].half_life = 1.0 changes the data. Setting an item copies the nuclide in.
        py::bind_map<data::DecayData>(m_core, "DecayData");

        m_core.def("decay_data_from_file", &endf::read_decay_data_from_file, py::arg("filename"),
                   "all the nuclides of an ENDF decay file (MF=8 MT=457), as DecayData");
        // matrices are sparse, so they are scipy.sparse.csc_matrix in python
        // with sf_yields (from fission_yields_from_file of a spontaneous fission yield file),
        // spontaneous fission makes the fission products, else it goes to the sink
        m_core.def("decay_mat", &data::construct_decay_matrix,
                   py::arg("decay_data"), py::arg("sf_yields") = nullptr);
        m_core.def("decay_energy_data", &data::construct_decay_energy_data, py::arg("decay_data"));
        m_core.def("decay_mat_from_file", &endf::construct_decay_matrix_from_file,
                   py::arg("filename"), py::arg("sf_yields") = nullptr);
        m_core.def("expand_rtypes", &data::expand_rtypes, py::arg("rtyp"));
        // (primary, [secondaries]), a primary of -1 is the sink
        m_core.def("decay_products", &data::compute_products,
                   py::arg("parent_zai"), py::arg("rtyp"), py::arg("final_state") = 0);

        // ---- reactions ----
        m_core.attr("ALL_STATES") = py::int_(data::ALL_STATES);
        m_core.attr("FISSION_MT") = py::int_(data::FISSION_MT);
        m_core.attr("INELASTIC_MT") = py::int_(data::INELASTIC_MT);
        m_core.attr("ANYTHING_MT") = py::int_(data::ANYTHING_MT);

        py::class_<data::CrossSection>(m_core, "CrossSection",
            "One group-wise cross section, for one parent and one reaction (ENDF MT number). "
            "energies are the group bounds in eV, values are in barns, values[g] is the group from energies[g] to "
            "energies[g+1], so the last value is not used. final_state is ALL_STATES for the total of the reaction, "
            "else the isomeric state of the residual. "
            "Note: energies and values are copied to python, so set them as a whole: xs.values = new_values.")
            .def(py::init([](int parent_zai, int mt_number, const std::vector<double>& energies,
                             const std::vector<double>& values, int final_state, int product_za){
                    return data::CrossSection{parent_zai, mt_number, energies, values, final_state, product_za};
                 }), py::arg("parent_zai"), py::arg("mt_number"), py::arg("energies"), py::arg("values"),
                 py::arg("final_state") = data::ALL_STATES, py::arg("product_za") = 0)
            .def_readwrite("parent_zai", &data::CrossSection::parent_zai)
            .def_readwrite("mt_number", &data::CrossSection::mt_number)
            .def_readwrite("energies", &data::CrossSection::energies)
            .def_readwrite("values", &data::CrossSection::values)
            .def_readwrite("final_state", &data::CrossSection::final_state)
            .def_readwrite("product_za", &data::CrossSection::product_za)
            .def("__repr__", [](const data::CrossSection& xs){
                std::ostringstream os;
                os << "CrossSection(parent_zai=" << xs.parent_zai << ", mt_number=" << xs.mt_number
                   << ", groups=" << (xs.energies.empty() ? 0 : xs.energies.size() - 1)
                   << ", final_state=" << xs.final_state << ", product_za=" << xs.product_za << ")";
                return os.str();
            });

        // A list of CrossSection that stays in C++ memory, as it is big. reaction_data[i] gives the cross section
        // in the data, not a copy. Do not keep these after you add to the list, as they can then become invalid.
        py::bind_vector<data::ReactionData>(m_core, "ReactionData")
            .def("find", [](py::object self, int parent_zai, std::optional<int> mt_number){
                    auto& reaction_data = self.cast<data::ReactionData&>();
                    py::list found;
                    for(auto& xs: reaction_data){
                        if((xs.parent_zai == parent_zai) && (!mt_number || (xs.mt_number == *mt_number))){
                            found.append(py::cast(&xs, py::return_value_policy::reference_internal, self));
                        }
                    }
                    return found;
                 }, py::arg("parent_zai"), py::arg("mt_number") = py::none(),
                 "the cross sections of the parent (and the reaction, if given), not copies")
            .def("remove_reactions", [](data::ReactionData& reaction_data, int parent_zai, std::optional<int> mt_number){
                    const std::size_t size_before = reaction_data.size();
                    reaction_data.erase(std::remove_if(reaction_data.begin(), reaction_data.end(),
                        [&](const data::CrossSection& xs){
                            return (xs.parent_zai == parent_zai) && (!mt_number || (xs.mt_number == *mt_number));
                        }), reaction_data.end());
                    return size_before - reaction_data.size();
                 }, py::arg("parent_zai"), py::arg("mt_number") = py::none(),
                 "removes the cross sections of the parent (and the reaction, if given), returns how many");

        // a cache or an ENDF tape
        m_core.def("reaction_data_from_file", &endf::read_reaction_data, py::arg("filename"));
        m_core.def("write_reaction_cache", &endf::write_reaction_cache, py::arg("reaction_data"), py::arg("filename"));
        m_core.def("is_reaction_cache", &endf::is_reaction_cache, py::arg("filename"));
        // with fission_yields (from fission_yields_from_file), fission makes the fission products,
        // else it goes to the sink
        m_core.def("reaction_mat", &data::construct_reaction_matrix,
                   py::arg("reaction_data"), py::arg("spectrum"), py::arg("fission_yields") = nullptr);
        m_core.def("reaction_mat_from_file", [](const std::string& filename, const core::Spectrum<double>& spectrum,
                                                const data::FissionYieldData* fission_yields){
            return data::construct_reaction_matrix(endf::read_reaction_data(filename), spectrum, fission_yields);
        }, py::arg("filename"), py::arg("spectrum"), py::arg("fission_yields") = nullptr);
        m_core.def("is_reaction_mt", &data::is_reaction_mt, py::arg("mt_number"),
                   "true if the MT is used for the matrix as a total, see REACTION_MAP. MT=4 and MT=5 are used per state");
        // (primary, [secondaries]), a primary of -1 is the sink
        m_core.def("reaction_products", &data::compute_reaction_products,
                   py::arg("parent_zai"), py::arg("mt_number"), py::arg("final_state") = 0);
        m_core.def("collapse_cross_section", &data::collapse_cross_section);

        // ---- fission yields ----
        py::class_<data::FissionYields>(m_core, "FissionYields",
            "The independent fission yields of one parent. energies are the incident neutron energies in eV, in order. "
            "yields has one dict (product ZAI -> yield) for each energy. "
            "Note: energies and yields are copied to python, so set them as a whole.")
            .def(py::init<>())
            .def(py::init([](int parent_zai, const std::vector<double>& energies,
                             const std::vector<std::map<int, double>>& yields){
                    if(energies.size() != yields.size()){
                        throw std::invalid_argument("There must be one set of yields for each energy");
                    }
                    return data::FissionYields{parent_zai, energies, yields};
                 }), py::arg("parent_zai"), py::arg("energies"), py::arg("yields"))
            .def_readwrite("parent_zai", &data::FissionYields::parent_zai)
            .def_readwrite("energies", &data::FissionYields::energies)
            .def_readwrite("yields", &data::FissionYields::yields);

        // a dict of FissionYields by parent ZAI, fission_yields[zai] gives the yields in the data, not a copy
        py::bind_map<data::FissionYieldData>(m_core, "FissionYieldData");

        m_core.def("fission_yields_from_file", &endf::read_fission_yields_from_file, py::arg("filename"));
        m_core.def("fission_yield_energy_weights", &data::energy_weights);
        m_core.def("fission_yields_at_energy", &data::yields_at_energy);
        m_core.def("weighted_fission_yields", &data::weighted_yields);
    }
}

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

    using TimeStepPyType = core::TimeStep<double>;
    py::class_<TimeStepPyType>(m_core, "TimeStep")
        .def(py::init<>())
        .def(py::init<double, std::optional<double>>(), py::arg("dt"), py::arg("flux") = py::none())
        .def_readwrite("dt", &TimeStepPyType::dt)
        .def_readwrite("flux", &TimeStepPyType::flux);

    using BasicModelPyType = core::BasicModel<double>;
    py::class_<BasicModelPyType>(m_core, "BasicModel")
        .def(py::init<>())
        .def_readwrite("nuclides", &BasicModelPyType::nuclides)
        .def_readwrite("dt", &BasicModelPyType::dt)
        .def_readwrite("spectrum", &BasicModelPyType::spectrum)
        .def_readwrite("steps", &BasicModelPyType::steps);

    using ModelResultPyType = core::ModelResult<double>;
    py::class_<ModelResultPyType>(m_core, "ModelResult")
        .def(py::init<>())
        .def_readwrite("nuclides", &ModelResultPyType::nuclides);

    m_core.def("filter_nuclides", &core::filter_nuclides<double>);


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

    m_core.def("decay_energy_data_from_file", &endf::read_decay_energy_data_from_file, py::arg("filename"));
    m_core.def("compute_activities", &core::compute_activities<double>);
    m_core.def("compute_activity", &core::compute_activity<double>);
    m_core.def("compute_decay_heat", &core::compute_decay_heat<double>);
    m_core.def("vector_from_inventory", &io::helpers::vector_from_inventory<double>);

    bind_nuclear_data(m_core);

    // solvers
    m_core.def("solver_cram48", &core::cram48::solve<double>);
    m_core.def("reachable_indices", &solver::reachable_indices<double>);
    m_core.def("solve_reachable", &solver::solve_reachable<double>);
    m_core.def("compute_inventory", py::overload_cast<const BasicModelPyType&, const Eigen::SparseMatrix<double>&, const solver::solver_func<double>&>(&solver::compute_inventory<double>));

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
