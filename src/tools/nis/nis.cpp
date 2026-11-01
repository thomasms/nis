#include <iostream>
#include <map>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "nis/decayheat.hpp"
#include "nis/model.hpp"
#include "nis/solver.hpp"
#include "nis/endf/core.hpp"
#include "nis/endf/decay.hpp"
#include "nis/endf/reaction.hpp"
#include "nis/endf/reaction_cache.hpp"
#include "nis/io/helpers.hpp"

using namespace nis;

namespace {

    struct Options {
        std::string decay_data = "data.dec";
        std::string reaction_data = "";
        std::string input = "run.inp";
        double threshold = 1e10;
    };

    struct DecayData {
        Eigen::SparseMatrix<double> matrix;
        core::DecayEnergyData<double> energy_data;
    };

    // the reaction matrix is linear in the flux, so make it once from the spectrum
    // as it is given, then scale it for the flux of each step
    struct ReactionData {
        Eigen::SparseMatrix<double> matrix;
        double spectrum_total = 0.0;
    };

    void add_options(CLI::App& app, Options& options){
        app.add_option("-d,--decay-data", options.decay_data, "An single ENDF file containing all decay data")
          ->required()
          ->check(CLI::ExistingFile);
        app.add_option("-r,--reaction-data", options.reaction_data, "An single ENDF file containing all reaction data, or a cache of it from nis_cache, needed if the input has a spectrum")
          ->check(CLI::ExistingFile);
        app.add_option("-i,--input", options.input, "Input file for inventory simulation")
          ->required()
          ->check(CLI::ExistingFile);
        app.add_option("-t,--threshold", options.threshold, "Threshold number of atoms for inventory simulation")
          ->check(CLI::PositiveNumber);
    }

    // read the decay tape once, for the matrix and the decay energies
    DecayData read_decay_data(const std::string& filename){
        const auto tape = njoy::ENDFtk::tree::fromFile(filename);
        return DecayData{endf::construct_decay_matrix_from_tape(tape),
                         endf::read_decay_energy_data_from_tape(tape)};
    }

    // no value if the input has no spectrum
    std::optional<ReactionData> read_reaction_data(const core::BasicModel<double>& model,
                                                   const std::string& filename){
        if(model.spectrum == nullptr){
            return std::nullopt;
        }
        if(filename.empty()){
            throw std::invalid_argument("The input has a spectrum, so reaction data (--reaction-data) is required.\n");
        }

        ReactionData data;
        data.matrix = endf::construct_reaction_matrix(endf::read_reaction_data(filename), *model.spectrum.get());
        for(const auto v: model.spectrum->y()){
            data.spectrum_total += v;
        }
        return data;
    }

    // an input with no dt line is one step of dt
    std::vector<core::TimeStep<double>> time_steps(const core::BasicModel<double>& model){
        if(model.steps.empty()){
            return {core::TimeStep<double>{model.dt, std::nullopt}};
        }
        return model.steps;
    }

    // the factor to scale the reaction matrix by for this step,
    // the spectrum is in n/cm^2/s, the step flux in n/m^2/s
    double flux_scale(const core::TimeStep<double>& step, double spectrum_total){
        if(!step.flux.has_value()){
            return 1.0;
        }
        if(step.flux.value() <= 0.0){
            return 0.0;
        }
        if(spectrum_total <= 0.0){
            throw std::invalid_argument("A step has a flux, so the input needs a spectrum with a total above 0.\n");
        }
        return step.flux.value() * 1e-4 / spectrum_total;
    }

    Eigen::SparseMatrix<double> step_matrix(const DecayData& decay,
                                            const std::optional<ReactionData>& reaction,
                                            double scale){
        Eigen::SparseMatrix<double> mat = decay.matrix;
        if(reaction.has_value() && scale > 0.0){
            mat += reaction->matrix * scale;
        }
        return mat;
    }

    core::ModelResult<double> solve_step(const core::ModelResult<double>& previous,
                                         double dt,
                                         const Eigen::SparseMatrix<double>& mat){
        // do computation using CRAM48 as default
        // permit other solvers in the future
        const solver::solver_func<double> solver_function = &core::cram48::solve<double>;

        core::BasicModel<double> step_model;
        step_model.dt = dt;
        step_model.nuclides = previous.nuclides;
        return solver::compute_inventory<double>(step_model, mat, solver_function);
    }

    void print_result(const core::ModelResult<double>& result,
                      const core::DecayEnergyData<double>& energy_data,
                      double threshold){
        // filter out nuclides below a certain threshold
        std::cout << core::filter_nuclides<double>(result, threshold);

        // activity and decay heat use the full inventory, not the filtered one
        const Eigen::VectorX<double> inventory = io::helpers::vector_from_inventory<double>(result.nuclides);
        std::cout << "\nactivity\n    total " << core::compute_activity<double>(energy_data, inventory) << "\n";
        std::cout << core::compute_decay_heat<double>(energy_data, inventory);
    }

    void run(const Options& options){
        const auto model = io::helpers::from_file<double, core::BasicModel>(options.input);
        // dump input to std::out
        // std::cout << model;

        const DecayData decay = read_decay_data(options.decay_data);
        const std::optional<ReactionData> reaction = read_reaction_data(model, options.reaction_data);
        const double spectrum_total = reaction.has_value() ? reaction->spectrum_total : 0.0;
        const auto steps = time_steps(model);

        core::ModelResult<double> result{model.nuclides};
        double time = 0.0;
        for(std::size_t i = 0; i < steps.size(); ++i){
            const auto& step = steps[i];
            const double scale = flux_scale(step, spectrum_total);
            result = solve_step(result, step.dt, step_matrix(decay, reaction, scale));
            time += step.dt;

            if(steps.size() > 1){
                std::cout << "\nstep " << i + 1 << "\n    time " << time << "\n    flux " << step.flux.value_or(scale * spectrum_total * 1e4) << "\n";
            }
            print_result(result, decay.energy_data, options.threshold);
        }
    }

}

int main(int argc, char **argv){

    Options options;
    CLI::App app{"NIS"};
    add_options(app, options);
    CLI11_PARSE(app, argc, argv);

    try{
      run(options);
    }
    catch(const std::exception& ex){
      std::cerr << ex.what();
      return 1;
    }
    catch(...){
      std::cerr << "There was an error - check the logs.\n";
      return 1;
    }

    return 0;
}
