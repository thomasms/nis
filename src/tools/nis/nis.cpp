#include <iostream>
#include <map>
#include <cmath>

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

int main(int argc, char **argv){

    // cli args
    std::string decay_data = "data.dec";
    std::string reaction_data = "";
    std::string input = "run.inp";
    double threshold = 1e10;

    CLI::App app{"NIS"};
    app.add_option("-d,--decay-data", decay_data, "An single ENDF file containing all decay data")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-r,--reaction-data", reaction_data, "An single ENDF file containing all reaction data, or a cache of it from nis_cache, needed if the input has a spectrum")
      ->check(CLI::ExistingFile);
    app.add_option("-i,--input", input, "Input file for inventory simulation")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-t,--threshold", threshold, "Threshold number of atoms for inventory simulation")
      ->check(CLI::PositiveNumber);
    CLI11_PARSE(app, argc, argv);

    try{
      // read input data
      const auto model = io::helpers::from_file<double, core::BasicModel>(input);
      // dump input to std::out
      // std::cout << model;

      // read the decay tape once, for the matrix and the decay energies
      const auto decay_tape = njoy::ENDFtk::tree::fromFile(decay_data);
      const Eigen::SparseMatrix<double> decay_mat = endf::construct_decay_matrix_from_tape(decay_tape);
      const auto decay_energy_data = endf::read_decay_energy_data_from_tape(decay_tape);

      // the reaction matrix is linear in the flux, so make it once from the spectrum
      // as it is given, then scale it for the flux of each step
      Eigen::SparseMatrix<double> reaction_mat;
      double spectrum_total = 0.0;
      if(model.spectrum != nullptr){
        if(reaction_data.empty()){
          throw std::invalid_argument("The input has a spectrum, so reaction data (--reaction-data) is required.\n");
        }
        reaction_mat = endf::construct_reaction_matrix(endf::read_reaction_data(reaction_data), *model.spectrum.get());
        for(const auto v: model.spectrum->y()){
          spectrum_total += v;
        }
      }

      // an input with no dt line is one step of dt
      std::vector<core::TimeStep<double>> steps = model.steps;
      if(steps.empty()){
        steps.push_back(core::TimeStep<double>{model.dt, std::nullopt});
      }

      const auto print_step = [&](const core::ModelResult<double>& step_result){
        // filter out nuclides below a certain threshold
        std::cout << core::filter_nuclides<double>(step_result, threshold);

        // activity and decay heat use the full inventory, not the filtered one
        const Eigen::VectorX<double> inventory = io::helpers::vector_from_inventory<double>(step_result.nuclides);
        std::cout << "\nactivity\n    total " << core::compute_activity<double>(decay_energy_data, inventory) << "\n";
        std::cout << core::compute_decay_heat<double>(decay_energy_data, inventory);
      };

      // do computation using CRAM48 as default
      // permit other solvers in the future
      solver::solver_func<double> solver_function = &core::cram48::solve<double>;
      core::ModelResult<double> result{model.nuclides};
      double time = 0.0;
      for(std::size_t i = 0; i < steps.size(); ++i){
        const auto& step = steps[i];
        // spectrum in n/cm^2/s, the step flux in n/m^2/s
        double scale = 1.0;
        if(step.flux.has_value()){
          if(step.flux.value() > 0.0 && spectrum_total <= 0.0){
            throw std::invalid_argument("A step has a flux, so the input needs a spectrum with a total above 0.\n");
          }
          scale = (step.flux.value() > 0.0) ? step.flux.value() * 1e-4 / spectrum_total : 0.0;
        }

        Eigen::SparseMatrix<double> mat = decay_mat;
        if(model.spectrum != nullptr && scale > 0.0){
          mat += reaction_mat * scale;
        }

        core::BasicModel<double> step_model;
        step_model.dt = step.dt;
        step_model.nuclides = result.nuclides;
        result = solver::compute_inventory<double>(step_model, mat, solver_function);
        time += step.dt;

        if(steps.size() > 1){
          std::cout << "\nstep " << i + 1 << "\n    time " << time << "\n    flux " << step.flux.value_or(scale * spectrum_total * 1e4) << "\n";
        }
        print_step(result);
      }
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
