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
    app.add_option("-r,--reaction-data", reaction_data, "An single ENDF file containing all reaction data, needed if the input has a spectrum")
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

      Eigen::SparseMatrix<double> mat = decay_mat;
      if(model.spectrum != nullptr){
        if(reaction_data.empty()){
          throw std::invalid_argument("The input has a spectrum, so reaction data (--reaction-data) is required.\n");
        }
        mat += endf::construct_reaction_matrix_from_file(reaction_data, *model.spectrum.get());
      }

      // do computation using CRAM48 as default
      // permit other solvers in the future
      solver::solver_func<double> solver_function = &core::cram48::solve<double>;
      const auto result = solver::compute_inventory<double>(model, mat, solver_function);

      // filter out nuclides below a certain threshold
      const auto filtered_result = core::filter_nuclides<double>(result, threshold);

      // dump output to std::out
      std::cout << filtered_result;

      // activity and decay heat use the full inventory, not the filtered one
      const Eigen::VectorX<double> final_inventory = io::helpers::vector_from_inventory<double>(result.nuclides);
      std::cout << "\nactivity\n    total " << core::compute_activity<double>(decay_energy_data, final_inventory) << "\n";
      std::cout << core::compute_decay_heat<double>(decay_energy_data, final_inventory);
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
