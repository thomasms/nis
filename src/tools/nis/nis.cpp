#include <iostream>
#include <map>
#include <cmath>

#include <CLI/CLI.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

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

      const Eigen::MatrixX<double> dense_decay_mat = endf::construct_decay_matrix_from_file(decay_data);

      Eigen::MatrixX<double> dense_reaction_mat = Eigen::MatrixXd::Zero(dense_decay_mat.rows(), dense_decay_mat.cols());
      if(model.spectrum != nullptr){
        if(reaction_data.empty()){
          throw std::invalid_argument("The input has a spectrum, so reaction data (--reaction-data) is required.\n");
        }
        dense_reaction_mat = endf::construct_reaction_matrix_from_file(reaction_data, *model.spectrum.get());
      }

      Eigen::MatrixX<double> dense_mat = dense_decay_mat + dense_reaction_mat;

      // do computation using CRAM48 as default
      // permit other solvers in the future
      solver::solver_func<double> solver_function = &core::cram48::solve<double>;
      const auto result = solver::compute_inventory<double>(model, dense_mat, solver_function);

      // filter out nuclides below a certain threshold
      const auto filtered_result = core::filter_nuclides<double>(result, threshold);

      // dump output to std::out
      std::cout << filtered_result;
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
