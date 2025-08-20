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
#include "nis/io/helpers.hpp"

using namespace nis;

int main(int argc, char **argv){

    // cli args
    std::string data = "data.dec";
    std::string input = "run.inp";
    double threshold = 1e10;

    CLI::App app{"NIS"};
    app.add_option("-d,--data", data, "An single ENDF file containing all decay data")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-i,--input", input, "Input file for inventory simulation")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-t,--threshold", threshold, "Threshold number of atoms for inventory simulation")
      ->check(CLI::PositiveNumber);
    CLI11_PARSE(app, argc, argv);

    try{
      // construct a decay matrix - straight from ENDF file(s)
      const Eigen::MatrixX<double> dense_decay_mat = endf::construct_decay_matrix_from_file(data);

      // read input data
      const auto model = io::helpers::from_file<double, core::BasicModel>(input);

      // do computation using CRAM48 as default
      // permit other solvers in the future
      solver::solver_func<double> solver_function = &core::cram48::solve<double>;
      const auto result = solver::compute_inventory<double>(model, dense_decay_mat, solver_function);

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
