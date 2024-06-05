#include <iostream>
#include <map>
#include <cmath>

#include <CLI/CLI.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "nis/io/helpers.hpp"
#include "nis/io/input.hpp"
#include "nis/io/output.hpp"
#include "nis/endf/core.hpp"
#include "nis/solver.hpp"

using namespace nis;

int main(int argc, char **argv){

    // cli args
    std::string data = "data.dec";
    std::string input = "run.inp";

    CLI::App app{"NIS"};
    app.add_option("-d,--data", data, "An single ENDF file containing all decay data")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-i,--input", input, "Input file for inventory simulation")
      ->required()
      ->check(CLI::ExistingFile);
    CLI11_PARSE(app, argc, argv);

    try{
      // construct a decay matrix - straight from ENDF file(s)
      const Eigen::MatrixX<double> dense_decay_mat = endf::construct_decay_matrix_from_file(data);

      // read input data
      const auto input_data = io::helpers::from_file<double, io::input::InputData>(input);

      // do computation
      const auto output_data = solver::compute_inventory<double>(input_data, dense_decay_mat);

      // dump output to std::out
      std::cout << output_data;
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
