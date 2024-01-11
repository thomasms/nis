#include <iostream>
#include <map>
#include <cmath>

#include <CLI/CLI.hpp>
#include <ENDFtk.hpp>
#include <Eigen/Dense>

#include "core.hpp"

using namespace nis;

int main(int argc, char **argv){
    CLI::App app{"NIS"};

    std::string filename = "default";
    app.add_option("-f,--file", filename, "An ENDF file")
      ->required()
      ->check(CLI::ExistingFile);

    CLI11_PARSE(app, argc, argv);

    auto tape = njoy::ENDFtk::tree::fromFile(filename);
    Eigen::MatrixXd decay_mat = core::construct_decay_matrix_from_tape(tape);

    return 0;
}
