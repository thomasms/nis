#include <iostream>
#include <map>
#include <cmath>

#include <CLI/CLI.hpp>
#include <ENDFtk.hpp>
#include <Eigen/Dense>

#include "core.hpp"
#include "util.hpp"

using namespace nis;

int main(int argc, char **argv){
    CLI::App app{"NIS"};

    std::string filename = "default";
    int nuclide_zai = 0;
    app.add_option("-f,--file", filename, "An ENDF file")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-n,--nuclide", nuclide_zai, "Nuclide to examine")
      ->required();

    CLI11_PARSE(app, argc, argv);

    auto tape = njoy::ENDFtk::tree::fromFile(filename);

    // construct a decay matrix
    Eigen::MatrixXd decay_mat = core::construct_decay_matrix_from_tape(tape);

    // for testing....
    const auto zai_map = util::make_map_index_lookup(data::ALL_ZAIS);
    if(zai_map.count(nuclide_zai) == 1){
        const auto zai_index = zai_map.at(nuclide_zai);
        int vec_size = static_cast<int>(decay_mat.cols());
        // let's examine one nuclide
        Eigen::VectorXd vec = Eigen::VectorXd::Zero(vec_size);
        vec(zai_index) = 1;
        Eigen::VectorXd result = decay_mat.transpose()*vec;
        for(int i=0;i<vec_size;++i){
          double val = result[i];
          if(val != 0){
            std::cout << nuclide_zai << " -> " << data::ALL_ZAIS[i] << " (" << val << ")\n";
          }
        }
    }

    // provide an input inventory
    // ...

    // solve the system for a given time
    // ...

    // output the inventory

    return 0;
}
