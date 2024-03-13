#include <iostream>
#include <map>
#include <cmath>

#include <CLI/CLI.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "nis/cram.hpp"
#include "nis/nuclides.hpp"
#include "nis/util.hpp"
#include "nis/zai.hpp"
#include "nis/io/input.hpp"
#include "nis/endf/core.hpp"

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
      const auto input_data = io::input::from_file<double>(input);

      // set initial vector (inventory)
      const int vec_size = static_cast<int>(dense_decay_mat.cols());
      Eigen::VectorX<double> n0 = Eigen::VectorX<double>::Zero(vec_size);
      const auto zai_map = util::make_map_index_lookup(data::ALL_ZAIS);

      for (auto it = input_data.nuclides.begin(); it != input_data.nuclides.end(); ++it) {
        const int nuclide_zai = zai::get_from_name(it->first);
        if(zai_map.count(nuclide_zai) == 1){
            const auto zai_index = zai_map.at(nuclide_zai);
            n0(zai_index) = it->second;
        }
      }

      const Eigen::VectorX<double> nt = core::cram48::solve<double>(dense_decay_mat, n0, input_data.dt);

      // output
      for(int i=0;i<vec_size;++i){
        double val = nt[i];
        int zai = (i == vec_size-1) ? -1 : data::ALL_ZAIS[i];
        std::string nuclide_name = zai::get_name_from(zai);
        if(val > 1e4){
          std::cout << nuclide_name << " " << val << "\n";
        }
      }
    }
    catch(...){
      std::cerr << "There was an error - check the logs.\n";
      return 1;
    }

    return 0;
}
