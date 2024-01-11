#include <iostream>
#include <map>
#include <cmath>

#include <CLI/CLI.hpp>
#include <ENDFtk.hpp>
#include <Eigen/Dense>

#include "nuclides.hpp"
#include "utils.hpp"

using namespace nis;

const double ln_2 = std::log(2);

int main(int argc, char **argv){
  CLI::App app{"NIS"};

  std::string filename = "default";
  app.add_option("-f,--file", filename, "An ENDF file")
    ->required()
    ->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);

  auto tape = njoy::ENDFtk::tree::fromFile(filename);
  const std::vector<int> mat_numbers = tape.materialNumbers();

  // this could be done at compile time using std::array
  const auto zai_map = util::make_map_index_lookup(data::ALL_ZAIS);

  // loop over all materials in the file, most likely just one
	for(auto mat_number: mat_numbers){
    // this can return multiple instances of the same MAT number
    // due to different temperatures
    // loop over temperatures for each MAT
    auto materials = tape.material(mat_number);
    for(auto material: materials){

      // add 1 for the sink
      const int mat_nr_rows = zai_map.size() + 1;
      const int mat_nr_cols = zai_map.size() + 1;
      Eigen::MatrixXd decay_mat(mat_nr_rows, mat_nr_cols);

      // loop over file numbers for the given material
      auto file_numbers = material.fileNumbers();
      for(auto file_number: file_numbers){
        // we only want file MF=8 for decay data
        if(file_number == 8){
          njoy::ENDFtk::tree::File file = material.file(file_number);
          // we want MT=457 for decay data
          if(file.hasMT(457)){
            njoy::ENDFtk::tree::Section decay_data = file.MT(457);
            auto parsed_section = decay_data.parse<8, 457>();
            const int ZA = parsed_section.ZA();
            const int za_index = zai_map.at(ZA);
            
            if( parsed_section.isStable() ){
              decay_mat(za_index, za_index) = 0;
            }
            else{
              const double half_life = parsed_section.halfLife()[0];
              const double decay_const = ln_2/half_life;
              decay_mat(za_index, za_index) = -decay_const;
              std::cout << ZA << " " << za_index << " " << decay_const << "\n";
              // todo: loop over decay states and compute the outgoing
            }
          }
        }
      }
      // std::cout << decay_mat;
    }
  }

  return 0;
}
