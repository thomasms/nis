#include <ENDFtk.hpp>
#include <Eigen/Dense>

#include "nuclides.hpp"
#include "utils.hpp"

namespace nis {
    namespace core {


        const double LN_2 = std::log(2);

        /*
        Needs to be cleaned up and move somewhere else.
        */
        void _update_decay_matrix(const njoy::ENDFtk::tree::Material& material, const std::map<int, int>& zai_map, Eigen::MatrixXd& decay_mat){

            // loop over file numbers for the given material
            auto file_numbers = material.fileNumbers();
            for(auto file_number: file_numbers){
                // we only want file MF=8 for decay data
                // TODO: It would be better to use tape.MF(8) instead of this loop.
                if(file_number == 8){
                    njoy::ENDFtk::tree::File file = material.file(file_number);
                    // we want MT=457 for decay data
                    if(file.hasMT(457)){
                        njoy::ENDFtk::tree::Section decay_data = file.MT(457);
                        auto parsed_section = decay_data.parse<8, 457>();
                        const int ZA = parsed_section.ZA();
                        const int ZAI = ZA*10;
                        if (zai_map.count(ZAI) == 1){
                            const int za_index = zai_map.at(ZAI);
            
                            if( !parsed_section.isStable() ){
                                const double half_life = parsed_section.halfLife()[0];
                                const double decay_const = LN_2/half_life;
                                decay_mat(za_index, za_index) = -decay_const;
                                std::cout << ZA << " " << za_index << " " << decay_const << "\n";
                                // todo: loop over decay states and compute the outgoing
                            }
                        }
                    }
                }
            }
        }

        Eigen::MatrixXd construct_decay_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape){

            // this could be done at compile time using std::array
            const auto zai_map = util::make_map_index_lookup(data::ALL_ZAIS);
            // add 1 for the sink
            const int mat_nr_rows = zai_map.size() + 1;
            const int mat_nr_cols = mat_nr_rows;

            Eigen::MatrixXd decay_mat = Eigen::MatrixXd::Zero(mat_nr_rows, mat_nr_cols);

            const std::vector<int> mat_numbers = tape.materialNumbers();
            // loop over all materials in the file, most likely just one
            for(auto mat_number: mat_numbers){
                // this can return multiple instances of the same MAT number
                // due to different temperatures
                // loop over temperatures for each MAT
                auto materials = tape.material(mat_number);
                for(auto material: materials){
                    _update_decay_matrix(material, zai_map, decay_mat);
                }
            }

            return decay_mat;
        }


    }
}