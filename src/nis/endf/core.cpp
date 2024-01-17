#include <cmath>
#include <memory>

#include "nis/physics.hpp"
#include "nis/nuclides.hpp"
#include "nis/util.hpp"

#include "nis/endf/core.hpp"
#include "nis/endf/decay.hpp"

namespace nis {
    namespace endf {

        /*
        Needs to be cleaned up and move somewhere else.
        */
        void _update_decay_matrix(const njoy::ENDFtk::tree::Material& material,
                                  const std::map<int, int>& zai_map, Eigen::MatrixXd& decay_mat){
            constexpr int file_number = 8;
            constexpr int mt_number = 457;
            if(material.hasMF(file_number)){
                njoy::ENDFtk::tree::File file = material.file(file_number);
                if(file.hasMT(mt_number)){
                    njoy::ENDFtk::tree::Section decay_data = file.MT(mt_number);
                    auto parsed_section = decay_data.parse<file_number, mt_number>();
                    const int za = parsed_section.ZA();
                    // ENDFtk gives us the isomeric state as an unsigned int
                    const int i = static_cast<int>(parsed_section.isomericState());
                    const int zai = za*10 + i;
                    if (zai_map.count(zai) == 1){
                        const int za_index = zai_map.at(zai);
        
                        if( !parsed_section.isStable() ){
                            const double half_life = parsed_section.halfLife()[0];
                            const double decay_const = physics::compute_decay_const(half_life);

                            // the diagonal is the rate lost from the parent
                            decay_mat(za_index, za_index) -= decay_const;

                            // loop over decay states and compute the outgoing
                            const auto decay_modes_collection = parsed_section.decayModes();
                            for(auto decay_mode: decay_modes_collection.decayModes()){
                                const double reaction_type = decay_mode.RTYP();
                                const int final_isomeric_state = std::floor(decay_mode.finalIsomericState());
                                const auto products = compute_products(zai, reaction_type, final_isomeric_state);
                                const int primary = std::get<0>(products);
                                std::vector<int> secondaries = std::get<1>(products);
                                secondaries.push_back(primary);
                                const double br = decay_mode.branchingRatio()[0];
                                for(auto daughter_zai: secondaries){
                                    // the sink is the last row
                                    int dza_index = decay_mat.cols() -1;
                                    if(zai_map.count(daughter_zai) == 1){
                                        dza_index = zai_map.at(daughter_zai);
                                    }
                                    decay_mat(za_index, dza_index) += br*decay_const;
                                }
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