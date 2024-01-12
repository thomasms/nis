#include <memory>

#include <ENDFtk.hpp>
#include <Eigen/Dense>

#include "physics.hpp"
#include "nuclides.hpp"
#include "util.hpp"
#include "zai.hpp"


namespace nis {
    namespace core {

        struct DecayProductComputer{
                virtual std::vector<int> get_products(int parent_zai) const = 0;
                virtual ~DecayProductComputer() {};
        };

        struct GammaRayReactionProducts : public DecayProductComputer{

            std::vector<int> get_products(int parent_zai) const override{
                return std::vector<int>{parent_zai};
            };

        };

        struct BetaMinusDecayReactionProducts : public DecayProductComputer{

            std::vector<int> get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // Z increases by one, as a neutron goes to a proton
                // thus A remains the same.
                // i.e. C14 (60140) -> N14 (70140) + e- + mu_e
                const int daughter_zai = zai::make(std::get<0>(parent_triple)+1,
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                // we don't care about the electron or neutrino
                std::vector<int> products = {daughter_zai};
                return products;
            };
        };

        struct BetaPlusDecayReactionProducts : public DecayProductComputer{

            std::vector<int> get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // Z decreases by one, as a proton goes to a neutron
                // thus A remains the same.
                // i.e. A26 (130260) + e- -> Mg26 (120260) + mu_e
                const int daughter_zai = zai::make(std::get<0>(parent_triple)-1,
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                // we don't care about the electron or neutrino
                std::vector<int> products = {daughter_zai};
                return products;
            };
        };

        struct ITDecayReactionProducts : public DecayProductComputer{

            std::vector<int> get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // TODO: need to deal with isometric state
                const int daughter_zai = zai::make(std::get<0>(parent_triple),
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                // we don't care about the electron or neutrino
                std::vector<int> products = {daughter_zai};
                return products;
            };
        };

        std::map<int, std::function<std::unique_ptr<DecayProductComputer>()>> decay_map = {
            {0, []() { return std::make_unique<GammaRayReactionProducts>(); }},
            {1, []() { return std::make_unique<BetaMinusDecayReactionProducts>(); }},
            {2, []() { return std::make_unique<BetaPlusDecayReactionProducts>(); }},
            {3, []() { return std::make_unique<ITDecayReactionProducts>(); }}
        };

        std::vector<int> compute_products(int parent_zai, int rtyp){
            if (decay_map.count(rtyp) == 1){
                auto create_computer = decay_map.at(rtyp);
                std::unique_ptr<DecayProductComputer> dpc = create_computer();
                return dpc->get_products(parent_zai);
            }
            return std::vector<int>();
        }


        /*
        Needs to be cleaned up and move somewhere else.
        */
        void _update_decay_matrix(const njoy::ENDFtk::tree::Material& material, const std::map<int, int>& zai_map, Eigen::MatrixXd& decay_mat){
            constexpr int file_number = 8;
            constexpr int mt_number = 457;
            if(material.hasMF(file_number)){
                njoy::ENDFtk::tree::File file = material.file(file_number);
                if(file.hasMT(mt_number)){
                    njoy::ENDFtk::tree::Section decay_data = file.MT(mt_number);
                    auto parsed_section = decay_data.parse<file_number, mt_number>();
                    const int za = parsed_section.ZA();
                    const int zai = za*10;
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
                                const auto products = compute_products(zai, decay_mode.RTYP());
                                const double br = decay_mode.branchingRatio()[0];
                                for(auto daughter_zai: products){
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