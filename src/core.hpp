#ifndef CORE_HPP
#define CORE_HPP

#include <cmath>
#include <memory>

#include <ENDFtk.hpp>
#include <Eigen/Dense>

#include "physics.hpp"
#include "nuclides.hpp"
#include "util.hpp"
#include "zai.hpp"


namespace nis {
    namespace core {

        // (primary, {secondaries})
        using DecayProductT = std::tuple<int, std::vector<int>>;

        struct DecayProductComputer{
            virtual DecayProductT get_products(int parent_zai) const = 0;
            virtual ~DecayProductComputer() {};
        };

        struct GammaRayReactionProducts : public DecayProductComputer{

            DecayProductT get_products(int parent_zai) const override{
                return DecayProductT(parent_zai, {});
            };

        };

        struct BetaMinusDecayReactionProducts : public DecayProductComputer{

            DecayProductT get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // Z increases by one, as a neutron goes to a proton
                // thus A remains the same.
                // i.e. C14 (60140) -> N14 (70140) + e- + mu_e
                const int daughter_zai = zai::make(std::get<0>(parent_triple)+1,
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                // we don't care about the electron or neutrino
                return DecayProductT(daughter_zai, {});
            };
        };

        struct BetaPlusDecayReactionProducts : public DecayProductComputer{

            DecayProductT get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // Z decreases by one, as a proton goes to a neutron
                // thus A remains the same.
                // i.e. A26 (130260) + e- -> Mg26 (120260) + mu_e
                const int daughter_zai = zai::make(std::get<0>(parent_triple)-1,
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                // we don't care about the electron or neutrino
                return DecayProductT(daughter_zai, {});
            };
        };

        struct ITDecayReactionProducts : public DecayProductComputer{

            DecayProductT get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // TODO: need to deal with isometric state
                const int daughter_zai = zai::make(std::get<0>(parent_triple),
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                return DecayProductT(daughter_zai, {});
            };
        };

        struct AlphaDecayReactionProducts : public DecayProductComputer{

            DecayProductT get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                const int daughter_zai = zai::make(std::get<0>(parent_triple)-2,
                                                   std::get<1>(parent_triple)-4,
                                                   std::get<2>(parent_triple));
                return DecayProductT(daughter_zai, {zai::ALPHA});
            };
        };

        struct NeutronEmissionDecayReactionProducts : public DecayProductComputer{

            DecayProductT get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                const int daughter_zai = zai::make(std::get<0>(parent_triple),
                                                   std::get<1>(parent_triple)-1,
                                                   std::get<2>(parent_triple));
                return DecayProductT(daughter_zai, {});
            };
        };

        // this one is tricky :) - we need more information/data from MT=454,459
        struct SFDecayReactionProducts : public DecayProductComputer{

            DecayProductT get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                // todo: ....
                return DecayProductT(-1, {});
            };
        };

        struct ProtonEmissionDecayReactionProducts : public DecayProductComputer{

            DecayProductT get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                const int daughter_zai = zai::make(std::get<0>(parent_triple)-1,
                                                   std::get<1>(parent_triple),
                                                   std::get<2>(parent_triple));
                return DecayProductT(daughter_zai, {});
            };
        };

        struct ToSinkDecayReactionProducts : public DecayProductComputer{

            DecayProductT get_products(int parent_zai) const override{
                const zai::triple parent_triple = zai::decompose(parent_zai);
                return DecayProductT(-1, {});
            };
        };

        const std::map<int, std::function<std::unique_ptr<DecayProductComputer>()>> DECAY_MAP = {
            {0, []() { return std::make_unique<GammaRayReactionProducts>(); }},
            {1, []() { return std::make_unique<BetaMinusDecayReactionProducts>(); }},
            {2, []() { return std::make_unique<BetaPlusDecayReactionProducts>(); }},
            {3, []() { return std::make_unique<ITDecayReactionProducts>(); }},
            {4, []() { return std::make_unique<AlphaDecayReactionProducts>(); }},
            {5, []() { return std::make_unique<NeutronEmissionDecayReactionProducts>(); }},
            {6, []() { return std::make_unique<SFDecayReactionProducts>(); }},
            {7, []() { return std::make_unique<ProtonEmissionDecayReactionProducts>(); }},
            {8, []() { return std::make_unique<BetaMinusDecayReactionProducts>(); }},
            {9, []() { return std::make_unique<GammaRayReactionProducts>(); }},
            // unknown origin
            {10, []() { return std::make_unique<ToSinkDecayReactionProducts>(); }}
        };

        std::vector<int> _expand_rtypes(double rtyp){
            // rtype is actually a double which can contain multiple reaction modes
            // i.e. 1.5 is beta minus followed by neutron emission
            // We assume that at most we can have 7 subsequent reactions, hence the 1e6 factor.
            // Since this is stored in the ENDF file, with limited precision, it should hold

            // if we are zero first then we are gamma
            const int initial_rtype = std::floor(rtyp);
            if(initial_rtype == 0){
                return std::vector<int>(0);
            }

            const int expanded_rtype = static_cast<int>(rtyp*1'000'000);
            std::list<int> rtypes = util::get_digits(expanded_rtype);
            // this is to account for initial procedures bigger than 9
            rtypes.pop_front();
            rtypes.push_front(initial_rtype);
            // we also need to remove trailing zeros, as this does not mean GammaRay emission
            // it is the end
            std::vector<int> reaction_types;
            for(int v: rtypes){
                if((v == 0) || (DECAY_MAP.count(v) == 0)){
                    break;
                }
                reaction_types.push_back(v);
            }
            return reaction_types;
        }

        DecayProductT compute_products(int parent_zai, double rtyp, int final_isomeric_state){
            using Comp_ptr_t = std::unique_ptr<DecayProductComputer>;
            const std::vector<int> rtypes = _expand_rtypes(rtyp);

            // we need to trim remaining zeros, as there is no way to tell between
            // 1.0 (beta -> gamma) and 1.0 (beta)
            Comp_ptr_t dpc = std::make_unique<ToSinkDecayReactionProducts>();
            int primary_zai = parent_zai;
            std::vector<int> secondaries;
            for(auto reaction_type: rtypes){
                if (DECAY_MAP.count(reaction_type) == 1){
                    auto create_computer = DECAY_MAP.at(reaction_type);
                    dpc = create_computer();
                    // we assume - and because the data tells us this -
                    // that SF (6) doesn't have any modes following it
                    // to make this simpler with the primary
                    // Perhaps a tree is a better data structure here....
                    DecayProductT prods = dpc->get_products(primary_zai);
                    primary_zai = std::get<0>(prods);
                    for(auto secondary: std::get<1>(prods)){
                        secondaries.push_back(secondary);
                    }
                }
            }
            // finally check the final isomeric state
            const zai::triple primary_triple = zai::decompose(primary_zai);
            primary_zai = zai::make(std::get<0>(primary_triple),
                                    std::get<1>(primary_triple),
                                    final_isomeric_state);

            return DecayProductT(primary_zai, secondaries);
        }


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

#endif //CORE_HPP