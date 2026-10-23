#include <cmath>
#include <memory>

#include "nis/physics.hpp"
#include "nis/nuclides.hpp"
#include "nis/util.hpp"
#include "nis/zai.hpp"

#include "nis/endf/decay.hpp"


namespace nis {
    namespace endf {

        std::vector<int> expand_rtypes(double rtyp){
            // rtype is actually a double which can contain multiple reaction modes
            // i.e. 1.5 is beta minus followed by neutron emission
            // We assume that at most we can have 7 subsequent reactions, hence the 1e6 factor.
            // Since this is stored in the ENDF file, with limited precision, it should hold
            constexpr int _rtyp_factor = 1'000'000;
            const int initial_rtype = static_cast<int>(std::floor(rtyp));

            // if we are zero first then we are gamma
            if(initial_rtype == 0){
                return std::vector<int>{0};
            }

            // if the initial rtype is not valid
            // then return an empty vector
            // regardless of what else comes next
            if(DECAY_MAP.count(initial_rtype) == 0){
                return std::vector<int>();
            }

            // all other cases we process the subsequent values
            const int rtyp_as_int = static_cast<int>(rtyp*_rtyp_factor);
            const int expanded_rtype = rtyp_as_int - (initial_rtype*_rtyp_factor);
            std::list<int> rtypes = util::get_digits(expanded_rtype);
            // we also need to remove trailing zeros,
            // as this does not mean GammaRay emission
            // it is the end
            std::vector<int> reaction_types{initial_rtype};
            for(int v: rtypes){
                if((v == 0) || (DECAY_MAP.count(v) == 0)){
                    break;
                }
                reaction_types.emplace_back(v);
            }
            return reaction_types;
        }

        DecayProductT compute_products(int parent_zai, double rtyp, int final_isomeric_state){
            const std::vector<int> rtypes = expand_rtypes(rtyp);

            // we need to trim remaining zeros, as there is no way to tell between
            // 1.0 (beta -> gamma) and 1.0 (beta)
            DecayProductComputer dpc = ToSinkDecayReactionProducts();
            int primary_zai = parent_zai;
            std::vector<int> secondaries;
            for(auto reaction_type: rtypes){
                if (DECAY_MAP.count(reaction_type) == 1){
                    dpc = DECAY_MAP.at(reaction_type);
                    // we assume - and because the data tells us this -
                    // that SF (6) doesn't have any modes following it
                    // to make this simpler with the primary
                    // Perhaps a tree is a better data structure here....
                    DecayProductT prods = std::visit([&primary_zai](auto const& obj) {
                        return obj.get_products(primary_zai);
                    }, dpc);
                    primary_zai = std::get<0>(prods);
                    for(auto& secondary: std::get<1>(prods)){
                        secondaries.emplace_back(secondary);
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

        constexpr int file_number = 8;
        constexpr int mt_number = 457;

        int _get_zai_from_section(const njoy::ENDFtk::section::Type<file_number, mt_number>& section){
            const int za = section.ZA();
            // ENDFtk gives us the isomeric state as an unsigned int
            const int i = static_cast<int>(section.LISO());
            const int zai = za*10 + i;
            return zai;
        }

        std::vector<int> _get_secondaries(int zai, const njoy::ENDFtk::section::Type<file_number, mt_number>::DecayMode& decay_mode){
            const double reaction_type = decay_mode.RTYP();
            const int final_isomeric_state = std::floor(decay_mode.finalIsomericState());
            const auto products = compute_products(zai, reaction_type, final_isomeric_state);
            const int primary = std::get<0>(products);

            std::vector<int> secondaries = std::get<1>(products);
            // add the primary as the last element
            secondaries.emplace_back(primary);
            return secondaries;
        }

        void _update_decay_matrix(const njoy::ENDFtk::tree::Material& material,
                                  const std::map<int, int>& zai_map, Eigen::SparseMatrix<double>& decay_mat){

            if(material.hasMF(file_number)){
                njoy::ENDFtk::tree::File file = material.file(file_number);
                if(file.hasMT(mt_number)){

                    njoy::ENDFtk::tree::Section decay_data = file.MT(mt_number);
                    auto parsed_section = decay_data.parse<file_number, mt_number>();
                    const int zai = _get_zai_from_section(parsed_section);

                    if (zai_map.count(zai) == 1){
                        const int za_index = zai_map.at(zai);
                        const double half_life = parsed_section.halfLife()[0];

                        if( (!parsed_section.isStable()) && (half_life > 0.) ){
                            const double decay_const = physics::compute_decay_const(half_life);

                            // the diagonal is the rate lost from the parent
                            decay_mat.coeffRef(za_index, za_index) -= decay_const;

                            // loop over decay states and compute the outgoing
                            const auto decay_modes_collection = parsed_section.decayModes();
                            for(auto decay_mode: decay_modes_collection.decayModes()){
                                const double br = decay_mode.branchingRatio()[0];
                                const auto secondaries = _get_secondaries(zai, decay_mode);
                                for(auto daughter_zai: secondaries){
                                    // the sink is the last row
                                    int dza_index = zai_map.count(daughter_zai) == 1 ? zai_map.at(daughter_zai) : data::SINK_INDEX;
                                    decay_mat.coeffRef(za_index, dza_index) += br*decay_const;
                                }
                            }
                        }
                    }
                }
            }
        }

        Eigen::SparseMatrix<double> construct_decay_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape){

            // add 1 for the sink
            const int mat_nr_rows = data::ZAI_SIZE;
            const int mat_nr_cols = data::ZAI_SIZE;

            Eigen::SparseMatrix<double> decay_mat(mat_nr_rows, mat_nr_cols);

            // loop over all materials in the file, most likely just one
            for(auto mat_number: tape.materialNumbers()){
                // this can return multiple instances of the same MAT number
                // due to different temperatures
                // loop over temperatures for each MAT
                //
                // How do we account for multiple temperatures?
                // We don't want to double count...
                auto materials = tape.material(mat_number);
                for(auto material: materials){
                    _update_decay_matrix(material, data::ZAI_MAP, decay_mat);
                }
            }

            decay_mat.makeCompressed();
            return decay_mat;
        }

        Eigen::SparseMatrix<double> construct_decay_matrix_from_file(std::string filename){
            auto tape = njoy::ENDFtk::tree::fromFile(filename);
            return endf::construct_decay_matrix_from_tape(tape);
        }

        void _update_decay_energy_data(const njoy::ENDFtk::tree::Material& material,
                                       const std::map<int, int>& zai_map,
                                       core::DecayEnergyData<double>& energy_data){

            if(!material.hasMF(file_number)){
                return;
            }
            njoy::ENDFtk::tree::File file = material.file(file_number);
            if(!file.hasMT(mt_number)){
                return;
            }

            njoy::ENDFtk::tree::Section decay_data = file.MT(mt_number);
            auto parsed_section = decay_data.parse<file_number, mt_number>();
            const int zai = _get_zai_from_section(parsed_section);
            if(zai_map.count(zai) == 0){
                return;
            }

            const double half_life = parsed_section.halfLife()[0];
            // same rule as the decay matrix, so both use the same nuclides
            if(parsed_section.isStable() || (half_life <= 0.)){
                return;
            }

            const int za_index = zai_map.at(zai);
            const auto& energies = parsed_section.averageDecayEnergies();
            energy_data.decay_constants(za_index) = physics::compute_decay_const(half_life);
            energy_data.light_particle(za_index) = energies.ELP()[0];
            energy_data.electromagnetic(za_index) = energies.EEM()[0];
            energy_data.heavy_particle(za_index) = energies.EHP()[0];
        }

        core::DecayEnergyData<double> read_decay_energy_data_from_tape(const njoy::ENDFtk::tree::Tape& tape){
            core::DecayEnergyData<double> energy_data;
            for(auto mat_number: tape.materialNumbers()){
                for(auto material: tape.material(mat_number)){
                    _update_decay_energy_data(material, data::ZAI_MAP, energy_data);
                }
            }
            return energy_data;
        }

        core::DecayEnergyData<double> read_decay_energy_data_from_file(std::string filename){
            auto tape = njoy::ENDFtk::tree::fromFile(filename);
            return endf::read_decay_energy_data_from_tape(tape);
        }


    }
}
