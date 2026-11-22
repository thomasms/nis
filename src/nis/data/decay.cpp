#include <cmath>
#include <list>

#include "nis/nuclides.hpp"
#include "nis/physics.hpp"
#include "nis/util.hpp"
#include "nis/zai.hpp"

#include "nis/data/decay.hpp"


namespace nis {
    namespace data {

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

        namespace {

            bool _is_spontaneous_fission(double rtyp){
                const std::vector<int> rtypes = expand_rtypes(rtyp);
                return !rtypes.empty() && (rtypes.front() == SPONTANEOUS_FISSION_RTYP);
            }

            // the products of the mode, with the primary as the last element
            std::vector<int> _daughters(int zai, const DecayMode& mode){
                const auto products = compute_products(zai, mode.rtyp, mode.final_state);
                std::vector<int> daughters = std::get<1>(products);
                daughters.emplace_back(std::get<0>(products));
                return daughters;
            }

            // the index of the ZAI in the matrix, or the sink if we don't track it
            int _index(int zai){
                return ZAI_MAP.count(zai) == 1 ? ZAI_MAP.at(zai) : SINK_INDEX;
            }

            void _add_nuclide(int zai, const DecayNuclide& nuclide, Eigen::SparseMatrix<double>& decay_mat,
                              const FissionYieldData* sf_yields){
                if((ZAI_MAP.count(zai) == 0) || nuclide.is_stable()){
                    return;
                }
                const int za_index = ZAI_MAP.at(zai);
                const double decay_const = nuclide.decay_constant();

                // the diagonal is the rate lost from the parent
                decay_mat.coeffRef(za_index, za_index) -= decay_const;

                // the outgoing rate of each decay mode
                for(const auto& mode: nuclide.modes){
                    const double br = mode.branching_ratio;
                    if((sf_yields != nullptr) && (sf_yields->count(zai) == 1) && _is_spontaneous_fission(mode.rtyp)){
                        // the yields at the first (lowest) energy, a spontaneous fission file has only one
                        const FissionYields& yields = sf_yields->at(zai);
                        for(const auto& [product, y]: yields_at_energy(yields, yields.energies.front())){
                            // if we don't track this isomer, the atoms go to the ground state
                            const int tracked = ZAI_MAP.count(product) == 1 ? product : (product/10)*10;
                            decay_mat.coeffRef(za_index, _index(tracked)) += br*decay_const*y;
                        }
                        continue;
                    }
                    for(const auto daughter_zai: _daughters(zai, mode)){
                        // the sink is the last row
                        decay_mat.coeffRef(za_index, _index(daughter_zai)) += br*decay_const;
                    }
                }
            }
        }

        Eigen::SparseMatrix<double> construct_decay_matrix(const DecayData& decay_data,
                                                           const FissionYieldData* sf_yields){
            // ZAI_SIZE includes the sink
            Eigen::SparseMatrix<double> decay_mat(ZAI_SIZE, ZAI_SIZE);
            for(const auto& [zai, nuclide]: decay_data){
                _add_nuclide(zai, nuclide, decay_mat, sf_yields);
            }
            decay_mat.makeCompressed();
            return decay_mat;
        }

        core::DecayEnergyData<double> construct_decay_energy_data(const DecayData& decay_data){
            core::DecayEnergyData<double> energy_data;
            for(const auto& [zai, nuclide]: decay_data){
                // same rule as the decay matrix, so both use the same nuclides
                if((ZAI_MAP.count(zai) == 0) || nuclide.is_stable()){
                    continue;
                }
                const int za_index = ZAI_MAP.at(zai);
                energy_data.decay_constants(za_index) = nuclide.decay_constant();
                energy_data.light_particle(za_index) = nuclide.light_particle;
                energy_data.electromagnetic(za_index) = nuclide.electromagnetic;
                energy_data.heavy_particle(za_index) = nuclide.heavy_particle;
            }
            return energy_data;
        }

    }
}
