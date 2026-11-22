#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>

#include "nis/maths.hpp"
#include "nis/nuclides.hpp"
#include "nis/spectrum.hpp"
#include "nis/zai.hpp"

#include "nis/data/reaction.hpp"

namespace nis {
    namespace data {

        bool is_reaction_mt(int mt_number){
            return (mt_number == FISSION_MT) || (REACTION_MAP.count(mt_number) == 1);
        }

        bool is_light_particle(int za){
            return (za == 1) || (za == 1001) || (za == 1002) || (za == 1003) || (za == 2003) || (za == 2004);
        }

        ReactionProductT compute_reaction_products(int parent_zai, int mt_number, int final_state){
            if(mt_number == FISSION_MT){
                // the fission products come from the yields (construct_reaction_matrix), not from here
                return ReactionProductT(-1, {});
            }

            if(REACTION_MAP.count(mt_number) == 0){
                throw std::invalid_argument("MT=" + std::to_string(mt_number) + " is not a supported reaction");
            }

            const ReactionChange& change = REACTION_MAP.at(mt_number);
            const zai::triple parent_triple = zai::decompose(parent_zai);
            const int z = std::get<0>(parent_triple) + change.dz;
            const int a = std::get<1>(parent_triple) + change.da;

            // unphysical residual, i.e. (n,a) on a very light target
            if((z < 0) || (a < z) || (a <= 0)){
                return ReactionProductT(-1, change.secondaries);
            }

            return ReactionProductT(zai::make(z, a, final_state), change.secondaries);
        }

        std::map<int, double> split_rate_by_state(double total_rate, const std::map<int, double>& state_rates){
            double states_sum = 0.;
            for(const auto& state_rate: state_rates){
                states_sum += state_rate.second;
            }

            // no split, all to the ground state
            if(states_sum <= 0.){
                return {{0, total_rate}};
            }

            std::map<int, double> split;
            const bool has_ground = state_rates.count(0) == 1;
            if(!has_ground && (states_sum <= total_rate)){
                // the isomers as given, the rest to the ground state
                split = state_rates;
                split[0] = total_rate - states_sum;
                return split;
            }

            // use MF=10 for the fractions only, so the total stays the same
            for(const auto& state_rate: state_rates){
                split[state_rate.first] = total_rate*state_rate.second/states_sum;
            }
            return split;
        }

        void _check_cross_section(const std::vector<double>& xs_energies, const std::vector<double>& xs_values){
            if((xs_energies.size() < 2) || (xs_energies.size() != xs_values.size())){
                throw std::invalid_argument("Cross section must have at least 2 points and the same number of energies and values");
            }
        }

        // the flux of each group of the cross section, dotted with its values
        double _collapse_on_groups(const std::vector<double>& group_flux, const std::vector<double>& xs_values){
            // because it is INT=1 (constant) the last value has no bin, so we drop it
            const std::vector<double> group_values(xs_values.begin(), xs_values.end()-1);
            return maths::dot_product<double>(group_flux, group_values)*BARN_TO_CM2;
        }

        double collapse_cross_section(const core::Spectrum<double>& in_spectrum,
                                      const std::vector<double>& xs_energies,
                                      const std::vector<double>& xs_values){
            _check_cross_section(xs_energies, xs_values);
            const core::Spectrum<double> data_spectrum = core::convert_by_lethargy<double>(in_spectrum, xs_energies);
            return _collapse_on_groups(data_spectrum.y(), xs_values);
        }

        void add_reaction_rate(Eigen::SparseMatrix<double>& reaction_mat, const std::map<int, int>& zai_map,
                               int parent_zai, int mt_number, double rate, int final_state){
            if((zai_map.count(parent_zai) == 0) || (rate <= 0.)){
                return;
            }
            const int za_index = zai_map.at(parent_zai);

            // the diagonal is the rate lost from the parent
            reaction_mat.coeffRef(za_index, za_index) -= rate;

            const ReactionProductT products = compute_reaction_products(parent_zai, mt_number, final_state);
            std::vector<int> daughters = std::get<1>(products);
            daughters.emplace_back(std::get<0>(products));
            for(auto daughter_zai: daughters){
                // the sink is the last row
                const int dza_index = zai_map.count(daughter_zai) == 1 ? zai_map.at(daughter_zai) : SINK_INDEX;
                reaction_mat.coeffRef(za_index, dza_index) += rate;
            }
        }

        // adds the rate from the parent to the product (or the sink), and the loss of the parent if asked
        void _add_transfer(Eigen::SparseMatrix<double>& reaction_mat, int parent_zai, int product_zai,
                           double rate, bool parent_loss){
            if((ZAI_MAP.count(parent_zai) == 0) || (rate <= 0.)){
                return;
            }
            const int za_index = ZAI_MAP.at(parent_zai);
            if(parent_loss){
                reaction_mat.coeffRef(za_index, za_index) -= rate;
            }
            const int dza_index = ZAI_MAP.count(product_zai) == 1 ? ZAI_MAP.at(product_zai) : SINK_INDEX;
            reaction_mat.coeffRef(za_index, dza_index) += rate;
        }

        Eigen::SparseMatrix<double> construct_reaction_matrix(const ReactionData& reaction_data,
                                                  const core::Spectrum<double>& in_spectrum,
                                                  const FissionYieldData* fission_yields){
            // collapse all first, then split each reaction by residual state
            struct Rates {
                double total = 0.;
                std::map<int, double> states;
            };
            std::map<std::pair<int, int>, Rates> rates;
            // (parent, product ZAI) for MT=4 and MT=5, these come from MF=10 only
            std::map<std::pair<int, int>, double> inelastic_rates;
            std::map<std::pair<int, int>, double> anything_rates;
            // fission of a parent with yields: the fission rate in each group of the cross section
            struct GroupRates {
                std::vector<double> bounds;
                std::vector<double> rates;
            };
            std::map<int, GroupRates> fission_rates;

            // the data has only a few energy grids, so we convert the spectrum once for each
            std::map<std::vector<double>, std::vector<double>> grid_flux;
            const auto collapse = [&](const CrossSection& xs){
                _check_cross_section(xs.energies, xs.values);
                auto it = grid_flux.find(xs.energies);
                if(it == grid_flux.end()){
                    it = grid_flux.emplace(xs.energies, core::convert_by_lethargy<double>(in_spectrum, xs.energies).ycopy()).first;
                }
                return _collapse_on_groups(it->second, xs.values);
            };

            for(const auto& xs: reaction_data){
                if(xs.mt_number == INELASTIC_MT){
                    if(xs.final_state > 0){
                        const int product = zai::make(zai::get_z(xs.parent_zai), zai::get_a(xs.parent_zai), xs.final_state);
                        inelastic_rates[{xs.parent_zai, product}] += collapse(xs);
                    }
                    continue;
                }
                if(xs.mt_number == ANYTHING_MT){
                    if(xs.product_za > 0){
                        const int product = xs.product_za*10 + std::max(xs.final_state, 0);
                        anything_rates[{xs.parent_zai, product}] += collapse(xs);
                    }
                    continue;
                }
                // the data can be made by hand, so check again
                if(!is_reaction_mt(xs.mt_number)){
                    continue;
                }
                if((xs.mt_number == FISSION_MT) && (fission_yields != nullptr)
                   && (fission_yields->count(xs.parent_zai) == 1) && (xs.final_state == ALL_STATES)){
                    collapse(xs);
                    const std::vector<double>& flux = grid_flux.at(xs.energies);
                    GroupRates& f = fission_rates[xs.parent_zai];
                    f.bounds = xs.energies;
                    f.rates.assign(flux.size(), 0.0);
                    for(std::size_t g = 0; g < flux.size(); ++g){
                        f.rates[g] = flux[g]*xs.values[g]*BARN_TO_CM2;
                    }
                    continue;
                }
                const double rate = collapse(xs);
                Rates& r = rates[{xs.parent_zai, xs.mt_number}];
                if(xs.final_state == ALL_STATES){
                    r.total += rate;
                }
                else{
                    r.states[xs.final_state] += rate;
                }
            }

            Eigen::SparseMatrix<double> reaction_mat(ZAI_SIZE, ZAI_SIZE);
            for(const auto& [key, r]: rates){
                const int parent_zai = key.first;
                const int mt_number = key.second;
                for(const auto& [state, rate]: split_rate_by_state(r.total, r.states)){
                    // if we don't track this isomer, the atoms go to the ground state
                    const int residual = std::get<0>(compute_reaction_products(parent_zai, mt_number, state));
                    const int final_state = (residual >= 0) && (ZAI_MAP.count(residual) == 0) ? 0 : state;
                    add_reaction_rate(reaction_mat, ZAI_MAP, parent_zai, mt_number, rate, final_state);
                }
            }

            // fission: the parent is lost at the fission rate, each product is made at the rate times its yield
            for(const auto& [parent_zai, f]: fission_rates){
                double total = 0.;
                for(const auto r: f.rates){
                    total += r;
                }
                if((total <= 0.) || (ZAI_MAP.count(parent_zai) == 0)){
                    continue;
                }
                const int parent_index = ZAI_MAP.at(parent_zai);
                reaction_mat.coeffRef(parent_index, parent_index) -= total;
                for(const auto& [product, y]: weighted_yields(fission_yields->at(parent_zai), f.bounds, f.rates)){
                    // if we don't track this isomer, the atoms go to the ground state
                    const int tracked = ZAI_MAP.count(product) == 1 ? product : (product/10)*10;
                    _add_transfer(reaction_mat, parent_zai, tracked, total*y, false);
                }
            }

            // (n,n') to an isomer we don't track gives the target again, so there is no change
            for(const auto& [key, rate]: inelastic_rates){
                if(ZAI_MAP.count(key.second) == 1){
                    _add_transfer(reaction_mat, key.first, key.second, rate, true);
                }
            }

            // MT=5: each residual is one reaction, so the loss of the parent is the sum of the residuals.
            // MF=3 MT=5 is not used, as it does not always agree with the sum in MF=10.
            // The emitted light particles are added, but are not a loss.
            for(const auto& [key, rate]: anything_rates){
                const int product = key.second;
                const bool residual = !is_light_particle(product/10);
                // the free neutron is not tracked
                if(product/10 == 1){
                    continue;
                }
                // if we don't track this isomer, the atoms go to the ground state
                const int tracked = ZAI_MAP.count(product) == 1 ? product : (product/10)*10;
                _add_transfer(reaction_mat, key.first, tracked, rate, residual);
            }
            reaction_mat.makeCompressed();
            return reaction_mat;
        }

    }
}
