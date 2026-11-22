#ifndef DATA_DECAY_HPP
#define DATA_DECAY_HPP

#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "nis/decayheat.hpp"
#include "nis/physics.hpp"
#include "nis/data/decay_products.hpp"
#include "nis/data/fission_yields.hpp"

// The decay data, with no dependency on the file format.
// You can read it from ENDF (nis/endf/decay.hpp), or make and change it by hand,
// then make the decay matrix and the decay energy data from it.

namespace nis {
    namespace data {

        // The rules that give the products of each decay mode (RTYP), with the ENDF numbers:
        // 0 gamma, 1 beta-, 2 beta+ (and EC), 3 IT, 4 alpha, 5 neutron, 6 spontaneous fission,
        // 7 proton, 8 beta- (other), 9 gamma (other), 10 unknown (goes to the sink).
        const std::map<int, DecayProductComputer> DECAY_MAP = {
            {0, GammaRayReactionProducts() },
            {1, BetaMinusDecayReactionProducts() },
            {2, BetaPlusDecayReactionProducts() },
            {3, ITDecayReactionProducts() },
            {4, AlphaDecayReactionProducts() },
            {5, NeutronEmissionDecayReactionProducts() },
            {6, SFDecayReactionProducts() },
            {7, ProtonEmissionDecayReactionProducts() },
            {8, BetaMinusDecayReactionProducts() },
            {9, GammaRayReactionProducts() },
            // unknown origin
            {10, ToSinkDecayReactionProducts() }
        };

        // An RTYP can be more than one decay, one after the other: 1.5 is beta- then neutron emission.
        std::vector<int> expand_rtypes(double rtyp);

        DecayProductT compute_products(int parent_zai, double rtyp, int final_isomeric_state);

        // spontaneous fission is decay mode 6 (RTYP)
        constexpr int SPONTANEOUS_FISSION_RTYP = 6;

        // one decay mode of a nuclide
        struct DecayMode {
            // the decay type, with the ENDF numbers (see DECAY_MAP), i.e. 1 is beta-, 1.5 is beta- then neutron
            double rtyp = 0.0;
            // the fraction of the decays that use this mode
            double branching_ratio = 0.0;
            // the isomeric state of the daughter
            int final_state = 0;
        };

        // the decay data of one nuclide
        struct DecayNuclide {
            // in s, infinity (or 0 or less) is stable
            double half_life = std::numeric_limits<double>::infinity();
            std::vector<DecayMode> modes;
            // average energy per decay in eV, for decay heat (see core::DecayEnergyData)
            double light_particle = 0.0;
            double electromagnetic = 0.0;
            double heavy_particle = 0.0;

            bool is_stable() const {
                return !((half_life > 0.0) && std::isfinite(half_life));
            }

            // 1/s, 0 for a stable nuclide
            double decay_constant() const {
                return is_stable() ? 0.0 : physics::compute_decay_const(half_life);
            }

            // a decay constant of 0 or less makes the nuclide stable
            void set_decay_constant(double decay_const){
                half_life = (decay_const > 0.0) ? physics::LN_2/decay_const : std::numeric_limits<double>::infinity();
            }

            // the sum of the branching ratios, 1 for good data
            double total_branching_ratio() const {
                double total = 0.0;
                for(const auto& mode: modes){
                    total += mode.branching_ratio;
                }
                return total;
            }
        };

        // by ZAI. Nuclides that are not in ZAI_MAP are not in the matrix.
        using DecayData = std::map<int, DecayNuclide>;

        // The matrix is sparse, (parent, daughter), with the sink as the last index.
        // With sf_yields, spontaneous fission makes the fission products, at the partial decay rate
        // times the independent yield. Without them, or for a nuclide with no yields,
        // spontaneous fission goes to the sink.
        Eigen::SparseMatrix<double> construct_decay_matrix(const DecayData& decay_data,
                                                           const FissionYieldData* sf_yields = nullptr);

        // decay constants and average decay energies for decay heat, stable nuclides are zero
        core::DecayEnergyData<double> construct_decay_energy_data(const DecayData& decay_data);

    }
}

#endif // DATA_DECAY_HPP
