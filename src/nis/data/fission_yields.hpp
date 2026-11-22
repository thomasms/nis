#ifndef DATA_FISSION_YIELDS_HPP
#define DATA_FISSION_YIELDS_HPP

#include <map>
#include <vector>

// The fission yields, with no dependency on the file format.
// You can read them from ENDF (nis/endf/fission_yields.hpp), or make and change them by hand.

namespace nis {
    namespace data {

        // the yields of one fissile nuclide, a set of yields at each incident neutron energy
        struct FissionYields {
            int parent_zai = 0;
            // incident neutron energies in eV, in order. A spontaneous fission file has one set, at 0 eV.
            std::vector<double> energies;
            // one set for each energy: product ZAI -> independent yield (products per fission)
            std::vector<std::map<int, double>> yields;
        };

        // by parent ZAI
        using FissionYieldData = std::map<int, FissionYields>;

        // The weight of each yield energy for a fission at this energy: linear interpolation in energy
        // between the 2 energies around it, and the nearest set outside the energies.
        std::vector<double> energy_weights(const std::vector<double>& energies, double energy);

        // the yields for a fission at this energy
        std::map<int, double> yields_at_energy(const FissionYields& yields, double energy);

        // The yields for fissions spread over energy groups: group_rates[g] is the fission rate in the group
        // from bounds[g] to bounds[g+1] (eV). Each group uses the yields at its mid energy (in lethargy),
        // and the groups are weighted by their fission rate. With no fissions, the yields at the first energy.
        std::map<int, double> weighted_yields(const FissionYields& yields,
                                              const std::vector<double>& bounds,
                                              const std::vector<double>& group_rates);
    }
}

#endif // DATA_FISSION_YIELDS_HPP
