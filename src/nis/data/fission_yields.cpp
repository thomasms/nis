#include <cmath>
#include <stdexcept>

#include "nis/data/fission_yields.hpp"

namespace nis {
    namespace data {

        std::vector<double> energy_weights(const std::vector<double>& energies, double energy){
            if(energies.empty()){
                throw std::invalid_argument("No fission yield energies");
            }
            std::vector<double> weights(energies.size(), 0.0);
            if(energy <= energies.front()){
                weights.front() = 1.0;
                return weights;
            }
            if(energy >= energies.back()){
                weights.back() = 1.0;
                return weights;
            }
            for(std::size_t i = 1; i < energies.size(); ++i){
                if(energy <= energies[i]){
                    const double f = (energy - energies[i-1])/(energies[i] - energies[i-1]);
                    weights[i-1] = 1.0 - f;
                    weights[i] = f;
                    break;
                }
            }
            return weights;
        }

        namespace {
            std::map<int, double> _combine(const FissionYields& yields, const std::vector<double>& weights){
                std::map<int, double> combined;
                for(std::size_t k = 0; k < weights.size(); ++k){
                    if(weights[k] == 0.0){
                        continue;
                    }
                    for(const auto& [product, y]: yields.yields[k]){
                        combined[product] += weights[k]*y;
                    }
                }
                return combined;
            }
        }

        std::map<int, double> yields_at_energy(const FissionYields& yields, double energy){
            return _combine(yields, energy_weights(yields.energies, energy));
        }

        std::map<int, double> weighted_yields(const FissionYields& yields,
                                              const std::vector<double>& bounds,
                                              const std::vector<double>& group_rates){
            if(bounds.size() != group_rates.size() + 1){
                throw std::invalid_argument("There must be one more group bound than group rates");
            }
            std::vector<double> weights(yields.energies.size(), 0.0);
            double total = 0.0;
            for(std::size_t g = 0; g < group_rates.size(); ++g){
                if(group_rates[g] <= 0.0){
                    continue;
                }
                // the middle of the group in lethargy, or in energy for a group from 0
                const double low = bounds[g];
                const double high = bounds[g+1];
                const double mid = (low > 0.0) ? std::sqrt(low*high) : 0.5*high;
                const auto w = energy_weights(yields.energies, mid);
                for(std::size_t k = 0; k < w.size(); ++k){
                    weights[k] += group_rates[g]*w[k];
                }
                total += group_rates[g];
            }
            if(total <= 0.0){
                return yields_at_energy(yields, yields.energies.front());
            }
            for(auto& w: weights){
                w /= total;
            }
            return _combine(yields, weights);
        }
    }
}
