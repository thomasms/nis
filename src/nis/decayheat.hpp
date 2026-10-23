#ifndef DECAY_HEAT_HPP
#define DECAY_HEAT_HPP

#include <cassert>
#include <iostream>

#include <Eigen/Core>

#include "nis/nuclides.hpp"
#include "nis/physics.hpp"

namespace nis {
    namespace core {

        // All vectors are indexed the same as the inventory (data::ZAI_MAP)
        // and are zero for stable nuclides and the sink.
        template<typename T=double>
        struct DecayEnergyData {
            DecayEnergyData() :
                decay_constants(Eigen::VectorX<T>::Zero(data::ZAI_SIZE)),
                light_particle(Eigen::VectorX<T>::Zero(data::ZAI_SIZE)),
                electromagnetic(Eigen::VectorX<T>::Zero(data::ZAI_SIZE)),
                heavy_particle(Eigen::VectorX<T>::Zero(data::ZAI_SIZE)) {}

            // 1/s
            Eigen::VectorX<T> decay_constants;
            // average energy per decay in eV
            // beta, electrons and positrons, ENDF does not include neutrinos
            Eigen::VectorX<T> light_particle;
            // gamma and x-rays
            Eigen::VectorX<T> electromagnetic;
            // alpha, protons, neutrons, fission fragments and recoil
            Eigen::VectorX<T> heavy_particle;
        };

        // decay heat in W
        template<typename T=double>
        struct DecayHeat {
            T light_particle = 0;
            T electromagnetic = 0;
            T heavy_particle = 0;

            T total() const {
                return light_particle + electromagnetic + heavy_particle;
            }
        };

        // activity per nuclide in Bq
        template<typename T=double>
        Eigen::VectorX<T> compute_activities(const DecayEnergyData<T>& energy_data,
                                             const Eigen::VectorX<T>& inventory){
            assert(inventory.size() == energy_data.decay_constants.size());
            return energy_data.decay_constants.cwiseProduct(inventory);
        }

        // total activity in Bq
        template<typename T=double>
        T compute_activity(const DecayEnergyData<T>& energy_data,
                           const Eigen::VectorX<T>& inventory){
            return compute_activities<T>(energy_data, inventory).sum();
        }

        // total decay heat in W
        template<typename T=double>
        DecayHeat<T> compute_decay_heat(const DecayEnergyData<T>& energy_data,
                                        const Eigen::VectorX<T>& inventory){
            const Eigen::VectorX<T> activities = compute_activities<T>(energy_data, inventory);
            DecayHeat<T> heat;
            heat.light_particle = activities.dot(energy_data.light_particle)*physics::EV_TO_J;
            heat.electromagnetic = activities.dot(energy_data.electromagnetic)*physics::EV_TO_J;
            heat.heavy_particle = activities.dot(energy_data.heavy_particle)*physics::EV_TO_J;
            return heat;
        }

        template<typename T=double>
        inline std::ostream& operator<<(std::ostream& os, const DecayHeat<T>& heat) {
            os << "\ndecay_heat\n";
            os << "    total " << heat.total() << "\n";
            os << "    light_particle " << heat.light_particle << "\n";
            os << "    electromagnetic " << heat.electromagnetic << "\n";
            os << "    heavy_particle " << heat.heavy_particle << "\n";
            return os;
        }

    }
}

#endif // DECAY_HEAT_HPP
