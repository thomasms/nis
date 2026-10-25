#ifndef CORE_MODEL_HPP
#define CORE_MODEL_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <memory>
#include <optional>
#include <vector>

#include "nis/spectrum.hpp"

namespace nis{
    namespace core{

        // one irradiation or cooling step of the input, a "dt <seconds> [flux]" line
        template<typename T=double>
        struct TimeStep {
            T dt;
            // total flux in n/m^2/s, the spectrum is normalised to it
            // 0 is a cooling step (decay only), no value uses the spectrum as it is given
            std::optional<T> flux;
        };

        template<typename T=double>
        struct BasicModel {
            // the step length for compute_inventory, the first step of the input
            T dt;
            std::map<std::string, T> nuclides;
            // we use a pointer here so we can
            // leave it null for decay only
            // otherwise we do a collapse for no reason
            std::shared_ptr<Spectrum<T>> spectrum = nullptr;
            // all the steps of the input, in order
            std::vector<TimeStep<T>> steps;
        };

        template<typename T=double>
        struct ModelResult {
            std::map<std::string, T> nuclides;
        };

        template<typename T=double>
        ModelResult<T> filter_nuclides(const ModelResult<T>& unfiltered, T threshold){
            std::map<std::string, T> filtered_nuclides;
            std::copy_if(unfiltered.nuclides.begin(),
                         unfiltered.nuclides.end(),
                         std::inserter(filtered_nuclides, filtered_nuclides.begin()),
                         [threshold](const auto& nuclide){
                            return nuclide.second >= threshold;
                        }
            );
            return core::ModelResult<double>{filtered_nuclides};
        }

    }
}

#endif