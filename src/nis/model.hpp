#ifndef CORE_MODEL_HPP
#define CORE_MODEL_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <memory>
#include <vector>

#include "nis/spectrum.hpp"

namespace nis{
    namespace core{

        template<typename T=double>
        struct BasicModel {
            T dt;
            std::map<std::string, T> nuclides;
            // we use a pointer here so we can
            // leave it null for decay only
            // otherwise we do a collapse for no reason
            std::shared_ptr<Spectrum<T>> spectrum = nullptr;
            // decay only steps after dt, each step length in seconds
            std::vector<T> cooling;
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