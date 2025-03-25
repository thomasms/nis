#ifndef CORE_MODEL_HPP
#define CORE_MODEL_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

namespace nis{
    namespace core{

        template<typename T=double>
        struct BasicModel {
            T dt;
            std::map<std::string, T> nuclides;
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