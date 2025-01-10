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
    }
}

#endif