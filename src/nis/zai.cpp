#include <iostream>
#include <sstream>

#include "nis/nuclides.hpp"
#include "nis/util.hpp"

#include "nis/zai.hpp"

namespace nis {
    namespace zai {

        int get_from_name(std::string name){
            // need better error checking for strings not well formatted.
            int z = 0;
            int a = 0;
            int i = 0;

            // Remove numeric values from the string and lowercase it
            std::string trimmed_string = util::string_rtrim(name);
            std::string name_lower = util::string_lower(trimmed_string);
            // if no metastable then this will be an integer
            unsigned char i_char = name_lower.back();
            name_lower.pop_back();
            name_lower.erase(std::remove_if(name_lower.begin(), name_lower.end(), [](char c) {
                return std::isdigit(c);
            }), name_lower.end());

            // find Z
            for (const auto& pair: data::Z_ELEMENT_MAP) {
                if(util::string_lower(pair.second) == name_lower){
                    z = pair.first;
                    break;
                }
            }

            // extract A
            trimmed_string.erase(std::remove_if(trimmed_string.begin(), trimmed_string.end(), [](char c) {
                return !std::isdigit(c);
            }), trimmed_string.end());
            if(!trimmed_string.empty()){
                a = std::stoi(trimmed_string);
            }

            // find I
            for (const auto& pair: META_STABLE_MAP) {
                if(std::tolower(pair.second) == i_char){
                    i = pair.first;
                    break;
                }
            }

            return make(z, a, i);
        }

        std::string get_name_from(int zai){
            triple z_a_i = decompose(zai);
            const int z = std::get<0>(z_a_i);
            const int a = std::get<1>(z_a_i);
            const int i = std::get<2>(z_a_i);

            // no matching Z then return null string
            if((data::Z_ELEMENT_MAP.count(z) == 0) || (META_STABLE_MAP.count(i) == 0)){
                return "";
            }

            std::string element_name = data::Z_ELEMENT_MAP.at(z);
            std::stringstream ss;
            ss << element_name << std::to_string(a) << META_STABLE_MAP.at(i);
            std::string name = util::string_rtrim(ss.str());

            return name.empty() ? "SINK" : name;
        }

    }
}
