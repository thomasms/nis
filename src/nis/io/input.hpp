#ifndef IO_INPUT_HPP
#define IO_INPUT_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <memory>

#include "nis/spectrum.hpp"
#include "nis/model.hpp"
#include "nis/util.hpp"

namespace nis{
    namespace io {
        namespace input {

            template<typename T=double>
            static void Serialize(std::ostream& os, const core::BasicModel<T>& data) {
                os << "dt " << data.dt << "\n\nnuclides\n";
                for (const auto& nuclide : data.nuclides) {
                    os << "    " << nuclide.first << " " << nuclide.second << "\n";
                }
                if(data.spectrum != nullptr){
                    os << "\nspectrum\n";
                    os << "    ";
                    for(auto v: data.spectrum->x()){
                        os << v << " ";
                    }
                    os << "\n" << "    ";
                    for(auto v: data.spectrum->y()){
                        os << v << " ";
                    }
                    os << "\n";
                }
            }

            template<typename T=double>
            static core::BasicModel<T> Deserialize(std::istream& is) {
                core::BasicModel<T> data;

                std::string line;
                std::string key, value;
                std::istringstream iss;
                
                while (std::getline(is, line)) {
                    if (line.empty() || line[0] == '#') {
                        continue;
                    }
                    if (line.substr(0, 2) == "dt") {
                        iss.str(line.substr(3));
                        iss >> data.dt;
                        iss.clear();
                    } else if (line == "nuclides") {
                        while (std::getline(is, line)) {
                            if (line.empty() || line[0] == '#') {
                                break;
                            }
                            iss.str(line.substr(4));
                            iss >> key >> value;
                            data.nuclides[key] = std::stod(value);
                            iss.clear();
                        }
                    } else if (line == "spectrum") {
                        std::vector<double> energies;
                        std::vector<double> values;
                        int line_count = 0;
                        while (std::getline(is, line)) {
                            if (line.empty() || line[0] == '#') {
                                break;
                            }
                            line_count++;
                            iss.str(line.substr(4));
                            std::vector<std::string> contents;
                            util::string_split<' '>(iss, contents);
                            if(line_count == 1){
                                for(auto s: contents){
                                    energies.emplace_back(std::stod(s));
                                }
                            }
                            if(line_count == 2){
                                for(auto s: contents){
                                    values.emplace_back(std::stod(s));
                                }
                                break;
                            }
                            iss.clear();
                        }
                        data.spectrum = std::make_shared<core::Spectrum<double>>(energies, values);
                    }
                }

                return data;
            }

            template<typename T=double>
            inline std::ostream& operator<<(std::ostream& os, const core::BasicModel<T>& data) {
                Serialize<T>(os, data);
                return os;
            }

            template<typename T=double>
            inline std::istream& operator>>(std::istream& is, core::BasicModel<T>& data) {
                data = Deserialize<T>(is);
                return is;
            }

            // Explicit instantiation for common types
            template std::ostream& operator<<<double>(std::ostream& os, const core::BasicModel<double>& data);
            template std::istream& operator>><double>(std::istream& is, core::BasicModel<double>& data);
        }
    }
}

#endif