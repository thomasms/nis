#ifndef IO_INPUT_HPP
#define IO_INPUT_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

namespace nis{
    namespace io {
        namespace input{

            template<typename T>
            class InputData;

            template<typename T=double>
            void serialize(std::ostream& os, const InputData<T>& data) {
                os << "dt " << data.dt << "\n\nnuclides\n";
                for (const auto& nuclide : data.nuclides) {
                    os << "    " << nuclide.first << " " << nuclide.second << "\n";
                }
            }

            template<typename T=double>
            InputData<T> deserialize(std::istream& is) {
                InputData<T> data;

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
                    }
                }

                return data;
            }

            template<typename T=double>
            struct InputData {
                T dt;
                std::map<std::string, T> nuclides;
                
                inline friend std::ostream& operator<<(std::ostream& os, const InputData<T>& data) {
                    serialize<T>(os, data);
                    return os;
                }

                inline friend std::istream& operator>>(std::istream& is, InputData<T>& data) {
                    data = deserialize<T>(is);
                    return is;
                }
            };

        }
    }
}

#endif