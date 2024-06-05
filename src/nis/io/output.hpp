#ifndef IO_OUTPUT_HPP
#define IO_OUTPUT_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

namespace nis{
    namespace io {
        namespace output{

            template<typename T>
            class OutputData;

            template<typename T=double>
            void serialize(std::ostream& os, const OutputData<T>& data) {
                os << "\nnuclides\n";
                for (const auto& nuclide : data.nuclides) {
                    os << "    " << nuclide.first << " " << nuclide.second << "\n";
                }
            }

            template<typename T=double>
            OutputData<T> deserialize(std::istream& is) {
                OutputData<T> data;

                std::string line;
                std::string key, value;
                std::istringstream iss;
                
                while (std::getline(is, line)) {
                    if (line.empty() || line[0] == '#') {
                        continue;
                    }else if (line == "nuclides") {
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
            struct OutputData {
                std::map<std::string, T> nuclides;
                
                inline friend std::ostream& operator<<(std::ostream& os, const OutputData<T>& data) {
                    serialize<T>(os, data);
                    return os;
                }

                inline friend std::istream& operator>>(std::istream& is, OutputData<T>& data) {
                    data = deserialize<T>(is);
                    return is;
                }
            };

        }
    }
}

#endif