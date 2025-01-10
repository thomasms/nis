#ifndef IO_HELPERS_HPP
#define IO_HELPERS_HPP

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <Eigen/Core>
#include "nis/io/input.hpp"
#include "nis/io/output.hpp"
#include "nis/nuclides.hpp"
#include "nis/util.hpp"
#include "nis/zai.hpp"

using namespace nis::io::input;
using namespace nis::io::output;

namespace nis{
    namespace io {
        namespace helpers{
            
            template<typename T>
            Eigen::VectorX<T> vector_from_inventory(const std::map<std::string, T>& inventory){
                Eigen::VectorX<T> vec = Eigen::VectorX<T>::Zero(data::ZAI_SIZE);

                for (auto it = inventory.begin(); it != inventory.end(); ++it) {
                    const int nuclide_zai = zai::get_from_name(it->first);
                    if(data::ZAI_MAP.count(nuclide_zai) == 1){
                        const auto zai_index = data::ZAI_MAP.at(nuclide_zai);
                        vec(zai_index) = it->second;
                    }
                }

                return vec;
            }

            // vector must be of size data::ZAI_SIZE
            template<typename T>
            std::map<std::string, T> inventory_from_vector(const Eigen::VectorX<T>& vec, const T& threshold=T()){
                assert(data::ZAI_SIZE == vec.size());
                std::map<std::string, T> inv;
                // this might be slow, better to do this once in one place
                // or at compile time - how?
                const auto reversed = util::invert_map<int, int>(data::ZAI_MAP);

                const int vec_size = vec.size();
                for(int i=0; i<vec_size;++i){
                    int zai = (reversed.count(i) == 0) ? -1 : reversed.at(i);
                    std::string nuclide_name = zai::get_name_from(zai);
                    const auto val = vec(i);
                    if(val > threshold){
                        inv.insert(std::make_pair(nuclide_name, val));
                    }
                }

                return inv;
            }

            template<typename T, template <typename> class Container>
            Container<T> from_file(std::string filename){
                Container<T> data;

                std::ifstream file(filename);
                if (!file) {
                    std::cerr << "Error opening file: " << filename << std::endl;
                    return data;
                }

                file >> data;

                return data;
            }

            template<typename T, template <typename> class Container>
            void to_file(const Container<T>& data, std::string filename){
                std::ofstream file(filename);

                if (!file) {
                    std::cerr << "Error opening file: " << filename << std::endl;
                    return;
                }

                file << data;
            }
        }
    }
}

#endif