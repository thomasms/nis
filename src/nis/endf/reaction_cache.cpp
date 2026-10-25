#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <stdexcept>

#include "nis/endf/reaction_cache.hpp"

namespace nis {
    namespace endf {

        // the header: 8 characters, the version, then a number to check the byte order
        constexpr char CACHE_MAGIC[8] = {'N', 'I', 'S', 'R', 'E', 'A', 'C', '\0'};
        constexpr std::uint32_t CACHE_VERSION = 1;
        constexpr std::uint32_t CACHE_BYTE_ORDER = 0x01020304;

        template<typename T>
        void _write(std::ofstream& out, const T& value){
            out.write(reinterpret_cast<const char*>(&value), sizeof(T));
        }

        void _write(std::ofstream& out, const double* values, std::uint64_t size){
            out.write(reinterpret_cast<const char*>(values), static_cast<std::streamsize>(size*sizeof(double)));
        }

        template<typename T>
        T _read(std::ifstream& in){
            T value;
            in.read(reinterpret_cast<char*>(&value), sizeof(T));
            if(!in){
                throw std::runtime_error("The reaction cache ends too soon");
            }
            return value;
        }

        void _read(std::ifstream& in, double* values, std::uint64_t size){
            in.read(reinterpret_cast<char*>(values), static_cast<std::streamsize>(size*sizeof(double)));
            if(!in){
                throw std::runtime_error("The reaction cache ends too soon");
            }
        }

        void write_reaction_cache(const ReactionData& reaction_data, const std::string& filename){
            std::ofstream out(filename, std::ios::binary);
            if(!out){
                throw std::runtime_error("Could not open file '" + filename + "' to write");
            }

            // each energy grid once, in the order we first see it
            std::map<std::vector<double>, std::uint32_t> grid_index;
            std::vector<const std::vector<double>*> grids;
            for(const auto& xs: reaction_data){
                if(grid_index.emplace(xs.energies, static_cast<std::uint32_t>(grids.size())).second){
                    grids.push_back(&xs.energies);
                }
            }

            out.write(CACHE_MAGIC, sizeof(CACHE_MAGIC));
            _write(out, CACHE_VERSION);
            _write(out, CACHE_BYTE_ORDER);

            _write(out, static_cast<std::uint64_t>(grids.size()));
            for(const auto* grid: grids){
                _write(out, static_cast<std::uint64_t>(grid->size()));
                _write(out, grid->data(), grid->size());
            }

            _write(out, static_cast<std::uint64_t>(reaction_data.size()));
            for(const auto& xs: reaction_data){
                if(xs.values.size() != xs.energies.size()){
                    throw std::invalid_argument("A cross section must have the same number of energies and values");
                }
                // only the values from the first to the last non-zero, the rest are 0
                std::uint64_t first = 0;
                std::uint64_t last = xs.values.size();
                while((first < last) && (xs.values[first] == 0.0)){
                    ++first;
                }
                while((last > first) && (xs.values[last-1] == 0.0)){
                    --last;
                }

                _write(out, static_cast<std::int32_t>(xs.parent_zai));
                _write(out, static_cast<std::int32_t>(xs.mt_number));
                _write(out, static_cast<std::int32_t>(xs.final_state));
                _write(out, static_cast<std::int32_t>(xs.product_za));
                _write(out, grid_index.at(xs.energies));
                _write(out, first);
                _write(out, last - first);
                _write(out, xs.values.data() + first, last - first);
            }

            if(!out){
                throw std::runtime_error("Could not write the reaction cache '" + filename + "'");
            }
        }

        bool is_reaction_cache(const std::string& filename){
            std::ifstream in(filename, std::ios::binary);
            char magic[sizeof(CACHE_MAGIC)] = {};
            in.read(magic, sizeof(magic));
            return in && (std::memcmp(magic, CACHE_MAGIC, sizeof(CACHE_MAGIC)) == 0);
        }

        ReactionData read_reaction_cache(const std::string& filename){
            std::ifstream in(filename, std::ios::binary);
            if(!in){
                throw std::runtime_error("Could not open file '" + filename + "'");
            }

            char magic[sizeof(CACHE_MAGIC)] = {};
            in.read(magic, sizeof(magic));
            if(!in || (std::memcmp(magic, CACHE_MAGIC, sizeof(CACHE_MAGIC)) != 0)){
                throw std::runtime_error("'" + filename + "' is not a reaction cache");
            }
            if(_read<std::uint32_t>(in) != CACHE_VERSION){
                throw std::runtime_error("'" + filename + "' is a reaction cache of another version, make it again");
            }
            if(_read<std::uint32_t>(in) != CACHE_BYTE_ORDER){
                throw std::runtime_error("'" + filename + "' was made on a machine with another byte order, make it again");
            }

            std::vector<std::vector<double>> grids(_read<std::uint64_t>(in));
            for(auto& grid: grids){
                grid.resize(_read<std::uint64_t>(in));
                _read(in, grid.data(), grid.size());
            }

            ReactionData reaction_data(_read<std::uint64_t>(in));
            for(auto& xs: reaction_data){
                xs.parent_zai = _read<std::int32_t>(in);
                xs.mt_number = _read<std::int32_t>(in);
                xs.final_state = _read<std::int32_t>(in);
                xs.product_za = _read<std::int32_t>(in);
                const auto grid = _read<std::uint32_t>(in);
                const auto first = _read<std::uint64_t>(in);
                const auto count = _read<std::uint64_t>(in);
                if((grid >= grids.size()) || (first + count > grids[grid].size())){
                    throw std::runtime_error("The reaction cache '" + filename + "' is not valid");
                }
                xs.energies = grids[grid];
                xs.values.assign(xs.energies.size(), 0.0);
                _read(in, xs.values.data() + first, count);
            }
            return reaction_data;
        }

        ReactionData read_reaction_data(const std::string& filename){
            if(is_reaction_cache(filename)){
                return read_reaction_cache(filename);
            }
            return read_reaction_data_from_file(filename);
        }
    }
}
