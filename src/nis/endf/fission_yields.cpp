#include <cmath>
#include <sstream>
#include <string>
#include <stdexcept>

#include "nis/endf/fission_yields.hpp"

namespace nis {
    namespace endf {

        namespace {

            // The isomeric state of the material: LISO, the 4th field (L2) of the 2nd line of MF=1 MT=451.
            // We read it from the text, as some yield files (JEFF-3.1.1) have a directory ENDFtk cannot parse.
            int _liso(const njoy::ENDFtk::tree::Material& material){
                if(!material.hasMF(1) || !material.file(1).hasMT(451)){
                    return 0;
                }
                std::istringstream lines(material.file(1).MT(451).content());
                std::string line;
                if(!std::getline(lines, line) || !std::getline(lines, line) || (line.size() < 44)){
                    return 0;
                }
                try{
                    return std::stoi(line.substr(33, 11));
                }
                catch(const std::exception&){
                    return 0;
                }
            }

            void _read_material(const njoy::ENDFtk::tree::Material& material, FissionYieldData& data){
                if(!material.hasMF(FISSION_YIELD_FILE) || !material.file(FISSION_YIELD_FILE).hasMT(INDEPENDENT_YIELDS_MT)){
                    return;
                }
                const auto section = material.file(FISSION_YIELD_FILE).MT(INDEPENDENT_YIELDS_MT)
                                         .parse<FISSION_YIELD_FILE, INDEPENDENT_YIELDS_MT>();

                FissionYields yields;
                yields.parent_zai = static_cast<int>(std::round(section.ZA()))*10 + _liso(material);
                for(const auto& entry: section.yields()){
                    std::map<int, double> products;
                    const auto zafp = entry.ZAFP();
                    const auto fps = entry.FPS();
                    const auto y = entry.Y();
                    for(std::size_t i = 0; i < zafp.size(); ++i){
                        const int product = static_cast<int>(std::round(zafp[i]))*10
                                            + static_cast<int>(std::round(fps[i]));
                        // the same product can be given more than once, add them
                        products[product] += y[i][0];
                    }
                    yields.energies.push_back(entry.E());
                    yields.yields.push_back(std::move(products));
                }
                // keep the energies in order
                for(std::size_t i = 1; i < yields.energies.size(); ++i){
                    if(yields.energies[i] < yields.energies[i-1]){
                        throw std::runtime_error("The fission yield energies of ZAI " + std::to_string(yields.parent_zai)
                                                 + " are not in order");
                    }
                }
                if(!yields.energies.empty()){
                    data[yields.parent_zai] = std::move(yields);
                }
            }
        }

        FissionYieldData read_fission_yields_from_tape(const njoy::ENDFtk::tree::Tape& tape){
            FissionYieldData data;
            for(auto mat_number: tape.materialNumbers()){
                for(const auto& material: tape.material(mat_number)){
                    _read_material(material, data);
                }
            }
            return data;
        }

        FissionYieldData read_fission_yields_from_file(const std::string& filename){
            return read_fission_yields_from_tape(njoy::ENDFtk::tree::fromFile(filename));
        }
    }
}
