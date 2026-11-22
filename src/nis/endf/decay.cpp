#include <cmath>
#include <limits>

#include "nis/endf/decay.hpp"


namespace nis {
    namespace endf {

        namespace {

            constexpr int file_number = 8;
            constexpr int mt_number = 457;

            using DecaySection = njoy::ENDFtk::section::Type<file_number, mt_number>;

            int _get_zai_from_section(const DecaySection& section){
                const int za = section.ZA();
                // ENDFtk gives us the isomeric state as an unsigned int
                const int i = static_cast<int>(section.LISO());
                return za*10 + i;
            }

            data::DecayNuclide _read_nuclide(const DecaySection& section){
                data::DecayNuclide nuclide;
                const double half_life = section.halfLife()[0];
                // a stable nuclide (or one with no half-life) keeps the infinite half-life, with no modes
                if(section.isStable() || (half_life <= 0.)){
                    return nuclide;
                }
                nuclide.half_life = half_life;
                for(const auto& decay_mode: section.decayModes().decayModes()){
                    data::DecayMode mode;
                    mode.rtyp = decay_mode.RTYP();
                    mode.branching_ratio = decay_mode.branchingRatio()[0];
                    mode.final_state = static_cast<int>(std::floor(decay_mode.finalIsomericState()));
                    nuclide.modes.push_back(mode);
                }
                const auto& energies = section.averageDecayEnergies();
                nuclide.light_particle = energies.ELP()[0];
                nuclide.electromagnetic = energies.EEM()[0];
                nuclide.heavy_particle = energies.EHP()[0];
                return nuclide;
            }

            void _read_material(const njoy::ENDFtk::tree::Material& material, data::DecayData& decay_data){
                if(!material.hasMF(file_number)){
                    return;
                }
                njoy::ENDFtk::tree::File file = material.file(file_number);
                if(!file.hasMT(mt_number)){
                    return;
                }
                const auto section = file.MT(mt_number).parse<file_number, mt_number>();
                // if the same nuclide comes again, the last one is used
                decay_data[_get_zai_from_section(section)] = _read_nuclide(section);
            }
        }

        data::DecayData read_decay_data_from_tape(const njoy::ENDFtk::tree::Tape& tape){
            data::DecayData decay_data;
            // loop over all materials in the file. One MAT number can come more than once,
            // i.e. at different temperatures, so loop over each instance too.
            for(auto mat_number: tape.materialNumbers()){
                for(auto material: tape.material(mat_number)){
                    _read_material(material, decay_data);
                }
            }
            return decay_data;
        }

        data::DecayData read_decay_data_from_file(const std::string& filename){
            return read_decay_data_from_tape(njoy::ENDFtk::tree::fromFile(filename));
        }

        Eigen::SparseMatrix<double> construct_decay_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape,
                                                                     const FissionYieldData* sf_yields){
            return data::construct_decay_matrix(read_decay_data_from_tape(tape), sf_yields);
        }

        Eigen::SparseMatrix<double> construct_decay_matrix_from_file(std::string filename,
                                                                     const FissionYieldData* sf_yields){
            return data::construct_decay_matrix(read_decay_data_from_file(filename), sf_yields);
        }

        core::DecayEnergyData<double> read_decay_energy_data_from_tape(const njoy::ENDFtk::tree::Tape& tape){
            return data::construct_decay_energy_data(read_decay_data_from_tape(tape));
        }

        core::DecayEnergyData<double> read_decay_energy_data_from_file(std::string filename){
            return data::construct_decay_energy_data(read_decay_data_from_file(filename));
        }

    }
}
