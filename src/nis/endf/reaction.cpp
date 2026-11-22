#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include "nis/physics.hpp"
#include "nis/nuclides.hpp"
#include "nis/util.hpp"
#include "nis/zai.hpp"

#include "nis/endf/core.hpp"
#include "nis/endf/reaction.hpp"
#include "nis/maths.hpp"
#include "nis/spectrum.hpp"


namespace nis {
    namespace endf {

        constexpr int file_number = 3;
        // production cross sections for each final state of the residual
        constexpr int isomer_file_number = 10;
        // 1= y is constant in x (constant, histogram)
        constexpr int histogram_interpolation = 1;

        std::map<std::pair<int, int>, int> isomer_states_from_levels(const std::map<int, std::set<int>>& levels){
            std::map<std::pair<int, int>, int> states;
            for(const auto& [za, level_numbers]: levels){
                // the set is sorted, so the levels come in order of energy
                int isomer = 0;
                for(auto lfs: level_numbers){
                    states[{za, lfs}] = (lfs == 0) ? 0 : ++isomer;
                }
            }
            return states;
        }

        template<typename Section>
        void _check_histogram(const Section& section, int parent_zai, int file, int mt_number){
            // we assume everything is already preprocessed into INT=1 (group-wise)
            for(auto interpolation_scheme: section.interpolants()){
                if(interpolation_scheme != histogram_interpolation){
                    throw std::runtime_error("Only INT=1 (histogram) cross sections are supported, found INT="
                                             + std::to_string(interpolation_scheme)
                                             + " for ZAI=" + std::to_string(parent_zai)
                                             + " MF=" + std::to_string(file)
                                             + " MT=" + std::to_string(mt_number));
                }
            }
        }

        using IsomerStates = std::map<std::pair<int, int>, int>;
        using MF10Sections = std::map<int, njoy::ENDFtk::section::Type<isomer_file_number>>;

        // all the MF=10 sections of the material, each parsed once
        MF10Sections _mf10_sections(const njoy::ENDFtk::tree::Material& material){
            MF10Sections sections;
            if(material.hasMF(isomer_file_number)){
                njoy::ENDFtk::tree::File file = material.file(isomer_file_number);
                for(auto mt_number: file.sectionNumbers()){
                    sections.emplace(mt_number, file.MT(mt_number).parse<isomer_file_number>());
                }
            }
            return sections;
        }

        // all the levels of each product in MF=10, to give the isomeric states
        IsomerStates _isomer_states(const MF10Sections& sections){
            std::map<int, std::set<int>> levels;
            for(const auto& [mt_number, section]: sections){
                // in MT=5 the LFS is already the isomer number (as in JEFF-4 and TENDL), so we leave it out
                if(mt_number == ANYTHING_MT){
                    continue;
                }
                for(const auto& product: section.reactionProducts()){
                    levels[product.IZAP()].insert(product.LFS());
                }
            }
            return isomer_states_from_levels(levels);
        }

        // calls func(product) for each product of the MF=10 section, if there is one
        template<typename Func>
        void _for_each_mf10_product(const MF10Sections& sections, int mt_number, Func func){
            const auto it = sections.find(mt_number);
            if(it == sections.end()){
                return;
            }
            for(const auto& product: it->second.reactionProducts()){
                func(product);
            }
        }

        void _read_isomer_data(const MF10Sections& sections, int parent_zai,
                               int mt_number, const IsomerStates& states, ReactionData& reaction_data){
            // the residual we expect for this MT, MF=10 can also list other products
            const int residual_zai = std::get<0>(compute_reaction_products(parent_zai, mt_number));
            if(residual_zai < 0){
                return;
            }
            const int residual_za = residual_zai/10;

            _for_each_mf10_product(sections, mt_number, [&](const auto& product){
                if(product.IZAP() != residual_za){
                    return;
                }
                _check_histogram(product, parent_zai, isomer_file_number, mt_number);
                reaction_data.push_back(CrossSection{
                    parent_zai,
                    mt_number,
                    {product.energies().begin(), product.energies().end()},
                    {product.crossSections().begin(), product.crossSections().end()},
                    states.at({residual_za, product.LFS()})
                });
            });
        }

        // MT=4 and MT=5 only have data we can use in MF=10
        void _read_mf10_only_data(const MF10Sections& sections, int parent_zai,
                                  const IsomerStates& states, ReactionData& reaction_data){
            const int parent_za = parent_zai/10;

            // (n,n') to the isomers of the target, the ground state does not change the nuclide
            // TODO: (n,n') on a target isomer can also go down to the ground state
            if(zai::get_i(parent_zai) == 0){
                _for_each_mf10_product(sections, INELASTIC_MT, [&](const auto& product){
                    const int state = states.at({product.IZAP(), product.LFS()});
                    if((product.IZAP() != parent_za) || (state == 0)){
                        return;
                    }
                    _check_histogram(product, parent_zai, isomer_file_number, INELASTIC_MT);
                    reaction_data.push_back(CrossSection{
                        parent_zai,
                        INELASTIC_MT,
                        {product.energies().begin(), product.energies().end()},
                        {product.crossSections().begin(), product.crossSections().end()},
                        state
                    });
                });
            }

            // (n,anything), each product with its own cross section
            _for_each_mf10_product(sections, ANYTHING_MT, [&](const auto& product){
                // IZAP=0 is photons, they are not tracked
                if(product.IZAP() <= 0){
                    return;
                }
                _check_histogram(product, parent_zai, isomer_file_number, ANYTHING_MT);
                reaction_data.push_back(CrossSection{
                    parent_zai,
                    ANYTHING_MT,
                    {product.energies().begin(), product.energies().end()},
                    {product.crossSections().begin(), product.crossSections().end()},
                    static_cast<int>(product.LFS()),
                    product.IZAP()
                });
            });
        }

        void _read_reaction_data(const njoy::ENDFtk::tree::Material& material,
                                 const std::map<int, int>& zai_map, ReactionData& reaction_data){

            const int parent_zai = get_parent_zai(material);
            if(zai_map.count(parent_zai) == 0){
                return;
            }

            // TODO: MF=10 can give isomer production when MF=3 has no data for it
            // TODO: MF=9 (multiplicities) is not used, JEFF-4 only has MF=10
            const MF10Sections sections = _mf10_sections(material);
            const IsomerStates states = _isomer_states(sections);
            njoy::ENDFtk::tree::File file = material.file(file_number);
            for(auto mt_number: file.sectionNumbers()){
                // we skip sums, elastic, inelastic levels and anything we can't map
                if(!is_reaction_mt(mt_number)){
                    continue;
                }

                njoy::ENDFtk::tree::Section section = file.MT(mt_number);
                auto parsed_section = section.parse<file_number>();

                _check_histogram(parsed_section, parent_zai, file_number, mt_number);

                reaction_data.push_back(CrossSection{
                    parent_zai,
                    mt_number,
                    {parsed_section.energies().begin(), parsed_section.energies().end()},
                    {parsed_section.crossSections().begin(), parsed_section.crossSections().end()}
                });

                _read_isomer_data(sections, parent_zai, mt_number, states, reaction_data);
            }

            _read_mf10_only_data(sections, parent_zai, states, reaction_data);
        }

        // adds the data of one material, if we have not seen its nuclide before
        void _read_material(const njoy::ENDFtk::tree::Material& material,
                            std::set<int>& done_zais, ReactionData& reaction_data){
            if(!material.hasMF(file_number)){
                return;
            }
            const int parent_zai = get_parent_zai(material);
            if(done_zais.count(parent_zai) == 0){
                _read_reaction_data(material, data::ZAI_MAP, reaction_data);
                done_zais.insert(parent_zai);
            }
        }

        ReactionData read_reaction_data_from_tape(const njoy::ENDFtk::tree::Tape& tape){
            ReactionData reaction_data;

            // the same nuclide can appear in more than one MAT, we only want it once
            std::set<int> done_zais;

            // loop over all materials in the file
            for(auto mat_number: tape.materialNumbers()){
                // this can return multiple instances of the same MAT number
                // due to different temperatures.
                // To not double count we only use the first temperature.
                // TODO: let the user select the temperature
                auto materials = tape.material(mat_number);
                for(auto material: materials){
                    if(material.hasMF(file_number)){
                        _read_material(material, done_zais, reaction_data);
                        break;
                    }
                }
            }

            return reaction_data;
        }

        // the MAT number is in columns 67-70 of each ENDF line
        int _mat_number(const std::string& line){
            if(line.size() < 70){
                throw std::runtime_error("ENDF line is too short: '" + line + "'");
            }
            return std::stoi(line.substr(66, 4));
        }

        // the MF and MT numbers are in columns 71-72 and 73-75
        std::pair<int, int> _mf_mt(const std::string& line){
            if(line.size() < 75){
                throw std::runtime_error("ENDF line is too short: '" + line + "'");
            }
            return {std::stoi(line.substr(70, 2)), std::stoi(line.substr(72, 3))};
        }

        // the sections we read: the ZA in MF=1 MT=451, the reactions in MF=3 and all of MF=10
        bool _is_used_section(int mf_number, int mt_number){
            return ((mf_number == 1) && (mt_number == 451))
                || ((mf_number == file_number) && is_reaction_mt(mt_number))
                || (mf_number == isomer_file_number);
        }

        ReactionData read_reaction_data_from_file(std::string filename){
            // A reaction tape can be ~1 GB, and parsing it all at once needs ~2x that.
            // So we read one material at a time and parse it as a small tape of its own.
            // We keep only the sections we use, because making the tape index of the
            // others (e.g. MF=33 covariances) takes most of the time.
            std::ifstream in(filename);
            if(!in){
                throw std::runtime_error("Could not open file '" + filename + "'");
            }

            // the TEND record that ends a tape
            const std::string tend = std::string(66, ' ') + "  -1 0  0    0\n";

            ReactionData reaction_data;
            std::set<int> done_zais;

            std::string line;
            // the first line is the tape id (TPID)
            std::getline(in, line);
            const std::string tpid = line + "\n";

            std::string material_lines;
            bool section_used = false;
            bool file_used = false;
            while(std::getline(in, line)){
                if(!line.empty() && line.back() == '\r'){
                    line.pop_back();
                }
                if(line.empty()){
                    continue;
                }
                const int mat_number = _mat_number(line);
                if(mat_number == -1){
                    break;
                }
                // keep the SEND and FEND records of the sections and files we keep
                const auto [mf_number, mt_number] = _mf_mt(line);
                bool keep = true;
                if(mat_number != 0){
                    if(mf_number == 0){
                        keep = file_used;
                        file_used = false;
                    }
                    else if(mt_number == 0){
                        keep = section_used;
                    }
                    else{
                        section_used = _is_used_section(mf_number, mt_number);
                        file_used = file_used || section_used;
                        keep = section_used;
                    }
                }
                if(keep){
                    material_lines += line;
                    material_lines += '\n';
                }
                // the MEND record (MAT=0) ends a material
                if(mat_number == 0){
                    const njoy::ENDFtk::tree::Tape tape(tpid + material_lines + tend);
                    // with more than one temperature, the same MAT comes again; we keep the first
                    for(auto mat: tape.materialNumbers()){
                        for(auto material: tape.material(mat)){
                            _read_material(material, done_zais, reaction_data);
                        }
                    }
                    material_lines.clear();
                }
            }

            return reaction_data;
        }

        Eigen::SparseMatrix<double> construct_reaction_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape,
                                                            const core::Spectrum<double>& in_spectrum){
            return construct_reaction_matrix(read_reaction_data_from_tape(tape), in_spectrum);
        }

        Eigen::SparseMatrix<double> construct_reaction_matrix_from_file(std::string filename,
                                                            const core::Spectrum<double>& in_spectrum){
            return construct_reaction_matrix(read_reaction_data_from_file(filename), in_spectrum);
        }

    }
}
