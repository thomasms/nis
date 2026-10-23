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

        bool is_reaction_mt(int mt_number){
            return (mt_number == FISSION_MT) || (REACTION_MAP.count(mt_number) == 1);
        }

        ReactionProductT compute_reaction_products(int parent_zai, int mt_number, int final_state){
            if(mt_number == FISSION_MT){
                // TODO: add fission yields, for now all goes to the sink
                return ReactionProductT(-1, {});
            }

            if(REACTION_MAP.count(mt_number) == 0){
                throw std::invalid_argument("MT=" + std::to_string(mt_number) + " is not a supported reaction");
            }

            const ReactionChange& change = REACTION_MAP.at(mt_number);
            const zai::triple parent_triple = zai::decompose(parent_zai);
            const int z = std::get<0>(parent_triple) + change.dz;
            const int a = std::get<1>(parent_triple) + change.da;

            // unphysical residual, i.e. (n,a) on a very light target
            if((z < 0) || (a < z) || (a <= 0)){
                return ReactionProductT(-1, change.secondaries);
            }

            return ReactionProductT(zai::make(z, a, final_state), change.secondaries);
        }

        std::map<int, double> split_rate_by_state(double total_rate, const std::map<int, double>& state_rates){
            double states_sum = 0.;
            for(const auto& state_rate: state_rates){
                states_sum += state_rate.second;
            }

            // no split, all to the ground state
            if(states_sum <= 0.){
                return {{0, total_rate}};
            }

            std::map<int, double> split;
            const bool has_ground = state_rates.count(0) == 1;
            if(!has_ground && (states_sum <= total_rate)){
                // the isomers as given, the rest to the ground state
                split = state_rates;
                split[0] = total_rate - states_sum;
                return split;
            }

            // use MF=10 for the fractions only, so the total stays the same
            for(const auto& state_rate: state_rates){
                split[state_rate.first] = total_rate*state_rate.second/states_sum;
            }
            return split;
        }

        double collapse_cross_section(const core::Spectrum<double>& in_spectrum,
                                      const std::vector<double>& xs_energies,
                                      const std::vector<double>& xs_values){
            if((xs_energies.size() < 2) || (xs_energies.size() != xs_values.size())){
                throw std::invalid_argument("Cross section must have at least 2 points and the same number of energies and values");
            }

            // because it is INT=1 (constant) the last value has no bin, so we drop it
            const std::vector<double> group_values(xs_values.begin(), xs_values.end()-1);

            const core::Spectrum<double> data_spectrum = core::convert_by_lethargy<double>(in_spectrum, xs_energies);
            return maths::dot_product<double>(data_spectrum.y(), group_values)*BARN_TO_CM2;
        }

        void add_reaction_rate(Eigen::SparseMatrix<double>& reaction_mat, const std::map<int, int>& zai_map,
                               int parent_zai, int mt_number, double rate, int final_state){
            if((zai_map.count(parent_zai) == 0) || (rate <= 0.)){
                return;
            }
            const int za_index = zai_map.at(parent_zai);

            // the diagonal is the rate lost from the parent
            reaction_mat.coeffRef(za_index, za_index) -= rate;

            const ReactionProductT products = compute_reaction_products(parent_zai, mt_number, final_state);
            std::vector<int> daughters = std::get<1>(products);
            daughters.emplace_back(std::get<0>(products));
            for(auto daughter_zai: daughters){
                // the sink is the last row
                const int dza_index = zai_map.count(daughter_zai) == 1 ? zai_map.at(daughter_zai) : data::SINK_INDEX;
                reaction_mat.coeffRef(za_index, dza_index) += rate;
            }
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

        void _read_isomer_data(const njoy::ENDFtk::tree::Material& material, int parent_zai,
                               int mt_number, ReactionData& reaction_data){
            if(!material.hasMF(isomer_file_number)){
                return;
            }
            njoy::ENDFtk::tree::File file = material.file(isomer_file_number);
            if(!file.hasMT(mt_number)){
                return;
            }

            // the residual we expect for this MT, MF=10 can also list other products
            const int residual_zai = std::get<0>(compute_reaction_products(parent_zai, mt_number));
            if(residual_zai < 0){
                return;
            }
            const int residual_za = residual_zai/10;

            njoy::ENDFtk::tree::Section section = file.MT(mt_number);
            auto parsed_section = section.parse<isomer_file_number>();
            for(const auto& product: parsed_section.reactionProducts()){
                if(product.IZAP() != residual_za){
                    continue;
                }
                _check_histogram(product, parent_zai, isomer_file_number, mt_number);
                // NOTE: we use the level number (LFS) as the isomeric state, as TENDL and JEFF-4 do
                reaction_data.push_back(CrossSection{
                    parent_zai,
                    mt_number,
                    {product.energies().begin(), product.energies().end()},
                    {product.crossSections().begin(), product.crossSections().end()},
                    product.LFS()
                });
            }
        }

        void _read_reaction_data(const njoy::ENDFtk::tree::Material& material,
                                 const std::map<int, int>& zai_map, ReactionData& reaction_data){

            const int parent_zai = get_parent_zai(material);
            if(zai_map.count(parent_zai) == 0){
                return;
            }

            // TODO: MF=10 can give isomer production when MF=3 has no data for it
            // TODO: MF=9 (multiplicities) is not used, JEFF-4 only has MF=10
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

                _read_isomer_data(material, parent_zai, mt_number, reaction_data);
            }
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

        ReactionData read_reaction_data_from_file(std::string filename){
            // A reaction tape can be ~1 GB, and parsing it all at once needs ~2x that.
            // So we read one material at a time and parse it as a small tape of its own.
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
                material_lines += line;
                material_lines += '\n';
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

        Eigen::SparseMatrix<double> construct_reaction_matrix(const ReactionData& reaction_data,
                                                  const core::Spectrum<double>& in_spectrum){
            // collapse all first, then split each reaction by residual state
            struct Rates {
                double total = 0.;
                std::map<int, double> states;
            };
            std::map<std::pair<int, int>, Rates> rates;
            for(const auto& xs: reaction_data){
                // the data can be made by hand, so check again
                if(!is_reaction_mt(xs.mt_number)){
                    continue;
                }
                const double rate = collapse_cross_section(in_spectrum, xs.energies, xs.values);
                Rates& r = rates[{xs.parent_zai, xs.mt_number}];
                if(xs.final_state == ALL_STATES){
                    r.total += rate;
                }
                else{
                    r.states[xs.final_state] += rate;
                }
            }

            Eigen::SparseMatrix<double> reaction_mat(data::ZAI_SIZE, data::ZAI_SIZE);
            for(const auto& [key, r]: rates){
                const int parent_zai = key.first;
                const int mt_number = key.second;
                for(const auto& [state, rate]: split_rate_by_state(r.total, r.states)){
                    // if we don't track this isomer, the atoms go to the ground state
                    const int residual = std::get<0>(compute_reaction_products(parent_zai, mt_number, state));
                    const int final_state = (residual >= 0) && (data::ZAI_MAP.count(residual) == 0) ? 0 : state;
                    add_reaction_rate(reaction_mat, data::ZAI_MAP, parent_zai, mt_number, rate, final_state);
                }
            }
            reaction_mat.makeCompressed();
            return reaction_mat;
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
