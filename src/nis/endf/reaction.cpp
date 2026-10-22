#include <cmath>
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
        // 1= y is constant in x (constant, histogram)
        constexpr int histogram_interpolation = 1;

        bool is_reaction_mt(int mt_number){
            return (mt_number == FISSION_MT) || (REACTION_MAP.count(mt_number) == 1);
        }

        ReactionProductT compute_reaction_products(int parent_zai, int mt_number){
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

            return ReactionProductT(zai::make(z, a, 0), change.secondaries);
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

        void add_reaction_rate(Eigen::MatrixXd& reaction_mat, const std::map<int, int>& zai_map,
                               int parent_zai, int mt_number, double rate){
            if((zai_map.count(parent_zai) == 0) || (rate <= 0.)){
                return;
            }
            const int za_index = zai_map.at(parent_zai);

            // the diagonal is the rate lost from the parent
            reaction_mat(za_index, za_index) -= rate;

            const ReactionProductT products = compute_reaction_products(parent_zai, mt_number);
            std::vector<int> daughters = std::get<1>(products);
            daughters.emplace_back(std::get<0>(products));
            for(auto daughter_zai: daughters){
                // the sink is the last row
                const int dza_index = zai_map.count(daughter_zai) == 1 ? zai_map.at(daughter_zai) : data::SINK_INDEX;
                reaction_mat(za_index, dza_index) += rate;
            }
        }

        void _update_reaction_matrix(const njoy::ENDFtk::tree::Material& material,
                                     const std::map<int, int>& zai_map, Eigen::MatrixXd& reaction_mat,
                                     const core::Spectrum<double>& in_spectrum){

            if(!material.hasMF(file_number)){
                return;
            }

            const int parent_zai = get_parent_zai(material);
            if(zai_map.count(parent_zai) == 0){
                return;
            }

            // TODO: MF=10 can give isomer production when MF=3 has no data for it
            njoy::ENDFtk::tree::File file = material.file(file_number);
            for(auto mt_number: file.sectionNumbers()){
                // we skip sums, elastic, inelastic levels and anything we can't map
                if(!is_reaction_mt(mt_number)){
                    continue;
                }

                njoy::ENDFtk::tree::Section section = file.MT(mt_number);
                auto parsed_section = section.parse<file_number>();

                // we assume everything is already preprocessed into INT=1 (group-wise)
                for(auto interpolation_scheme: parsed_section.interpolants()){
                    if(interpolation_scheme != histogram_interpolation){
                        throw std::runtime_error("Only INT=1 (histogram) cross sections are supported, found INT="
                                                 + std::to_string(interpolation_scheme)
                                                 + " for ZAI=" + std::to_string(parent_zai)
                                                 + " MT=" + std::to_string(mt_number));
                    }
                }

                const std::vector<double> xs_energies{parsed_section.energies().begin(), parsed_section.energies().end()};
                const std::vector<double> xs_values{parsed_section.crossSections().begin(), parsed_section.crossSections().end()};

                const double rate = collapse_cross_section(in_spectrum, xs_energies, xs_values);
                add_reaction_rate(reaction_mat, zai_map, parent_zai, mt_number, rate);
            }
        }

        Eigen::MatrixXd construct_reaction_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape,
                                                            const core::Spectrum<double>& in_spectrum){

            const int mat_nr_rows = data::ZAI_SIZE;
            const int mat_nr_cols = data::ZAI_SIZE;

            Eigen::MatrixXd reaction_mat = Eigen::MatrixXd::Zero(mat_nr_rows, mat_nr_cols);

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
                    if(!material.hasMF(file_number)){
                        continue;
                    }
                    const int parent_zai = get_parent_zai(material);
                    if(done_zais.count(parent_zai) == 0){
                        _update_reaction_matrix(material, data::ZAI_MAP, reaction_mat, in_spectrum);
                        done_zais.insert(parent_zai);
                    }
                    break;
                }
            }

            return reaction_mat;
        }

        Eigen::MatrixXd construct_reaction_matrix_from_file(std::string filename,
                                                            const core::Spectrum<double>& in_spectrum){
            auto tape = njoy::ENDFtk::tree::fromFile(filename);
            return endf::construct_reaction_matrix_from_tape(tape, in_spectrum);
        }

    }
}
