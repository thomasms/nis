#include <cmath>
#include <memory>
#include <iostream>

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

        // TODO: of course we need to add the spectrum so we can do the collapse
        // OR we could instead make an array of matricies, one per energy bin - much more memory intensive but
        // no dependence on the flux
        // this would mean 709*3352*3352*8/(1024**3) = 60 GB memory for 709 group (assuming non-sparse storage), this is not feasible
        // in practice the matricies are very sparse so it would be much lower but still above 1GB in memory.
        void _update_reaction_matrix(const njoy::ENDFtk::tree::Material& material,
                                     const std::map<int, int>& zai_map, Eigen::MatrixXd& reaction_mat,
                                     const core::Spectrum<double>& in_spectrum){

            const int parent_zai = get_parent_zai(material);
            // how do we deal with file 10?
            // do we need to look in file 10 if not in file 3?
            if(material.hasMF(file_number)){
                njoy::ENDFtk::tree::File file = material.file(file_number);
                // loop over MT numbers
                // need a lookup of MT numbers and reaction products
                // we ignore 200's because these are redundant i.e 201 is (z,Xn) total neutron production
                // we also ignore 300's because these are total heat production
                // same for 400's because these are kermas - we just want reaction data to make the graph/matrix
                // 301-450: Energy release parameters, for total and partial cross sections;
                // MT= 300 plus the reaction MT number, e.g., MT=302 is the elastic scattering kerma.
                for(auto mt_number: file.sectionNumbers()){
                    njoy::ENDFtk::tree::Section section = file.MT(mt_number);
                    auto parsed_section = section.parse<file_number>();
                    if (zai_map.count(parent_zai) == 1){
                        const int za_index = zai_map.at(parent_zai);
                        // we don't care about anything higher than 199 for neutrons
                        // actually 152-200 is unassigned so we can skip below 152
                        if(mt_number < 152){
                            // map MTs to ReactionProduct like what we do for the DecayProduct

                            // TODO: should we check if the total MT=1 is actually the sum of the others?
                            // It's commonly not the case for many nuclear data libraries

                            std::vector<int> interpolation_types{parsed_section.interpolants().begin(), parsed_section.interpolants().end()};
                            // NOTE: we assume 1 interpolation type - is this valid?
                            assert(interpolation_types.size() == 1);
                            const int interpolation_scheme = interpolation_types[0];
                            // NOTE: 1= y is constant in x (constant, histogram)
                            // we assume everything is already preprocessed into INT=1
                            assert(interpolation_scheme == 1);

                            std::vector<double> xs_energies{parsed_section.energies().begin(), parsed_section.energies().end()};
                            std::vector<double> xs_values{parsed_section.crossSections().begin(), parsed_section.crossSections().end()};
                            // HACK: because it is INT=1 (constant) we can just pop off the last entry
                            xs_values.pop_back();

                            core::Spectrum<double> data_spectrum = core::convert_by_lethargy<double>(in_spectrum, xs_energies);
                            // TODO: check units
                            double collapsed_xs = maths::dot_product<double>(data_spectrum.y(), xs_values);

                            // TODO: finish
                            // the diagonal is the production rate
                            // if(mt_number == 1){
                            //     reaction_mat(za_index, za_index) += collapsed_xs;
                            // }
                        }
                    }
                }
            }
        }

        Eigen::MatrixXd construct_reaction_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape,
                                                            const core::Spectrum<double>& in_spectrum){

            const int mat_nr_rows = data::ZAI_SIZE;
            const int mat_nr_cols = data::ZAI_SIZE;

            Eigen::MatrixXd reaction_mat = Eigen::MatrixXd::Zero(mat_nr_rows, mat_nr_cols);

            // loop over all materials in the file, most likely just one
            for(auto mat_number: tape.materialNumbers()){
                // this can return multiple instances of the same MAT number
                // due to different temperatures
                // loop over temperatures for each MAT
                //
                // How do we account for multiple temperatures?
                // We don't want to double count...
                auto materials = tape.material(mat_number);
                for(auto material: materials){
                    _update_reaction_matrix(material, data::ZAI_MAP, reaction_mat, in_spectrum);
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
