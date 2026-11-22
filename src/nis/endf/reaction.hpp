#ifndef ENDF_REACTION_HPP
#define ENDF_REACTION_HPP

#include <map>
#include <set>
#include <string>
#include <utility>

#include <ENDFtk.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "nis/spectrum.hpp"
#include "nis/data/reaction.hpp"
#include "nis/endf/fission_yields.hpp"

// Reads the cross sections from ENDF (MF=3 and MF=10) into data::ReactionData.
// The reaction matrix is made from data::ReactionData, see nis/data/reaction.hpp.

namespace nis {
    namespace endf {

        // the names before the data layer was made, so old code still compiles
        using data::ReactionProductT;
        using data::PROTON;
        using data::DEUTERON;
        using data::TRITON;
        using data::HELION;
        using data::BARN_TO_CM2;
        using data::FISSION_MT;
        using data::INELASTIC_MT;
        using data::ANYTHING_MT;
        using data::is_light_particle;
        using data::ReactionChange;
        using data::REACTION_MAP;
        using data::is_reaction_mt;
        using data::compute_reaction_products;
        using data::collapse_cross_section;
        using data::add_reaction_rate;
        using data::ALL_STATES;
        using data::CrossSection;
        using data::split_rate_by_state;
        using data::ReactionData;
        using data::construct_reaction_matrix;

        // The level number (LFS) in MF=10 is the excited level of the product, not the isomer number.
        // So for each product ZA we sort the levels in MF=10 of the material (but not MT=5,
        // where JEFF-4 and TENDL give the isomer number as the LFS):
        // LFS=0 is the ground state, the next level is the 1st isomer, and so on.
        // key is (product ZA, LFS), value is the isomeric state.
        std::map<std::pair<int, int>, int> isomer_states_from_levels(const std::map<int, std::set<int>>& levels);

        ReactionData read_reaction_data_from_tape(const njoy::ENDFtk::tree::Tape& tape);
        // reads the file one material at a time, so the whole tape is never in memory
        ReactionData read_reaction_data_from_file(std::string filename);

        Eigen::SparseMatrix<double> construct_reaction_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape,
                                                            const core::Spectrum<double>& in_spectrum);
        // convience function, reads the file with read_reaction_data_from_file
        Eigen::SparseMatrix<double> construct_reaction_matrix_from_file(std::string filename,
                                                            const core::Spectrum<double>& in_spectrum);
    }
}

#endif // ENDF_REACTION_HPP
