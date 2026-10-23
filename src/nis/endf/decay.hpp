#ifndef ENDF_DECAY_HPP
#define ENDF_DECAY_HPP

#include <functional>
#include <map>
#include <memory>

#include <ENDFtk.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "nis/decayheat.hpp"
#include "nis/endf/decay/product_computer.hpp"


namespace nis {
    namespace endf {

        const std::map<int, DecayProductComputer> DECAY_MAP = {
            {0, GammaRayReactionProducts() },
            {1, BetaMinusDecayReactionProducts() },
            {2, BetaPlusDecayReactionProducts() },
            {3, ITDecayReactionProducts() },
            {4, AlphaDecayReactionProducts() },
            {5, NeutronEmissionDecayReactionProducts() },
            {6, SFDecayReactionProducts() },
            {7, ProtonEmissionDecayReactionProducts() },
            {8, BetaMinusDecayReactionProducts() },
            {9, GammaRayReactionProducts() },
            // unknown origin
            {10, ToSinkDecayReactionProducts() }
        };

        std::vector<int> expand_rtypes(double rtyp);

        DecayProductT compute_products(int parent_zai, double rtyp, int final_isomeric_state);

        // the matrix is sparse, (parent, daughter), with the sink as the last index
        Eigen::SparseMatrix<double> construct_decay_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape);
        // convience function to avoid calling ENDFtk fromFile
        Eigen::SparseMatrix<double> construct_decay_matrix_from_file(std::string filename);

        // decay constants and average decay energies for decay heat
        core::DecayEnergyData<double> read_decay_energy_data_from_tape(const njoy::ENDFtk::tree::Tape& tape);
        // convience function to avoid calling ENDFtk fromFile
        core::DecayEnergyData<double> read_decay_energy_data_from_file(std::string filename);


    }
}

#endif //ENDF_DECAY_HPP