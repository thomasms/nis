#ifndef ENDF_DECAY_HPP
#define ENDF_DECAY_HPP

#include <string>

#include <ENDFtk.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "nis/decayheat.hpp"
#include "nis/data/decay.hpp"
#include "nis/endf/fission_yields.hpp"

// Reads the decay data from ENDF (MF=8 MT=457) into data::DecayData.
// The matrix and the decay energies are made from data::DecayData, see nis/data/decay.hpp.

namespace nis {
    namespace endf {

        // the names before the data layer was made, so old code still compiles
        using data::DecayProductT;
        using data::DecayProductComputer;
        using data::DECAY_MAP;
        using data::expand_rtypes;
        using data::compute_products;
        using data::SPONTANEOUS_FISSION_RTYP;

        // all the nuclides of the tape, also the stable ones and the ones not in ZAI_MAP
        data::DecayData read_decay_data_from_tape(const njoy::ENDFtk::tree::Tape& tape);
        data::DecayData read_decay_data_from_file(const std::string& filename);

        // the matrix is sparse, (parent, daughter), with the sink as the last index.
        // Same as data::construct_decay_matrix(read_decay_data_from_tape(tape), sf_yields).
        Eigen::SparseMatrix<double> construct_decay_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape,
                                                                     const FissionYieldData* sf_yields = nullptr);
        // convience function to avoid calling ENDFtk fromFile
        Eigen::SparseMatrix<double> construct_decay_matrix_from_file(std::string filename,
                                                                     const FissionYieldData* sf_yields = nullptr);

        // decay constants and average decay energies for decay heat
        // Same as data::construct_decay_energy_data(read_decay_data_from_tape(tape)).
        core::DecayEnergyData<double> read_decay_energy_data_from_tape(const njoy::ENDFtk::tree::Tape& tape);
        // convience function to avoid calling ENDFtk fromFile
        core::DecayEnergyData<double> read_decay_energy_data_from_file(std::string filename);


    }
}

#endif //ENDF_DECAY_HPP
