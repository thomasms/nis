#ifndef ENDF_FISSION_YIELDS_HPP
#define ENDF_FISSION_YIELDS_HPP

#include <string>

#include <ENDFtk.hpp>

#include "nis/data/fission_yields.hpp"

// Reads the fission yields from ENDF. The data types and the maths are in nis/data/fission_yields.hpp.

namespace nis {
    namespace endf {

        // the names before the data layer was made, so old code still compiles
        using data::FissionYields;
        using data::FissionYieldData;
        using data::energy_weights;
        using data::yields_at_energy;
        using data::weighted_yields;

        // MF=8 MT=454: the independent fission product yields. These are the yields before any decay,
        // the decay matrix makes the rest of each chain, so we do not use the cumulative yields (MT=459).
        constexpr int FISSION_YIELD_FILE = 8;
        constexpr int INDEPENDENT_YIELDS_MT = 454;

        // reads MF=8 MT=454 of each material, the parent isomeric state is the LISO in MF=1 MT=451
        FissionYieldData read_fission_yields_from_tape(const njoy::ENDFtk::tree::Tape& tape);
        FissionYieldData read_fission_yields_from_file(const std::string& filename);
    }
}

#endif // ENDF_FISSION_YIELDS_HPP
