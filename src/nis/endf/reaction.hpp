#ifndef ENDF_REACTION_HPP
#define ENDF_REACTION_HPP

#include <functional>
#include <map>
#include <memory>

#include <ENDFtk.hpp>

#include <Eigen/Core>

#include "nis/spectrum.hpp"

namespace nis {
    namespace endf {

        Eigen::MatrixXd construct_reaction_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape,
                                                            const core::Spectrum<double>& in_spectrum);
        // convience function to avoid calling ENDFtk fromFile
        Eigen::MatrixXd construct_reaction_matrix_from_file(std::string filename,
                                                            const core::Spectrum<double>& in_spectrum);
    }
}

#endif // ENDF_REACTION_HPP