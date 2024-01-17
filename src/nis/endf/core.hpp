#ifndef ENDF_CORE_HPP
#define ENDF_CORE_HPP

#include <ENDFtk.hpp>
#include <Eigen/Dense>

namespace nis {
    namespace endf {

        void _update_decay_matrix(const njoy::ENDFtk::tree::Material& material,
                                  const std::map<int, int>& zai_map, Eigen::MatrixXd& decay_mat);

        Eigen::MatrixXd construct_decay_matrix_from_tape(const njoy::ENDFtk::tree::Tape& tape);

    }
}

#endif //ENDF_CORE_HPP