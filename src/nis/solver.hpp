#ifndef NIS_SOLVER_HPP
#define NIS_SOLVER_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "nis/io/input.hpp"

namespace nis {
    namespace solver {

        Eigen::VectorX<double> compute_inventory(const io::input::InputData<double>& input_data, const Eigen::MatrixXd& mat);

    }
}

#endif