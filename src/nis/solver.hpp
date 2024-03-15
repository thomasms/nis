#ifndef NIS_SOLVER_HPP
#define NIS_SOLVER_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "nis/cram.hpp"
#include "nis/io/input.hpp"
#include "nis/io/output.hpp"
#include "nis/io/helpers.hpp"

namespace nis {
    namespace solver {

        template<typename T>
        inline io::output::OutputData<T> compute_inventory(const io::input::InputData<T>& input_data, const Eigen::MatrixX<T>& mat){
            const Eigen::VectorX<T> n0 = io::helpers::vector_from_inventory<T>(input_data.nuclides);
            const Eigen::VectorX<T> nt = core::cram48::solve<T>(mat, n0, input_data.dt);
            io::output::OutputData<T> output_data;
            const auto final_inv = io::helpers::inventory_from_vector<T>(nt);
            output_data.nuclides = final_inv;
            return output_data;
        }

    }
}

#endif