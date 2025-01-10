#ifndef NIS_SOLVER_HPP
#define NIS_SOLVER_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "nis/cram.hpp"
#include "nis/model.hpp"
#include "nis/io/output.hpp"
#include "nis/io/helpers.hpp"

namespace nis {
    namespace solver {

        template<typename T>
        inline core::ModelResult<T> compute_inventory(const core::BasicModel<T>& model, const Eigen::MatrixX<T>& mat){
            const Eigen::VectorX<T> n0 = io::helpers::vector_from_inventory<T>(model.nuclides);
            const Eigen::VectorX<T> nt = core::cram48::solve<T>(mat, n0, model.dt);
            const auto final_inv = io::helpers::inventory_from_vector<T>(nt);
            core::ModelResult<T> result;
            result.nuclides = final_inv;
            return result;
        }

    }
}

#endif