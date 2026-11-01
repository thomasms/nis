#ifndef NIS_SOLVER_HPP
#define NIS_SOLVER_HPP

#include <cassert>
#include <functional>
#include <optional>
#include <stdexcept>
#include <vector>

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
        using solver_func = std::function<Eigen::VectorX<T>(const Eigen::SparseMatrix<T>&, const Eigen::VectorX<T>&, double)>;

        // The indices that can get atoms from the non-zero entries of n0.
        // Atoms only move from a parent to a daughter where mat(parent, daughter) != 0,
        // so every other nuclide stays exactly zero, and we can leave it out of the solve.
        template<typename T>
        std::vector<int> reachable_indices(const Eigen::SparseMatrix<T>& mat, const Eigen::VectorX<T>& n0){
            assert(mat.rows() == n0.size());
            // the daughters of a parent are in its row, so use row major storage
            const Eigen::SparseMatrix<T, Eigen::RowMajor> by_parent = mat;
            std::vector<char> reached(n0.size(), 0);
            std::vector<int> to_visit;
            for(int i=0; i<n0.size(); ++i){
                if(n0(i) != T(0)){
                    reached[i] = 1;
                    to_visit.push_back(i);
                }
            }
            while(!to_visit.empty()){
                const int parent = to_visit.back();
                to_visit.pop_back();
                for(typename Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(by_parent, parent); it; ++it){
                    if((it.value() != T(0)) && !reached[it.col()]){
                        reached[it.col()] = 1;
                        to_visit.push_back(it.col());
                    }
                }
            }
            std::vector<int> indices;
            for(int i=0; i<n0.size(); ++i){
                if(reached[i]){
                    indices.push_back(i);
                }
            }
            return indices;
        }

        // the rows and columns of mat at the indices, in the same order
        template<typename T>
        Eigen::SparseMatrix<T> submatrix(const Eigen::SparseMatrix<T>& mat, const std::vector<int>& indices){
            std::vector<int> position(mat.rows(), -1);
            for(int i=0; i<static_cast<int>(indices.size()); ++i){
                position[indices[i]] = i;
            }
            std::vector<Eigen::Triplet<T>> triplets;
            for(int col=0; col<mat.outerSize(); ++col){
                if(position[col] < 0){
                    continue;
                }
                for(typename Eigen::SparseMatrix<T>::InnerIterator it(mat, col); it; ++it){
                    if(position[it.row()] >= 0){
                        triplets.emplace_back(position[it.row()], position[col], it.value());
                    }
                }
            }
            Eigen::SparseMatrix<T> sub(indices.size(), indices.size());
            sub.setFromTriplets(triplets.begin(), triplets.end());
            return sub;
        }

        // solves only for the nuclides that can be reached from n0, the rest stay zero
        template<typename T>
        Eigen::VectorX<T> solve_reachable(const Eigen::SparseMatrix<T>& mat, const Eigen::VectorX<T>& n0,
                                          double dt, const solver_func<T>& solver){
            const std::vector<int> indices = reachable_indices<T>(mat, n0);
            if(indices.empty()){
                return n0;
            }
            if(static_cast<int>(indices.size()) == n0.size()){
                return solver(mat, n0, dt);
            }
            Eigen::VectorX<T> n0_reduced(indices.size());
            for(int i=0; i<static_cast<int>(indices.size()); ++i){
                n0_reduced(i) = n0(indices[i]);
            }
            const Eigen::VectorX<T> nt_reduced = solver(submatrix<T>(mat, indices), n0_reduced, dt);
            Eigen::VectorX<T> nt = Eigen::VectorX<T>::Zero(n0.size());
            for(int i=0; i<static_cast<int>(indices.size()); ++i){
                nt(indices[i]) = nt_reduced(i);
            }
            return nt;
        }

        template<typename T>
        inline core::ModelResult<T> compute_inventory(const core::BasicModel<T>& model,
                                                      const Eigen::SparseMatrix<T>& mat,
                                                      const solver_func<T>& solver){
            const Eigen::VectorX<T> n0 = io::helpers::vector_from_inventory<T>(model.nuclides);
            const Eigen::VectorX<T> nt = solve_reachable<T>(mat, n0, model.dt, solver);
            // we should make sure that atoms are conserved
            // assert(n0.sum() == nt.sum());
            const auto final_inv = io::helpers::inventory_from_vector<T>(nt);
            core::ModelResult<T> result;
            result.nuclides = final_inv;
            return result;
        }

        // the steps of the model, an input with no dt line is one step of dt
        template<typename T>
        std::vector<core::TimeStep<T>> time_steps(const core::BasicModel<T>& model){
            if(model.steps.empty()){
                return {core::TimeStep<T>{model.dt, std::nullopt}};
            }
            return model.steps;
        }

        // the factor to scale the reaction matrix by for a step. The reaction matrix is made from
        // the spectrum as it is given (n/cm^2/s, with a total of spectrum_total), the step flux is in n/m^2/s.
        // No flux uses the spectrum as it is given (1), a flux of 0 is a cooling step (0).
        template<typename T>
        T flux_scale(const core::TimeStep<T>& step, T spectrum_total){
            if(!step.flux.has_value()){
                return T(1);
            }
            if(step.flux.value() <= T(0)){
                return T(0);
            }
            if(spectrum_total <= T(0)){
                throw std::invalid_argument("A step has a flux, so the input needs a spectrum with a total above 0.\n");
            }
            return step.flux.value() * T(1e-4) / spectrum_total;
        }

        // Solves each step of the model, each from the result of the step before.
        // The matrix of a step is decay_mat + reaction_mat * flux_scale(step, spectrum_total).
        // An empty reaction_mat (0 x 0) is no reactions, then all the steps are decay only.
        template<typename T>
        std::vector<core::ModelResult<T>> compute_steps(const core::BasicModel<T>& model,
                                                        const Eigen::SparseMatrix<T>& decay_mat,
                                                        const Eigen::SparseMatrix<T>& reaction_mat,
                                                        T spectrum_total,
                                                        const solver_func<T>& solver){
            const bool has_reactions = reaction_mat.size() > 0;
            if(has_reactions && ((reaction_mat.rows() != decay_mat.rows()) || (reaction_mat.cols() != decay_mat.cols()))){
                throw std::invalid_argument("The reaction and decay matrices must have the same size.\n");
            }

            std::vector<core::ModelResult<T>> results;
            core::BasicModel<T> step_model;
            step_model.nuclides = model.nuclides;
            for(const auto& step: time_steps(model)){
                const T scale = flux_scale(step, spectrum_total);
                step_model.dt = step.dt;
                if(has_reactions && (scale > T(0))){
                    const Eigen::SparseMatrix<T> mat = decay_mat + reaction_mat * scale;
                    results.push_back(compute_inventory<T>(step_model, mat, solver));
                }
                else{
                    results.push_back(compute_inventory<T>(step_model, decay_mat, solver));
                }
                step_model.nuclides = results.back().nuclides;
            }
            return results;
        }

    }
}

#endif