#ifndef NIS_SOLVER_HPP
#define NIS_SOLVER_HPP

#include <cassert>
#include <functional>
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
                                          double dt, solver_func<T>& solver){
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
                                                      solver_func<T>& solver){
            const Eigen::VectorX<T> n0 = io::helpers::vector_from_inventory<T>(model.nuclides);
            const Eigen::VectorX<T> nt = solve_reachable<T>(mat, n0, model.dt, solver);
            // we should make sure that atoms are conserved
            // assert(n0.sum() == nt.sum());
            const auto final_inv = io::helpers::inventory_from_vector<T>(nt);
            core::ModelResult<T> result;
            result.nuclides = final_inv;
            return result;
        }

    }
}

#endif