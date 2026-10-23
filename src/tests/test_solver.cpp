#include <algorithm>
#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Sparse>

#include "nis/cram.hpp"
#include "nis/solver.hpp"

using namespace nis;

using Catch::Matchers::WithinRel;

namespace {

    // two chains that do not mix: 0 -> 1 -> 2 (stable) and 3 -> 4 (stable)
    Eigen::SparseMatrix<double> two_chains(){
        const std::vector<Eigen::Triplet<double>> triplets = {
            {0, 0, -1e-3}, {0, 1, 1e-3},
            {1, 1, -2e-4}, {1, 2, 2e-4},
            {3, 3, -5e-3}, {3, 4, 5e-3},
        };
        Eigen::SparseMatrix<double> mat(5, 5);
        mat.setFromTriplets(triplets.begin(), triplets.end());
        return mat;
    }

    solver::solver_func<double> cram = &core::cram48::solve<double>;
}

TEST_CASE( "reachable_indices follows parent to daughter" ) {
    const auto mat = two_chains();
    Eigen::VectorXd n0 = Eigen::VectorXd::Zero(5);
    n0(0) = 1e20;
    REQUIRE( solver::reachable_indices<double>(mat, n0) == std::vector<int>{0, 1, 2} );

    // a daughter does not reach its parent
    n0 = Eigen::VectorXd::Zero(5);
    n0(1) = 1e20;
    REQUIRE( solver::reachable_indices<double>(mat, n0) == std::vector<int>{1, 2} );

    n0(4) = 1.0;
    REQUIRE( solver::reachable_indices<double>(mat, n0) == std::vector<int>{1, 2, 4} );

    REQUIRE( solver::reachable_indices<double>(mat, Eigen::VectorXd::Zero(5)).empty() );
}

TEST_CASE( "submatrix keeps the rows and columns in order" ) {
    const auto sub = solver::submatrix<double>(two_chains(), {1, 2});
    REQUIRE( sub.rows() == 2 );
    REQUIRE( sub.coeff(0, 0) == -2e-4 );
    REQUIRE( sub.coeff(0, 1) == 2e-4 );
    REQUIRE( sub.nonZeros() == 2 );
}

TEST_CASE( "solve_reachable gives the same as the full solve" ) {
    const auto mat = two_chains();
    for(int start: {0, 1, 3}){
        Eigen::VectorXd n0 = Eigen::VectorXd::Zero(5);
        n0(start) = 1e20;
        for(double dt: {1.0, 1e3, 1e5}){
            const Eigen::VectorXd full = cram(mat, n0, dt);
            const Eigen::VectorXd reduced = solver::solve_reachable<double>(mat, n0, dt, cram);
            const auto reachable = solver::reachable_indices<double>(mat, n0);
            for(int i=0; i<5; ++i){
                // the full solve has rounding noise of ~1e-16 of the total atoms,
                // also for the nuclides it cannot reach
                REQUIRE( std::abs(reduced(i) - full(i)) <= 1e-12*n0.sum() );
                if(std::find(reachable.begin(), reachable.end(), i) == reachable.end()){
                    REQUIRE( reduced(i) == 0.0 );
                }
                else if(full(i) > 1e-6*n0.sum()){
                    REQUIRE_THAT( reduced(i), WithinRel(full(i), 1e-10) );
                }
            }
        }
    }
}

TEST_CASE( "solve_reachable with no atoms" ) {
    const Eigen::VectorXd n0 = Eigen::VectorXd::Zero(5);
    REQUIRE( solver::solve_reachable<double>(two_chains(), n0, 1e3, cram).isZero() );
}
