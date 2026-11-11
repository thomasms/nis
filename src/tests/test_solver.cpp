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

namespace {

    // A long-lived parent P feeds a short-lived daughter D (half-life 17.6 ms, like As75m) with a small
    // branching ratio, and D decays to a stable nuclide S with very many atoms (like As75).
    // A very short-lived nuclide Q (lambda ~2e16, like Be8) decays to R. So lambda*dt of D is about
    // 600 to 1200 and that of Q about 1e17, and the atoms of the nuclides go from ~1 to 6e21.
    struct EquilibriumChain {
        static constexpr int P = 0, D = 1, S = 2, Q = 3, R = 4, SIZE = 5;
        static constexpr double lambda_p = 1.4531e-2;   // like Ge75m, half-life 47.7 s
        static constexpr double branching = 3e-4;      // P -> D, the rest of P goes to S
        static constexpr double lambda_d = 39.3387;     // like As75m, half-life 17.6 ms
        static constexpr double lambda_q = 1.98e16;     // like Be8

        // the nuclides in the matrix are in the order of the indices in order
        static Eigen::SparseMatrix<double> matrix(const std::vector<int>& order){
            const std::vector<Eigen::Triplet<double>> triplets = {
                {P, P, -lambda_p}, {P, D, branching*lambda_p}, {P, S, (1.0 - branching)*lambda_p},
                {D, D, -lambda_d}, {D, S, lambda_d},
                {Q, Q, -lambda_q}, {Q, R, 2.0*lambda_q},
            };
            std::vector<int> position(SIZE);
            for(int i=0; i<SIZE; ++i){
                position[order[i]] = i;
            }
            std::vector<Eigen::Triplet<double>> ordered;
            for(const auto& t: triplets){
                ordered.emplace_back(position[t.row()], position[t.col()], t.value());
            }
            Eigen::SparseMatrix<double> mat(SIZE, SIZE);
            mat.setFromTriplets(ordered.begin(), ordered.end());
            return mat;
        }

        static Eigen::VectorXd atoms(const std::vector<int>& order){
            Eigen::VectorXd n0 = Eigen::VectorXd::Zero(SIZE);
            const std::vector<double> values = {5.7e7, 0.0, 6.08792e21, 1.0, 8.86e9};
            for(int i=0; i<SIZE; ++i){
                n0(i) = values[order[i]];
            }
            return n0;
        }

        // the exact D at t, from D(0) = 0: in equilibrium with P after a few ms
        static double daughter(double t){
            return branching*lambda_p*atoms({0, 1, 2, 3, 4})(P)/(lambda_d - lambda_p)
                   * (std::exp(-lambda_p*t) - std::exp(-lambda_d*t));
        }
    };
}

TEST_CASE( "a short-lived daughter in equilibrium with a long-lived parent" ) {
    using C = EquilibriumChain;
    // all the orders of the nuclides in the matrix give the same result
    std::vector<int> order = {0, 1, 2, 3, 4};
    do {
        const auto mat = C::matrix(order);
        const auto n0 = C::atoms(order);
        int d = 0;
        while(order[d] != C::D){
            ++d;
        }
        for(double dt: {15.96, 16.02, 27.0, 33.96}){
            const double expected = C::daughter(dt);
            // the daughter is ~0.16 to 0.3 atoms in 6e21 atoms
            REQUIRE( dt*C::lambda_d > 600.0 );
            for(double factor: {1.0, 1.0 + 1e-12, 1.0 - 1e-12}){
                const Eigen::VectorXd full = cram(mat, n0, dt*factor);
                const Eigen::VectorXd reduced = solver::solve_reachable<double>(mat, n0, dt*factor, cram);
                REQUIRE_THAT( full(d), WithinRel(expected, 1e-9) );
                REQUIRE_THAT( reduced(d), WithinRel(expected, 1e-9) );
            }
        }
    } while(std::next_permutation(order.begin(), order.end()));
}
