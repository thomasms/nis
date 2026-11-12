// The accuracy of the CRAM48 solver against exact solutions, and the properties a solution must have:
// atoms are conserved, no atoms are made from nothing, the result does not depend on how the problem is put
// (the order of the nuclides, steps, nuclides with no atoms). No nuclear data is needed.
//
// The matrices are (parent, daughter), the same as the decay and reaction matrices of NIS.

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Sparse>

#include "nis/cram.hpp"
#include "nis/solver.hpp"

using namespace nis;

using Catch::Matchers::WithinRel;

namespace {

    using Triplets = std::vector<Eigen::Triplet<double>>;

    const solver::solver_func<double> cram = &core::cram48::solve<double>;

    Eigen::SparseMatrix<double> make_matrix(int size, const Triplets& triplets){
        Eigen::SparseMatrix<double> mat(size, size);
        mat.setFromTriplets(triplets.begin(), triplets.end());
        return mat;
    }

    // CRAM48 is a rational approximation of exp(-x) for x >= 0, with an error of ~2.3e-47
    // (of the atoms at the start), so a value that is exactly 0 can be ~1e-47 of the atoms at the start
    constexpr double CRAM_ERROR = 1e-45;

    // the exact solution of a linear chain 0 -> 1 -> ... -> n-1 with all the atoms in 0 at the start
    // (the Bateman equations). The decay constants must all be different, the last can be 0 (stable).
    // The sum has terms of both signs, so for a small lambda*t it loses digits. error[m] is the
    // rounding error of result[m] (long double).
    struct Bateman {
        std::vector<long double> result;
        std::vector<long double> error;
    };

    Bateman bateman(const std::vector<long double>& lambda, long double n0, long double t){
        const int n = static_cast<int>(lambda.size());
        Bateman b{std::vector<long double>(n, 0.0L), std::vector<long double>(n, 0.0L)};
        for(int m=0; m<n; ++m){
            long double production = 1.0L;
            for(int i=0; i<m; ++i){
                production *= lambda[i];
            }
            long double sum = 0.0L;
            long double sum_abs = 0.0L;
            for(int j=0; j<=m; ++j){
                long double denominator = 1.0L;
                for(int k=0; k<=m; ++k){
                    if(k != j){
                        denominator *= (lambda[k] - lambda[j]);
                    }
                }
                sum += std::exp(-lambda[j]*t)/denominator;
                sum_abs += std::abs(std::exp(-lambda[j]*t)/denominator);
            }
            b.result[m] = n0*production*sum;
            b.error[m] = 1e-18L*n0*production*sum_abs;
        }
        return b;
    }

    Eigen::SparseMatrix<double> chain_matrix(const std::vector<long double>& lambda){
        Triplets triplets;
        const int n = static_cast<int>(lambda.size());
        for(int i=0; i<n; ++i){
            if(lambda[i] > 0){
                triplets.emplace_back(i, i, -static_cast<double>(lambda[i]));
                if(i + 1 < n){
                    triplets.emplace_back(i, i + 1, static_cast<double>(lambda[i]));
                }
            }
        }
        return make_matrix(n, triplets);
    }

    // a random decay matrix with no loops: each nuclide decays to 1 to 3 of the nuclides after it,
    // with branching ratios that add up to 1. The last nuclide is stable. The half-lives are
    // log-uniform from 10^min_log10 to 10^max_log10 s.
    Eigen::SparseMatrix<double> random_decay_matrix(int size, std::mt19937& rng,
                                                    double min_log10 = -12.0, double max_log10 = 12.0){
        std::uniform_real_distribution<double> log_half_life(min_log10, max_log10);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        Triplets triplets;
        for(int i=0; i+1<size; ++i){
            const double lambda = std::log(2.0)/std::pow(10.0, log_half_life(rng));
            triplets.emplace_back(i, i, -lambda);
            std::uniform_int_distribution<int> daughter(i + 1, size - 1);
            const int n_daughters = std::min(1 + static_cast<int>(uniform(rng)*3), size - 1 - i);
            std::vector<double> weights(n_daughters);
            for(auto& w: weights){
                w = 0.01 + uniform(rng);
            }
            const double total = std::accumulate(weights.begin(), weights.end(), 0.0);
            for(int d=0; d<n_daughters; ++d){
                // two branches to the same daughter add up
                triplets.emplace_back(i, daughter(rng), lambda*weights[d]/total);
            }
        }
        return make_matrix(size, triplets);
    }

    Eigen::VectorXd random_atoms(int size, std::mt19937& rng){
        // from 1 to 1e22 atoms, and some nuclides with none
        std::uniform_real_distribution<double> log_atoms(0.0, 22.0);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        Eigen::VectorXd n0(size);
        for(int i=0; i<size; ++i){
            n0(i) = uniform(rng) < 0.3 ? 0.0 : std::pow(10.0, log_atoms(rng));
        }
        return n0;
    }

    // a and b agree to rel_tol, or the difference is below abs_tol
    bool close(double a, double b, double rel_tol, double abs_tol){
        return std::abs(a - b) <= std::max(rel_tol*std::max(std::abs(a), std::abs(b)), abs_tol);
    }
}

TEST_CASE( "CRAM48 gives the Bateman solution of a linear chain" ) {
    // half-lives from ~22 years to ~3e-17 s (like Be8), each a factor of 10 or more from the others
    const std::vector<long double> lambda = {1e-9L, 3e-6L, 1e-3L, 0.7L, 40.0L, 2e4L, 1e8L, 2e16L, 0.0L};
    const int n = static_cast<int>(lambda.size());
    const Eigen::SparseMatrix<double> mat = chain_matrix(lambda);
    const double n0_total = 1e20;
    Eigen::VectorXd n0 = Eigen::VectorXd::Zero(n);
    n0(0) = n0_total;

    int precise_checks = 0;
    int checks = 0;
    for(double dt: {1e-3, 1.0, 16.02, 1e3, 1e6, 1e9, 1e10}){
        const auto exact = bateman(lambda, n0_total, dt);
        const Eigen::VectorXd full = cram(mat, n0, dt);
        const Eigen::VectorXd reduced = solver::solve_reachable<double>(mat, n0, dt, cram);
        for(int i=0; i<n; ++i){
            const double expected = static_cast<double>(exact.result[i]);
            const double abs_tol = std::max(static_cast<double>(exact.error[i]), CRAM_ERROR*n0_total);
            INFO( "dt " << dt << ", nuclide " << i << ", exact " << expected << ", cram " << full(i)
                  << ", error of the exact value " << static_cast<double>(exact.error[i]) );
            // each nuclide to 1e-10 of its own value, it does not matter how small it is next to the others
            REQUIRE( close(full(i), expected, 1e-10, abs_tol) );
            REQUIRE( close(reduced(i), expected, 1e-10, abs_tol) );
            ++checks;
            if(abs_tol <= 1e-10*std::abs(expected)){
                ++precise_checks;
            }
        }
    }
    // most of the exact values are precise enough to check CRAM to 1e-10
    REQUIRE( precise_checks >= 0.7*checks );
}

TEST_CASE( "CRAM48 with branching: two paths to the same daughter" ) {
    // A -> B (30 %) and A -> C (70 %), B -> D, C -> D, D stable
    const long double la = 2e-4L, lb = 5e-2L, lc = 3e1L;
    const Eigen::SparseMatrix<double> mat = make_matrix(4, {
        {0, 0, -double(la)}, {0, 1, 0.3*double(la)}, {0, 2, 0.7*double(la)},
        {1, 1, -double(lb)}, {1, 3, double(lb)},
        {2, 2, -double(lc)}, {2, 3, double(lc)},
    });
    const double total = 6e21;
    Eigen::VectorXd n0 = Eigen::VectorXd::Zero(4);
    n0(0) = total;
    for(double dt: {0.01, 10.0, 1e3, 1e5}){
        const auto via_b = bateman({la, lb, 0.0L}, total, dt);
        const auto via_c = bateman({la, lc, 0.0L}, total, dt);
        const std::vector<double> exact = {
            double(via_b.result[0]), 0.3*double(via_b.result[1]), 0.7*double(via_c.result[1]),
            0.3*double(via_b.result[2]) + 0.7*double(via_c.result[2])
        };
        const Eigen::VectorXd nt = cram(mat, n0, dt);
        for(int i=0; i<4; ++i){
            INFO( "dt " << dt << ", nuclide " << i );
            const double abs_tol = std::max({double(via_b.error[2]), double(via_c.error[2]), CRAM_ERROR*total});
            REQUIRE( close(nt(i), exact[i], 1e-10, abs_tol) );
        }
        REQUIRE_THAT( nt.sum(), WithinRel(total, 1e-13) );
    }
}

TEST_CASE( "CRAM48 with a loop between two nuclides" ) {
    // A -> B at rate a, B -> A at rate b (for example a capture and a decay back), a stiff loop
    for(const auto& [a, b]: std::vector<std::pair<double, double>>{{1e-3, 50.0}, {50.0, 1e-3}, {2.0, 3.0}, {1e-8, 1e6}}){
        const Eigen::SparseMatrix<double> mat = make_matrix(2, {{0, 0, -a}, {0, 1, a}, {1, 1, -b}, {1, 0, b}});
        const Eigen::VectorXd n0 = Eigen::Vector2d(1e20, 3e5);
        for(double dt: {1e-4, 1.0, 100.0, 1e7}){
            const double total = n0.sum();
            const double decay = std::exp(-(a + b)*dt);
            // A goes to its equilibrium total*b/(a+b)
            const double exact_a = total*b/(a + b) + (n0(0) - total*b/(a + b))*decay;
            const Eigen::VectorXd nt = cram(mat, n0, dt);
            INFO( "a " << a << ", b " << b << ", dt " << dt );
            // In a loop the rounding error grows with (a+b)*dt: ~7e-10 for a = 2/s, b = 3/s and dt = 1e7 s.
            // Real loops (a capture and a decay back) have a reaction rate of ~1e-4/s or less, so this is small.
            const double rel_tol = std::max(1e-10, 1e-16*(a + b)*dt);
            REQUIRE( close(nt(0), exact_a, rel_tol, 1e-14*total) );
            REQUIRE( close(nt(1), total - exact_a, rel_tol, 1e-14*total) );
            REQUIRE_THAT( nt.sum(), WithinRel(total, std::max(1e-12, 1e-16*(a + b)*dt)) );
        }
    }
}

TEST_CASE( "CRAM48 with more than one product per decay" ) {
    // alpha decay: A -> B + He4. Fission: F -> two products, 2 atoms from 1.
    const double la = 1e-2, lf = 3e-5;
    const Eigen::SparseMatrix<double> mat = make_matrix(6, {
        {0, 0, -la}, {0, 1, la}, {0, 2, la},                          // A -> B + He4
        {3, 3, -lf}, {3, 4, 0.6*lf}, {3, 5, 1.4*lf},                  // F -> P1 (0.6) + P2 (1.4)
    });
    Eigen::VectorXd n0 = Eigen::VectorXd::Zero(6);
    n0(0) = 1e18;
    n0(3) = 5e21;
    for(double dt: {1.0, 1e3, 1e5, 1e7}){
        const Eigen::VectorXd nt = cram(mat, n0, dt);
        const double a_gone = n0(0)*(1.0 - std::exp(-la*dt));
        const double f_gone = n0(3)*(-std::expm1(-lf*dt));
        INFO( "dt " << dt );
        REQUIRE( close(nt(0), n0(0) - a_gone, 1e-10, CRAM_ERROR*n0.sum()) );
        REQUIRE( close(nt(1), a_gone, 1e-10, CRAM_ERROR*n0.sum()) );
        REQUIRE( close(nt(2), a_gone, 1e-10, CRAM_ERROR*n0.sum()) );
        REQUIRE( close(nt(3), n0(3) - f_gone, 1e-10, CRAM_ERROR*n0.sum()) );
        REQUIRE( close(nt(4), 0.6*f_gone, 1e-10, CRAM_ERROR*n0.sum()) );
        REQUIRE( close(nt(5), 1.4*f_gone, 1e-10, CRAM_ERROR*n0.sum()) );
        // each alpha decay and each fission adds one atom
        REQUIRE_THAT( nt.sum(), WithinRel(n0.sum() + a_gone + f_gone, 1e-12) );
    }
}

TEST_CASE( "CRAM48 conserves atoms and makes no negative atoms" ) {
    std::mt19937 rng(20261112);
    for(int trial=0; trial<50; ++trial){
        const int size = 5 + trial % 40;
        const Eigen::SparseMatrix<double> mat = random_decay_matrix(size, rng);
        const Eigen::VectorXd n0 = random_atoms(size, rng);
        if(n0.sum() == 0.0){
            continue;
        }
        for(double dt: {1e-6, 1.0, 3600.0, 3.15e7, 3.15e10}){
            const Eigen::VectorXd full = cram(mat, n0, dt);
            const Eigen::VectorXd reduced = solver::solve_reachable<double>(mat, n0, dt, cram);
            INFO( "trial " << trial << ", size " << size << ", dt " << dt );
            REQUIRE_THAT( full.sum(), WithinRel(n0.sum(), 1e-13) );
            REQUIRE_THAT( reduced.sum(), WithinRel(n0.sum(), 1e-13) );
            // a value can only be below zero by rounding noise
            REQUIRE( full.minCoeff() >= -1e-15*n0.sum() );
            REQUIRE( reduced.minCoeff() >= -1e-15*n0.sum() );
        }
    }
}

TEST_CASE( "CRAM48: many short steps give the same as one long step" ) {
    std::mt19937 rng(42);
    for(int trial=0; trial<20; ++trial){
        const int size = 10 + trial;
        // half-lives from 1 ms to 30 years, so the steps are not all much longer than the half-lives
        const Eigen::SparseMatrix<double> mat = random_decay_matrix(size, rng, -3.0, 9.0);
        const Eigen::VectorXd n0 = random_atoms(size, rng);
        for(double dt: {10.0, 1e4, 1e7}){
            const Eigen::VectorXd one = solver::solve_reachable<double>(mat, n0, dt, cram);
            for(int n_steps: {2, 7, 50}){
                Eigen::VectorXd many = n0;
                for(int s=0; s<n_steps; ++s){
                    many = solver::solve_reachable<double>(mat, many, dt/n_steps, cram);
                }
                for(int i=0; i<size; ++i){
                    INFO( "trial " << trial << ", dt " << dt << ", steps " << n_steps << ", nuclide " << i
                          << ", one step " << one(i) << ", many steps " << many(i) );
                    REQUIRE( close(one(i), many(i), 1e-9, 1e-14*n0.sum()) );
                }
            }
        }
    }
}

TEST_CASE( "CRAM48: nuclides with no atoms and other chains do not change the result" ) {
    std::mt19937 rng(7);
    for(int trial=0; trial<20; ++trial){
        const int size = 8 + trial;
        const Eigen::SparseMatrix<double> small = random_decay_matrix(size, rng);
        const Eigen::VectorXd n0_small = random_atoms(size, rng);

        // the same problem in a larger matrix: other chains with no atoms, some of them very short-lived
        // and some very long-lived, and the nuclides of the small problem are spread out in the large one
        const int extra = 3*size;
        const int big_size = size + extra;
        std::vector<int> position(big_size);
        std::iota(position.begin(), position.end(), 0);
        std::shuffle(position.begin(), position.end(), rng);
        Triplets triplets;
        for(int col=0; col<small.outerSize(); ++col){
            for(Eigen::SparseMatrix<double>::InnerIterator it(small, col); it; ++it){
                triplets.emplace_back(position[it.row()], position[it.col()], it.value());
            }
        }
        const Eigen::SparseMatrix<double> other = random_decay_matrix(extra, rng, -17.0, 17.0);
        for(int col=0; col<other.outerSize(); ++col){
            for(Eigen::SparseMatrix<double>::InnerIterator it(other, col); it; ++it){
                triplets.emplace_back(position[size + it.row()], position[size + it.col()], it.value());
            }
        }
        const Eigen::SparseMatrix<double> big = make_matrix(big_size, triplets);
        Eigen::VectorXd n0_big = Eigen::VectorXd::Zero(big_size);
        for(int i=0; i<size; ++i){
            n0_big(position[i]) = n0_small(i);
        }

        for(double dt: {1e-3, 30.0, 1e6}){
            const Eigen::VectorXd expected = cram(small, n0_small, dt);
            const Eigen::VectorXd full = cram(big, n0_big, dt);
            const Eigen::VectorXd reduced = solver::solve_reachable<double>(big, n0_big, dt, cram);
            for(int i=0; i<size; ++i){
                INFO( "trial " << trial << ", dt " << dt << ", nuclide " << i );
                REQUIRE( close(full(position[i]), expected(i), 1e-10, 1e-15*n0_small.sum()) );
                REQUIRE( close(reduced(position[i]), expected(i), 1e-10, 1e-15*n0_small.sum()) );
            }
            for(int i=size; i<big_size; ++i){
                // nothing can get to the other chains
                REQUIRE( reduced(position[i]) == 0.0 );
                REQUIRE( std::abs(full(position[i])) <= 1e-15*n0_small.sum() );
            }
        }
    }
}

TEST_CASE( "CRAM48: a small change of dt gives a small change of the result" ) {
    std::mt19937 rng(1234);
    for(int trial=0; trial<20; ++trial){
        const int size = 10 + trial;
        const Eigen::SparseMatrix<double> mat = random_decay_matrix(size, rng, -17.0, 12.0);
        const Eigen::VectorXd n0 = random_atoms(size, rng);
        for(double dt: {0.1, 16.02, 3600.0, 1e8}){
            const Eigen::VectorXd a = solver::solve_reachable<double>(mat, n0, dt, cram);
            const Eigen::VectorXd b = solver::solve_reachable<double>(mat, n0, dt*(1.0 + 1e-12), cram);
            for(int i=0; i<size; ++i){
                INFO( "trial " << trial << ", dt " << dt << ", nuclide " << i << ": " << a(i) << " " << b(i) );
                REQUIRE( close(a(i), b(i), 1e-9, 1e-15*n0.sum()) );
            }
        }
    }
}

TEST_CASE( "CRAM48 with no time or no decay gives the atoms at the start" ) {
    std::mt19937 rng(99);
    const Eigen::VectorXd n0 = random_atoms(12, rng);
    // dt = 0
    const Eigen::VectorXd no_time = cram(random_decay_matrix(12, rng), n0, 0.0);
    // all the nuclides stable, for a very long time
    const Eigen::VectorXd stable = cram(Eigen::SparseMatrix<double>(12, 12), n0, 1e15);
    for(int i=0; i<12; ++i){
        REQUIRE( close(no_time(i), n0(i), 1e-13, 1e-30) );
        REQUIRE( close(stable(i), n0(i), 1e-13, 1e-30) );
    }
}
