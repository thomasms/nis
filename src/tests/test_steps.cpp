#include <cmath>
#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Sparse>

#include "nis/cram.hpp"
#include "nis/nuclides.hpp"
#include "nis/solver.hpp"
#include "nis/zai.hpp"

using namespace nis;

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

namespace {

    const solver::solver_func<double> cram = &core::cram48::solve<double>;

    int index(const std::string& name){
        return data::ZAI_MAP.at(zai::get_from_name(name));
    }

    Eigen::SparseMatrix<double> matrix(const std::vector<Eigen::Triplet<double>>& triplets){
        const int size = data::ZAI_SIZE;
        Eigen::SparseMatrix<double> mat(size, size);
        mat.setFromTriplets(triplets.begin(), triplets.end());
        return mat;
    }

    // Co60 decays to Ni60 (stable), half life 5.27 y
    const double co60_lambda = std::log(2.0) / (5.2714 * 365.25 * 86400.0);

    Eigen::SparseMatrix<double> decay_matrix(){
        return matrix({{index("Co60"), index("Co60"), -co60_lambda},
                       {index("Co60"), index("Ni60"), co60_lambda}});
    }

    // Co59 (n,g) Co60, a reaction rate of 1e-9 /s for the spectrum as it is given
    const double co59_rate = 1e-9;

    Eigen::SparseMatrix<double> reaction_matrix(){
        return matrix({{index("Co59"), index("Co59"), -co59_rate},
                       {index("Co59"), index("Co60"), co59_rate}});
    }

    core::BasicModel<double> co59_model(const std::vector<core::TimeStep<double>>& steps){
        core::BasicModel<double> model;
        model.nuclides = {{"Co59", 1e20}};
        model.steps = steps;
        model.dt = steps.empty() ? 0.0 : steps.front().dt;
        return model;
    }

    double atoms(const core::ModelResult<double>& result, const std::string& name){
        const auto it = result.nuclides.find(name);
        return it == result.nuclides.end() ? 0.0 : it->second;
    }
}

TEST_CASE( "time_steps with no dt line is one step of dt" ) {
    core::BasicModel<double> model;
    model.dt = 10.0;
    const auto steps = solver::time_steps(model);
    REQUIRE( steps.size() == 1 );
    REQUIRE( steps[0].dt == 10.0 );
    REQUIRE_FALSE( steps[0].flux.has_value() );

    model.steps = {{5.0, 1e14}, {7.0, 0.0}};
    REQUIRE( solver::time_steps(model).size() == 2 );
}

TEST_CASE( "flux_scale" ) {
    // no flux: the spectrum as it is given
    REQUIRE( solver::flux_scale<double>({1.0, std::nullopt}, 2e10) == 1.0 );
    REQUIRE( solver::flux_scale<double>({1.0, std::nullopt}, 0.0) == 1.0 );
    // a flux of 0 is a cooling step
    REQUIRE( solver::flux_scale<double>({1.0, 0.0}, 2e10) == 0.0 );
    // the step flux is in n/m^2/s, the spectrum in n/cm^2/s
    REQUIRE_THAT( solver::flux_scale<double>({1.0, 1e14}, 2e10), WithinRel(1e14 * 1e-4 / 2e10, 1e-15) );
    // a flux with no spectrum
    REQUIRE_THROWS_AS( solver::flux_scale<double>({1.0, 1e14}, 0.0), std::invalid_argument );
}

TEST_CASE( "compute_steps with decay only" ) {
    const auto decay = decay_matrix();
    core::BasicModel<double> model;
    model.nuclides = {{"Co60", 1e20}};
    const double half_life = std::log(2.0) / co60_lambda;
    model.steps = {{half_life, std::nullopt}, {half_life, 0.0}};
    model.dt = half_life;

    const auto results = solver::compute_steps<double>(model, decay, Eigen::SparseMatrix<double>(), 0.0, cram);
    REQUIRE( results.size() == 2 );
    REQUIRE_THAT( atoms(results[0], "Co60"), WithinRel(0.5e20, 1e-10) );
    REQUIRE_THAT( atoms(results[1], "Co60"), WithinRel(0.25e20, 1e-10) );
    REQUIRE_THAT( atoms(results[1], "Ni60"), WithinRel(0.75e20, 1e-10) );
}

TEST_CASE( "compute_steps irradiation then cooling" ) {
    const auto decay = decay_matrix();
    const auto reaction = reaction_matrix();
    const double spectrum_total = 1e10;   // n/cm^2/s
    const double irradiation = 3600.0;
    const double cooling = 86400.0;

    // a flux of 2 times the spectrum total, in n/m^2/s
    const double flux = 2.0 * spectrum_total * 1e4;
    const auto results = solver::compute_steps<double>(
        co59_model({{irradiation, flux}, {cooling, 0.0}}), decay, reaction, spectrum_total, cram);
    REQUIRE( results.size() == 2 );

    // after the irradiation: Co59 burns with 2 times the rate, analytic solution of Co59 -> Co60 -> Ni60
    const double rate = 2.0 * co59_rate;
    const double co59 = 1e20 * std::exp(-rate * irradiation);
    const double co60 = 1e20 * rate / (co60_lambda - rate)
                        * (std::exp(-rate * irradiation) - std::exp(-co60_lambda * irradiation));
    REQUIRE_THAT( atoms(results[0], "Co59"), WithinRel(co59, 1e-10) );
    REQUIRE_THAT( atoms(results[0], "Co60"), WithinRel(co60, 1e-8) );

    // the cooling step is decay only: Co59 does not change, Co60 decays
    REQUIRE_THAT( atoms(results[1], "Co59"), WithinRel(atoms(results[0], "Co59"), 1e-12) );
    REQUIRE_THAT( atoms(results[1], "Co60"), WithinRel(co60 * std::exp(-co60_lambda * cooling), 1e-8) );

    // the cooling step is the same as a step with no reactions
    core::BasicModel<double> after;
    after.nuclides = results[0].nuclides;
    after.dt = cooling;
    const auto decay_only = solver::compute_inventory<double>(after, decay, cram);
    REQUIRE( atoms(results[1], "Co60") == atoms(decay_only, "Co60") );
}

TEST_CASE( "compute_steps with no flux uses the spectrum as it is given" ) {
    const auto decay = decay_matrix();
    const auto reaction = reaction_matrix();
    const double spectrum_total = 1e10;
    const auto given = solver::compute_steps<double>(
        co59_model({{3600.0, std::nullopt}}), decay, reaction, spectrum_total, cram);
    // the same as a flux equal to the spectrum total
    const auto same = solver::compute_steps<double>(
        co59_model({{3600.0, spectrum_total * 1e4}}), decay, reaction, spectrum_total, cram);
    REQUIRE_THAT( atoms(given[0], "Co60"), WithinRel(atoms(same[0], "Co60"), 1e-12) );
    REQUIRE_THAT( atoms(given[0], "Co59"), WithinRel(1e20 * std::exp(-co59_rate * 3600.0), 1e-10) );
}

TEST_CASE( "compute_steps: two half steps are the same as one step" ) {
    const auto decay = decay_matrix();
    const auto reaction = reaction_matrix();
    const double flux = 1e14;
    const auto one = solver::compute_steps<double>(
        co59_model({{7200.0, flux}}), decay, reaction, 1e10, cram);
    const auto two = solver::compute_steps<double>(
        co59_model({{3600.0, flux}, {3600.0, flux}}), decay, reaction, 1e10, cram);
    REQUIRE( two.size() == 2 );
    REQUIRE_THAT( atoms(two[1], "Co60"), WithinRel(atoms(one[0], "Co60"), 1e-10) );
    REQUIRE_THAT( atoms(two[1], "Co59"), WithinRel(atoms(one[0], "Co59"), 1e-12) );
}

TEST_CASE( "compute_steps bad input" ) {
    const auto decay = decay_matrix();
    // a flux, but no spectrum
    REQUIRE_THROWS_AS( solver::compute_steps<double>(co59_model({{1.0, 1e14}}), decay, reaction_matrix(), 0.0, cram),
                       std::invalid_argument );
    // matrices of different size
    Eigen::SparseMatrix<double> small(3, 3);
    small.insert(0, 0) = -1.0;
    REQUIRE_THROWS_AS( solver::compute_steps<double>(co59_model({{1.0, 1e14}}), decay, small, 1e10, cram),
                       std::invalid_argument );
}
