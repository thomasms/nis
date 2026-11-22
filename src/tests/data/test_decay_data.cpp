#include <cmath>
#include <limits>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// only the data layer, so these tests do not need ENDF
#include "nis/nuclides.hpp"
#include "nis/physics.hpp"
#include "nis/zai.hpp"
#include "nis/data/decay.hpp"
#include "nis/data/fission_yields.hpp"

using namespace nis;

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

namespace {

    constexpr int CO60 = 270600;
    constexpr int NI60 = 280600;
    constexpr int PO210 = 842100;
    constexpr int PB206 = 822060;
    constexpr int CF252 = 982520;
    constexpr int CM248 = 962480;

    int index(int zai){
        return data::ZAI_MAP.at(zai);
    }

    data::DecayNuclide beta_minus(double half_life){
        data::DecayNuclide nuclide;
        nuclide.half_life = half_life;
        nuclide.modes = {{1.0, 1.0, 0}};
        return nuclide;
    }
}

TEST_CASE( "DecayNuclide is stable by default" ) {
    const data::DecayNuclide nuclide;
    REQUIRE( nuclide.is_stable() );
    REQUIRE( nuclide.decay_constant() == 0.0 );
    REQUIRE( nuclide.modes.empty() );
}

TEST_CASE( "DecayNuclide decay constant and half-life" ) {
    data::DecayNuclide nuclide = beta_minus(100.0);
    REQUIRE( !nuclide.is_stable() );
    REQUIRE_THAT( nuclide.decay_constant(), WithinRel(physics::LN_2/100.0, 1e-15) );

    nuclide.set_decay_constant(0.01);
    REQUIRE_THAT( nuclide.half_life, WithinRel(physics::LN_2/0.01, 1e-15) );

    // 0, a negative half-life and infinity are all stable
    nuclide.set_decay_constant(0.0);
    REQUIRE( nuclide.is_stable() );
    nuclide.half_life = -1.0;
    REQUIRE( nuclide.is_stable() );
    nuclide.half_life = 0.0;
    REQUIRE( nuclide.is_stable() );
}

TEST_CASE( "DecayNuclide total branching ratio" ) {
    data::DecayNuclide nuclide;
    nuclide.half_life = 1.0;
    nuclide.modes = {{1.0, 0.25, 0}, {4.0, 0.7, 0}};
    REQUIRE_THAT( nuclide.total_branching_ratio(), WithinRel(0.95, 1e-15) );
}

TEST_CASE( "decay matrix from decay data made by hand" ) {
    data::DecayData decay_data;
    decay_data[CO60] = beta_minus(1.663415e8);
    // stable, no entries
    decay_data[NI60] = data::DecayNuclide();
    // alpha decay gives the daughter and a He4
    decay_data[PO210].half_life = 1.195569e7;
    decay_data[PO210].modes = {{4.0, 1.0, 0}};

    const Eigen::SparseMatrix<double> mat = data::construct_decay_matrix(decay_data);
    REQUIRE( mat.rows() == data::ZAI_SIZE );
    REQUIRE( mat.cols() == data::ZAI_SIZE );

    const double co60 = physics::LN_2/1.663415e8;
    REQUIRE_THAT( mat.coeff(index(CO60), index(CO60)), WithinRel(-co60, 1e-15) );
    REQUIRE_THAT( mat.coeff(index(CO60), index(NI60)), WithinRel(co60, 1e-15) );

    const double po210 = physics::LN_2/1.195569e7;
    REQUIRE_THAT( mat.coeff(index(PO210), index(PO210)), WithinRel(-po210, 1e-15) );
    REQUIRE_THAT( mat.coeff(index(PO210), index(PB206)), WithinRel(po210, 1e-15) );
    REQUIRE_THAT( mat.coeff(index(PO210), index(zai::ALPHA)), WithinRel(po210, 1e-15) );

    // 2 for Co60 and 3 for Po210
    REQUIRE( mat.nonZeros() == 5 );
}

TEST_CASE( "change a half-life and a branching ratio, then make the matrix again" ) {
    data::DecayData decay_data;
    decay_data[CO60] = beta_minus(1.663415e8);
    const Eigen::SparseMatrix<double> before = data::construct_decay_matrix(decay_data);

    // half the half-life, and 10% of the decays now go to the sink (unknown, RTYP=10)
    decay_data[CO60].half_life /= 2.0;
    decay_data[CO60].modes[0].branching_ratio = 0.9;
    decay_data[CO60].modes.push_back({10.0, 0.1, 0});
    const Eigen::SparseMatrix<double> after = data::construct_decay_matrix(decay_data);

    const double co60 = 2.0*physics::LN_2/1.663415e8;
    REQUIRE_THAT( after.coeff(index(CO60), index(CO60)), WithinRel(2.0*before.coeff(index(CO60), index(CO60)), 1e-15) );
    REQUIRE_THAT( after.coeff(index(CO60), index(NI60)), WithinRel(0.9*co60, 1e-15) );
    REQUIRE_THAT( after.coeff(index(CO60), data::SINK_INDEX), WithinRel(0.1*co60, 1e-15) );
}

TEST_CASE( "remove a nuclide from the decay data" ) {
    data::DecayData decay_data;
    decay_data[CO60] = beta_minus(1.663415e8);
    decay_data.erase(CO60);
    REQUIRE( data::construct_decay_matrix(decay_data).nonZeros() == 0 );
}

TEST_CASE( "decay data for a nuclide we do not track is not used" ) {
    data::DecayData decay_data;
    // Z=200 is not in ZAI_MAP
    decay_data[zai::make(200, 500, 0)] = beta_minus(1.0);
    REQUIRE( data::construct_decay_matrix(decay_data).nonZeros() == 0 );
    REQUIRE( data::construct_decay_energy_data(decay_data).decay_constants.sum() == 0.0 );
}

TEST_CASE( "a daughter we do not track goes to the sink" ) {
    data::DecayData decay_data;
    // beta- to the isomer Ni60m, it is not in ZAI_MAP
    REQUIRE( data::ZAI_MAP.count(280601) == 0 );
    decay_data[CO60].half_life = 1.0;
    decay_data[CO60].modes = {{1.0, 1.0, 1}};
    const Eigen::SparseMatrix<double> mat = data::construct_decay_matrix(decay_data);
    REQUIRE_THAT( mat.coeff(index(CO60), data::SINK_INDEX), WithinRel(physics::LN_2, 1e-15) );
    REQUIRE( mat.coeff(index(CO60), index(NI60)) == 0.0 );
}

TEST_CASE( "spontaneous fission with yields made by hand" ) {
    data::DecayData decay_data;
    decay_data[CF252].half_life = 8.3e7;
    decay_data[CF252].modes = {{4.0, 0.97, 0}, {6.0, 0.03, 0}};

    data::FissionYieldData sf_yields;
    sf_yields[CF252] = data::FissionYields{CF252, {0.0}, {{{561400, 1.5}, {421000, 0.5}}}};

    const double lambda = physics::LN_2/8.3e7;

    // without yields, spontaneous fission goes to the sink
    const Eigen::SparseMatrix<double> no_yields = data::construct_decay_matrix(decay_data);
    REQUIRE_THAT( no_yields.coeff(index(CF252), data::SINK_INDEX), WithinRel(0.03*lambda, 1e-15) );

    const Eigen::SparseMatrix<double> mat = data::construct_decay_matrix(decay_data, &sf_yields);
    REQUIRE( mat.coeff(index(CF252), data::SINK_INDEX) == 0.0 );
    REQUIRE_THAT( mat.coeff(index(CF252), index(561400)), WithinRel(0.03*lambda*1.5, 1e-15) );
    REQUIRE_THAT( mat.coeff(index(CF252), index(421000)), WithinRel(0.03*lambda*0.5, 1e-15) );
    REQUIRE_THAT( mat.coeff(index(CF252), index(CM248)), WithinRel(0.97*lambda, 1e-15) );
    REQUIRE_THAT( mat.coeff(index(CF252), index(CF252)), WithinRel(-lambda, 1e-15) );
}

TEST_CASE( "decay energy data from decay data made by hand" ) {
    data::DecayData decay_data;
    decay_data[CO60] = beta_minus(1.663415e8);
    decay_data[CO60].light_particle = 9.577e4;
    decay_data[CO60].electromagnetic = 2.5036e6;
    // a stable nuclide has no energy, even if it is given
    decay_data[NI60].electromagnetic = 1.0;

    const auto energy_data = data::construct_decay_energy_data(decay_data);
    REQUIRE( energy_data.decay_constants.size() == data::ZAI_SIZE );
    REQUIRE_THAT( energy_data.decay_constants(index(CO60)), WithinRel(physics::LN_2/1.663415e8, 1e-15) );
    REQUIRE( energy_data.light_particle(index(CO60)) == 9.577e4 );
    REQUIRE( energy_data.electromagnetic(index(CO60)) == 2.5036e6 );
    REQUIRE( energy_data.heavy_particle(index(CO60)) == 0.0 );
    REQUIRE( energy_data.decay_constants(index(NI60)) == 0.0 );
    REQUIRE( energy_data.electromagnetic(index(NI60)) == 0.0 );
}
