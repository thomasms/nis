#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <ENDFtk.hpp>

#include "nis/cram.hpp"
#include "nis/decayheat.hpp"
#include "nis/nuclides.hpp"
#include "nis/physics.hpp"
#include "nis/endf/decay.hpp"

using namespace nis;

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

namespace {

    // A small decay tape: Co60 (beta-), stable Ni60 and Po210 (alpha).
    // Values are close to ENDF/B-VIII.0, but this is test data only.
    const std::string test_tape_string =
            "nis test decay tape                                                  1 0  0     \n"
            " 2.706000+4 5.941900+1          0          0          0          02731 8457     \n"
            " 1.663415+8 0.000000+0          0          0          6          02731 8457     \n"
            " 9.577000+4 0.000000+0 2.503600+6 0.000000+0 0.000000+0 0.000000+02731 8457     \n"
            " 0.000000+0 1.000000+0          0          0          6          12731 8457     \n"
            " 1.000000+0 0.000000+0 2.823900+6 0.000000+0 1.000000+0 0.000000+02731 8457     \n"
            "                                                                  2731 8  0     \n"
            "                                                                  2731 0  0     \n"
            "                                                                     0 0  0     \n"
            " 2.806000+4 5.941600+1          0          0          1          02834 8457     \n"
            " 0.000000+0 0.000000+0          0          0          6          02834 8457     \n"
            " 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02834 8457     \n"
            " 0.000000+0 1.000000+0          0          0          6          02834 8457     \n"
            " 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02834 8457     \n"
            "                                                                  2834 8  0     \n"
            "                                                                  2834 0  0     \n"
            "                                                                     0 0  0     \n"
            " 8.421000+4 2.081840+2          0          0          0          08434 8457     \n"
            " 1.195569+7 0.000000+0          0          0          6          08434 8457     \n"
            " 0.000000+0 0.000000+0 8.500000+0 0.000000+0 5.407500+6 0.000000+08434 8457     \n"
            " 0.000000+0 1.000000+0          0          0          6          18434 8457     \n"
            " 4.000000+0 0.000000+0 5.407500+6 0.000000+0 1.000000+0 0.000000+08434 8457     \n"
            "                                                                  8434 8  0     \n"
            "                                                                  8434 0  0     \n"
            "                                                                     0 0  0     \n"
            "                                                                    -1 0  0     \n"
        ;

    constexpr double CO60_HALF_LIFE = 1.663415e8;
    constexpr double CO60_ELP = 9.577e4;
    constexpr double CO60_EEM = 2.5036e6;
    constexpr double PO210_HALF_LIFE = 1.195569e7;
    constexpr double PO210_EEM = 8.5;
    constexpr double PO210_EHP = 5.4075e6;

    constexpr double AVOGADRO = 6.02214076e23;
    // g/mol
    constexpr double CO60_MOLAR_MASS = 59.9338;

    njoy::ENDFtk::tree::Tape test_tape(){
        return njoy::ENDFtk::tree::Tape(test_tape_string);
    }

    int index(int zai){
        return data::ZAI_MAP.at(zai);
    }
}

TEST_CASE( "read_decay_energy_data_from_tape" ) {
    const auto energy_data = endf::read_decay_energy_data_from_tape(test_tape());

    const int co60 = index(270600);
    REQUIRE_THAT( energy_data.decay_constants(co60), WithinRel(physics::LN_2/CO60_HALF_LIFE, 1e-12) );
    REQUIRE_THAT( energy_data.light_particle(co60), WithinRel(CO60_ELP, 1e-12) );
    REQUIRE_THAT( energy_data.electromagnetic(co60), WithinRel(CO60_EEM, 1e-12) );
    REQUIRE( energy_data.heavy_particle(co60) == 0.0 );

    const int po210 = index(842100);
    REQUIRE_THAT( energy_data.decay_constants(po210), WithinRel(physics::LN_2/PO210_HALF_LIFE, 1e-12) );
    REQUIRE_THAT( energy_data.heavy_particle(po210), WithinRel(PO210_EHP, 1e-12) );

    // stable nuclides have no decay energy
    const int ni60 = index(280600);
    REQUIRE( energy_data.decay_constants(ni60) == 0.0 );
    REQUIRE( energy_data.electromagnetic(ni60) == 0.0 );

    // only 2 radioactive nuclides on the tape
    REQUIRE( (energy_data.decay_constants.array() > 0.0).count() == 2 );
}

TEST_CASE( "decay matrix from the test tape" ) {
    const Eigen::SparseMatrix<double> mat = endf::construct_decay_matrix_from_tape(test_tape());

    const double co60_lambda = physics::LN_2/CO60_HALF_LIFE;
    REQUIRE_THAT( mat.coeff(index(270600), index(270600)), WithinRel(-co60_lambda, 1e-12) );
    REQUIRE_THAT( mat.coeff(index(270600), index(280600)), WithinRel(co60_lambda, 1e-12) );

    const double po210_lambda = physics::LN_2/PO210_HALF_LIFE;
    REQUIRE_THAT( mat.coeff(index(842100), index(822060)), WithinRel(po210_lambda, 1e-12) );
    REQUIRE_THAT( mat.coeff(index(842100), index(zai::ALPHA)), WithinRel(po210_lambda, 1e-12) );
}

TEST_CASE( "decay heat of 1 g of Co60" ) {
    const auto energy_data = endf::read_decay_energy_data_from_tape(test_tape());

    Eigen::VectorXd inventory = Eigen::VectorXd::Zero(data::ZAI_SIZE);
    const double n = AVOGADRO/CO60_MOLAR_MASS;
    inventory(index(270600)) = n;

    const double lambda = physics::LN_2/CO60_HALF_LIFE;
    const double activity = core::compute_activity<double>(energy_data, inventory);
    REQUIRE_THAT( activity, WithinRel(lambda*n, 1e-12) );
    // the known specific activity of Co60 is 1131 Ci/g = 4.185e13 Bq/g
    REQUIRE_THAT( activity, WithinRel(4.185e13, 1e-3) );

    const auto heat = core::compute_decay_heat<double>(energy_data, inventory);
    REQUIRE_THAT( heat.light_particle, WithinRel(lambda*n*CO60_ELP*physics::EV_TO_J, 1e-12) );
    REQUIRE_THAT( heat.electromagnetic, WithinRel(lambda*n*CO60_EEM*physics::EV_TO_J, 1e-12) );
    REQUIRE( heat.heavy_particle == 0.0 );
    // about 17.4 W per gram
    REQUIRE_THAT( heat.total(), WithinRel(17.4, 1e-2) );
}

TEST_CASE( "decay heat by component for a mix" ) {
    const auto energy_data = endf::read_decay_energy_data_from_tape(test_tape());

    Eigen::VectorXd inventory = Eigen::VectorXd::Zero(data::ZAI_SIZE);
    inventory(index(270600)) = 1e20;
    inventory(index(842100)) = 1e18;
    // stable, so no heat
    inventory(index(280600)) = 1e24;

    const double co60_a = 1e20*physics::LN_2/CO60_HALF_LIFE;
    const double po210_a = 1e18*physics::LN_2/PO210_HALF_LIFE;

    const Eigen::VectorXd activities = core::compute_activities<double>(energy_data, inventory);
    REQUIRE( activities(index(280600)) == 0.0 );
    REQUIRE_THAT( activities.sum(), WithinRel(co60_a + po210_a, 1e-12) );

    const auto heat = core::compute_decay_heat<double>(energy_data, inventory);
    REQUIRE_THAT( heat.light_particle, WithinRel(co60_a*CO60_ELP*physics::EV_TO_J, 1e-12) );
    REQUIRE_THAT( heat.electromagnetic, WithinRel((co60_a*CO60_EEM + po210_a*PO210_EEM)*physics::EV_TO_J, 1e-12) );
    REQUIRE_THAT( heat.heavy_particle, WithinRel(po210_a*PO210_EHP*physics::EV_TO_J, 1e-12) );
}

TEST_CASE( "decay heat halves after one half life" ) {
    const auto tape = test_tape();
    const Eigen::SparseMatrix<double> mat = endf::construct_decay_matrix_from_tape(tape);
    const auto energy_data = endf::read_decay_energy_data_from_tape(tape);

    Eigen::VectorXd n0 = Eigen::VectorXd::Zero(data::ZAI_SIZE);
    n0(index(270600)) = 1e20;
    const Eigen::VectorXd nt = core::cram48::solve<double>(mat, n0, CO60_HALF_LIFE);

    REQUIRE_THAT( nt(index(280600)), WithinRel(0.5e20, 1e-10) );
    const double heat0 = core::compute_decay_heat<double>(energy_data, n0).total();
    const double heat1 = core::compute_decay_heat<double>(energy_data, nt).total();
    REQUIRE_THAT( heat1, WithinRel(0.5*heat0, 1e-10) );
}

TEST_CASE( "no inventory gives no heat" ) {
    const core::DecayEnergyData<double> energy_data;
    const Eigen::VectorXd inventory = Eigen::VectorXd::Zero(data::ZAI_SIZE);
    REQUIRE( core::compute_activity<double>(energy_data, inventory) == 0.0 );
    REQUIRE( core::compute_decay_heat<double>(energy_data, inventory).total() == 0.0 );
}
