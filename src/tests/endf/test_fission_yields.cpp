#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nis/nuclides.hpp"
#include "nis/endf/fission_yields.hpp"
#include "nis/endf/reaction.hpp"

using namespace nis;
using namespace nis::endf;

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

namespace {

    // MF=1 MT=451 (LISO=1, so the parent is Am242m) and MF=8 MT=454 with 2 energies:
    // Xe135, Cs137 and Ba137m at 0.0253 eV and 14 MeV
    const std::string test_yield_tape = R"(
 test fission yield tape                                             1 0  0    0
 9.524200+4 2.400000+2          0          0          0          09547 1451    1
 0.000000+0 0.000000+0          0          1          0          69547 1451    2
 1.000000+0 0.000000+0          0          0          0          09547 1451    3
 0.000000+0 0.000000+0          0          0          0          19547 1451    4
 test                                                             9547 1451    5
                                1        451          5          09547 1451    6
                                                                  9547 1  099999
                                                                  9547 0  0    0
 9.524200+4 2.400000+2          2          0          0          09547 8454    1
 2.530000-2 0.000000+0          0          0         12          39547 8454    2
 5.413500+4 0.000000+0 1.000000-2 1.000000-3 5.513700+4 0.000000+09547 8454    3
 6.000000-2 1.000000-3 5.613700+4 1.000000+0 2.000000-3 0.000000+09547 8454    4
 1.400000+7 0.000000+0          2          0         12          39547 8454    5
 5.413500+4 0.000000+0 3.000000-2 1.000000-3 5.513700+4 0.000000+09547 8454    6
 5.000000-2 1.000000-3 5.613700+4 1.000000+0 4.000000-3 0.000000+09547 8454    7
                                                                  9547 8  099999
                                                                  9547 0  0    0
                                                                     0 0  0    0
                                                                    -1 0  0    0
)";

    std::string tape_file(){
        const auto path = std::filesystem::temp_directory_path() / "nis_test_fission_yields.endf";
        std::ofstream out(path);
        // the raw string starts with a new line
        out << test_yield_tape.substr(1);
        return path.string();
    }

    FissionYields u235_yields(){
        FissionYields yields;
        yields.parent_zai = 922350;
        yields.energies = {0.0253, 4e5, 1.4e7};
        yields.yields = {{{541350, 0.06}, {551370, 0.06}},
                         {{541350, 0.05}, {551370, 0.07}},
                         {{541350, 0.02}, {551370, 0.04}, {9999990, 0.5}}};
        return yields;
    }
}

TEST_CASE( "read_fission_yields_from_file" ) {
    const auto filename = tape_file();
    const FissionYieldData data = read_fission_yields_from_file(filename);
    std::filesystem::remove(filename);

    REQUIRE( data.size() == 1 );
    // the isomer from LISO in MF=1 MT=451
    REQUIRE( data.count(952421) == 1 );
    const FissionYields& am = data.at(952421);
    REQUIRE( am.parent_zai == 952421 );
    REQUIRE( am.energies == std::vector<double>{0.0253, 1.4e7} );
    REQUIRE( am.yields.size() == 2 );
    REQUIRE( am.yields[0].size() == 3 );
    REQUIRE( am.yields[0].at(541350) == 0.01 );
    REQUIRE( am.yields[0].at(551370) == 0.06 );
    // the product isomer from FPS
    REQUIRE( am.yields[0].at(561371) == 0.002 );
    REQUIRE( am.yields[1].at(541350) == 0.03 );
}

TEST_CASE( "read_fission_yields_from_file with a missing file" ) {
    REQUIRE_THROWS( read_fission_yields_from_file("no_such_file.nfy") );
}

TEST_CASE( "fission yield energy_weights" ) {
    const std::vector<double> energies = {0.0253, 4e5, 1.4e7};
    // the nearest set outside the energies
    REQUIRE( energy_weights(energies, 0.0) == std::vector<double>{1.0, 0.0, 0.0} );
    REQUIRE( energy_weights(energies, 0.0253) == std::vector<double>{1.0, 0.0, 0.0} );
    REQUIRE( energy_weights(energies, 2e7) == std::vector<double>{0.0, 0.0, 1.0} );
    // linear in energy between 2 energies
    const auto w = energy_weights(energies, 7.2e6);
    REQUIRE_THAT( w[0], WithinAbs(0.0, 1e-15) );
    REQUIRE_THAT( w[1], WithinRel(0.5, 1e-12) );
    REQUIRE_THAT( w[2], WithinRel(0.5, 1e-12) );
    // one energy: always that set
    REQUIRE( energy_weights({0.0}, 1e6) == std::vector<double>{1.0} );
    REQUIRE_THROWS_AS( energy_weights({}, 1.0), std::invalid_argument );
}

TEST_CASE( "fission yields_at_energy" ) {
    const auto yields = u235_yields();
    // 1 eV is a little above 0.0253 eV, so a very small part of the 400 keV set
    REQUIRE_THAT( yields_at_energy(yields, 1.0).at(541350), WithinRel(0.06, 1e-5) );
    const auto mid = yields_at_energy(yields, 7.2e6);
    REQUIRE_THAT( mid.at(541350), WithinRel(0.035, 1e-12) );
    REQUIRE_THAT( mid.at(551370), WithinRel(0.055, 1e-12) );
    // a product only at 14 MeV
    REQUIRE_THAT( mid.at(9999990), WithinRel(0.25, 1e-12) );
    REQUIRE( yields_at_energy(yields, 1.0).count(9999990) == 0 );
}

TEST_CASE( "weighted fission yields" ) {
    const auto yields = u235_yields();
    // 3 fissions per s in a thermal group, 1 in a 14 MeV group
    const std::vector<double> bounds = {0.02, 0.03, 1.39e7, 1.41e7};
    const auto w = weighted_yields(yields, bounds, {3.0, 0.0, 1.0});
    // the middle of each group is near 0.0253 eV and 14 MeV, so only a very small part of the 400 keV set
    REQUIRE_THAT( w.at(541350), WithinRel((3.0*0.06 + 1.0*0.02)/4.0, 1e-5) );
    REQUIRE_THAT( w.at(551370), WithinRel((3.0*0.06 + 1.0*0.04)/4.0, 1e-5) );

    // a group from 0 uses the middle in energy, so all thermal
    REQUIRE_THAT( weighted_yields(yields, {0.0, 0.04}, {1.0}).at(541350), WithinRel(0.06, 1e-12) );
    // no fissions: the first set
    REQUIRE_THAT( weighted_yields(yields, bounds, {0.0, 0.0, 0.0}).at(541350), WithinRel(0.06, 1e-12) );
    REQUIRE_THROWS_AS( weighted_yields(yields, bounds, {1.0}), std::invalid_argument );
}

TEST_CASE( "construct_reaction_matrix with fission yields" ) {
    // fission of U235: 500 b in one group from 1 to 2 eV
    const ReactionData reaction_data = {
        CrossSection{922350, FISSION_MT, {1.0, 2.0}, {500.0, 0.0}},
    };
    const core::Spectrum<double> spectrum({1.0, 2.0}, {1e10});
    const double rate = 1e10*500.0*BARN_TO_CM2;
    const int u235 = data::ZAI_MAP.at(922350);
    const int xe135 = data::ZAI_MAP.at(541350);
    const int cs137 = data::ZAI_MAP.at(551370);

    SECTION( "without yields, fission goes to the sink" ) {
        const auto mat = construct_reaction_matrix(reaction_data, spectrum);
        REQUIRE_THAT( mat.coeff(u235, u235), WithinRel(-rate, 1e-12) );
        REQUIRE_THAT( mat.coeff(u235, data::SINK_INDEX), WithinRel(rate, 1e-12) );
        REQUIRE( mat.coeff(u235, xe135) == 0.0 );
    }

    SECTION( "with yields, fission makes the fission products" ) {
        const FissionYieldData yields = {{922350, u235_yields()}};
        const auto mat = construct_reaction_matrix(reaction_data, spectrum, &yields);
        REQUIRE_THAT( mat.coeff(u235, u235), WithinRel(-rate, 1e-12) );
        // the group is at 1.4 eV, so the thermal yields (linear in energy to the 400 keV set)
        const auto y = yields_at_energy(u235_yields(), std::sqrt(2.0));
        REQUIRE_THAT( mat.coeff(u235, xe135), WithinRel(rate*y.at(541350), 1e-12) );
        REQUIRE_THAT( mat.coeff(u235, cs137), WithinRel(rate*y.at(551370), 1e-12) );
        // nothing to the sink here, the 14 MeV only product has no weight
        REQUIRE( mat.coeff(u235, data::SINK_INDEX) == 0.0 );
    }

    SECTION( "a product we do not track goes to the sink" ) {
        const FissionYieldData yields = {{922350, u235_yields()}};
        const core::Spectrum<double> fast({1.0e7, 2.0e7}, {1e10});
        const ReactionData fast_data = {
            CrossSection{922350, FISSION_MT, {1.0e7, 2.0e7}, {2.0, 0.0}},
        };
        const auto mat = construct_reaction_matrix(fast_data, fast, &yields);
        const double fast_rate = 1e10*2.0*BARN_TO_CM2;
        REQUIRE_THAT( mat.coeff(u235, data::SINK_INDEX), WithinRel(fast_rate*0.5, 1e-12) );
        REQUIRE_THAT( mat.coeff(u235, xe135), WithinRel(fast_rate*0.02, 1e-12) );
    }

    SECTION( "a nuclide with no yields still goes to the sink" ) {
        const FissionYieldData yields = {{942390, u235_yields()}};
        const auto mat = construct_reaction_matrix(reaction_data, spectrum, &yields);
        REQUIRE_THAT( mat.coeff(u235, data::SINK_INDEX), WithinRel(rate, 1e-12) );
    }
}
