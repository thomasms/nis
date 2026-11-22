#include <algorithm>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// only the data layer, so these tests do not need ENDF
#include "nis/nuclides.hpp"
#include "nis/spectrum.hpp"
#include "nis/data/fission_yields.hpp"
#include "nis/data/reaction.hpp"

using namespace nis;

using Catch::Matchers::WithinRel;

namespace {

    constexpr int FE56 = 260560;
    constexpr int FE57 = 260570;
    constexpr int MN56 = 250560;
    constexpr int U235 = 922350;

    int index(int zai){
        return data::ZAI_MAP.at(zai);
    }

    data::ReactionData fe56_data(){
        return {
            {FE56, 102, {1.0, 10.0, 100.0}, {2.0, 3.0, 3.0}},
            {FE56, 103, {10.0, 100.0}, {0.5, 0.5}},
        };
    }

    const core::Spectrum<double> FLUX({1.0, 10.0, 100.0}, {1e12, 4e12});
}

TEST_CASE( "change a cross section, then make the reaction matrix again" ) {
    data::ReactionData reaction_data = fe56_data();
    const Eigen::SparseMatrix<double> before = data::construct_reaction_matrix(reaction_data, FLUX);

    // 2 times the (n,g) cross section
    for(auto& v: reaction_data[0].values){
        v *= 2.0;
    }
    const Eigen::SparseMatrix<double> after = data::construct_reaction_matrix(reaction_data, FLUX);

    REQUIRE_THAT( after.coeff(index(FE56), index(FE57)), WithinRel(2.0*before.coeff(index(FE56), index(FE57)), 1e-14) );
    REQUIRE_THAT( after.coeff(index(FE56), index(MN56)), WithinRel(before.coeff(index(FE56), index(MN56)), 1e-14) );
}

TEST_CASE( "add and remove a reaction" ) {
    data::ReactionData reaction_data = fe56_data();

    // add (n,a) to Cr53
    reaction_data.push_back({FE56, 107, {10.0, 100.0}, {0.1, 0.1}});
    Eigen::SparseMatrix<double> mat = data::construct_reaction_matrix(reaction_data, FLUX);
    const double na = 0.1*4e12*data::BARN_TO_CM2;
    REQUIRE_THAT( mat.coeff(index(FE56), index(240530)), WithinRel(na, 1e-12) );
    REQUIRE_THAT( mat.coeff(index(FE56), index(20040)), WithinRel(na, 1e-12) );

    // remove (n,g)
    reaction_data.erase(std::remove_if(reaction_data.begin(), reaction_data.end(),
                                       [](const data::CrossSection& xs){ return xs.mt_number == 102; }),
                        reaction_data.end());
    mat = data::construct_reaction_matrix(reaction_data, FLUX);
    REQUIRE( mat.coeff(index(FE56), index(FE57)) == 0.0 );
    REQUIRE( mat.coeff(index(FE56), index(MN56)) > 0.0 );
}

TEST_CASE( "fission with yields made by hand" ) {
    const data::ReactionData reaction_data = {
        {U235, data::FISSION_MT, {1.0, 10.0, 100.0}, {100.0, 100.0, 100.0}},
    };
    data::FissionYieldData fission_yields;
    fission_yields[U235] = data::FissionYields{U235, {0.0253}, {{{561400, 1.2}, {421000, 0.8}}}};

    const double rate = 100.0*5e12*data::BARN_TO_CM2;

    // without yields, fission goes to the sink
    const Eigen::SparseMatrix<double> to_sink = data::construct_reaction_matrix(reaction_data, FLUX);
    REQUIRE_THAT( to_sink.coeff(index(U235), data::SINK_INDEX), WithinRel(rate, 1e-12) );

    const Eigen::SparseMatrix<double> mat = data::construct_reaction_matrix(reaction_data, FLUX, &fission_yields);
    REQUIRE_THAT( mat.coeff(index(U235), index(U235)), WithinRel(-rate, 1e-12) );
    REQUIRE_THAT( mat.coeff(index(U235), index(561400)), WithinRel(1.2*rate, 1e-12) );
    REQUIRE_THAT( mat.coeff(index(U235), index(421000)), WithinRel(0.8*rate, 1e-12) );
    REQUIRE( mat.coeff(index(U235), data::SINK_INDEX) == 0.0 );
}
