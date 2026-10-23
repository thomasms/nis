#include <filesystem>
#include <fstream>
#include <stdexcept>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <ENDFtk.hpp>

#include "nis/cram.hpp"
#include "nis/nuclides.hpp"
#include "nis/endf/reaction.hpp"
#include "nis/zai.hpp"

using namespace nis;
using namespace nis::endf;

using Catch::Matchers::WithinRel;

void check_products(int parent_zai, int mt_number, int expected_primary, const std::vector<int>& expected_secondaries){
    const ReactionProductT p = compute_reaction_products(parent_zai, mt_number);
    REQUIRE ( std::get<0>(p) == expected_primary );
    REQUIRE ( std::get<1>(p) == expected_secondaries );
}

TEST_CASE( "compute_reaction_products Fe56" ) {
    const int fe56 = 260560;
    check_products(fe56, 102, 260570, {});
    check_products(fe56, 16, 260550, {});
    check_products(fe56, 17, 260540, {});
    check_products(fe56, 103, 250560, {PROTON});
    check_products(fe56, 104, 250550, {DEUTERON});
    check_products(fe56, 105, 250540, {TRITON});
    check_products(fe56, 106, 240540, {HELION});
    check_products(fe56, 107, 240530, {zai::ALPHA});
    check_products(fe56, 22, 240520, {zai::ALPHA});
    check_products(fe56, 28, 250550, {PROTON});
}

TEST_CASE( "compute_reaction_products residual is in the ground state" ) {
    // Am242m (n,g) -> Am243
    check_products(952421, 102, 952430, {});
}

TEST_CASE( "compute_reaction_products fission goes to sink" ) {
    check_products(922350, FISSION_MT, -1, {});
}

TEST_CASE( "compute_reaction_products unphysical residual goes to sink" ) {
    // H1 (n,a) is not possible
    check_products(10010, 107, -1, {zai::ALPHA});
}

TEST_CASE( "compute_reaction_products unsupported MT" ) {
    // elastic and total do not change the nuclide
    REQUIRE ( !is_reaction_mt(1) );
    REQUIRE ( !is_reaction_mt(2) );
    REQUIRE ( !is_reaction_mt(51) );
    REQUIRE ( is_reaction_mt(102) );
    REQUIRE ( is_reaction_mt(FISSION_MT) );
    REQUIRE_THROWS_AS( compute_reaction_products(260560, 2), std::invalid_argument );
}

TEST_CASE( "collapse_cross_section same group structure" ) {
    const core::Spectrum<double> flux({1.0, 10.0, 100.0}, {2.0, 3.0});
    // last value is not used for INT=1
    const double rate = collapse_cross_section(flux, {1.0, 10.0, 100.0}, {5.0, 7.0, 7.0});
    REQUIRE_THAT( rate, WithinRel(31.0e-24, 1e-12) );
}

TEST_CASE( "collapse_cross_section threshold reaction" ) {
    const core::Spectrum<double> flux({1.0, 10.0, 100.0}, {2.0, 3.0});
    // the cross section starts at 10, so only the upper flux group counts
    const double rate = collapse_cross_section(flux, {10.0, 100.0}, {4.0, 4.0});
    REQUIRE_THAT( rate, WithinRel(12.0e-24, 1e-12) );
}

TEST_CASE( "collapse_cross_section bad input" ) {
    const core::Spectrum<double> flux({1.0, 10.0}, {2.0});
    REQUIRE_THROWS_AS( collapse_cross_section(flux, {1.0}, {1.0}), std::invalid_argument );
    REQUIRE_THROWS_AS( collapse_cross_section(flux, {1.0, 10.0}, {1.0}), std::invalid_argument );
}

TEST_CASE( "add_reaction_rate Fe56 (n,p)" ) {
    Eigen::SparseMatrix<double> mat(data::ZAI_SIZE, data::ZAI_SIZE);
    const double rate = 2.5e-6;
    add_reaction_rate(mat, data::ZAI_MAP, 260560, 103, rate);

    const int fe56 = data::ZAI_MAP.at(260560);
    const int mn56 = data::ZAI_MAP.at(250560);
    const int h1 = data::ZAI_MAP.at(PROTON);
    REQUIRE ( mat.coeff(fe56, fe56) == -rate );
    REQUIRE ( mat.coeff(fe56, mn56) == rate );
    REQUIRE ( mat.coeff(fe56, h1) == rate );
    REQUIRE ( mat.sum() == rate );
}

TEST_CASE( "add_reaction_rate fission to sink" ) {
    Eigen::SparseMatrix<double> mat(data::ZAI_SIZE, data::ZAI_SIZE);
    const double rate = 1e-4;
    add_reaction_rate(mat, data::ZAI_MAP, 922350, FISSION_MT, rate);

    const int u235 = data::ZAI_MAP.at(922350);
    REQUIRE ( mat.coeff(u235, u235) == -rate );
    REQUIRE ( mat.coeff(u235, data::SINK_INDEX) == rate );
    REQUIRE ( mat.sum() == 0.0 );
}

TEST_CASE( "add_reaction_rate ignores unknown parent and zero rate" ) {
    Eigen::SparseMatrix<double> mat(data::ZAI_SIZE, data::ZAI_SIZE);
    add_reaction_rate(mat, data::ZAI_MAP, -1, 102, 1.0);
    add_reaction_rate(mat, data::ZAI_MAP, 260560, 102, 0.0);
    REQUIRE ( mat.nonZeros() == 0 );
}

TEST_CASE( "reaction matrix with cram48 - single (n,g)" ) {
    Eigen::SparseMatrix<double> mat(data::ZAI_SIZE, data::ZAI_SIZE);
    const double rate = 1e-5;
    add_reaction_rate(mat, data::ZAI_MAP, 260560, 102, rate);

    const int fe56 = data::ZAI_MAP.at(260560);
    const int fe57 = data::ZAI_MAP.at(260570);
    Eigen::VectorXd n0 = Eigen::VectorXd::Zero(data::ZAI_SIZE);
    n0(fe56) = 1e20;

    const double dt = 1e4;
    const Eigen::VectorXd nt = core::cram48::solve<double>(mat, n0, dt);
    const double expected_fe56 = 1e20*std::exp(-rate*dt);
    REQUIRE_THAT( nt(fe56), WithinRel(expected_fe56, 1e-10) );
    REQUIRE_THAT( nt(fe57), WithinRel(1e20 - expected_fe56, 1e-10) );
}

TEST_CASE( "construct_reaction_matrix from reaction data" ) {
    const ReactionData reaction_data = {
        {260560, 102, {1.0, 10.0, 100.0}, {2.0, 3.0, 3.0}},
        {260560, 103, {10.0, 100.0}, {0.5, 0.5}},
        // elastic is not used
        {260560, 2, {1.0, 100.0}, {10.0, 10.0}},
    };
    const core::Spectrum<double> flux({1.0, 10.0, 100.0}, {1e12, 4e12});
    const Eigen::SparseMatrix<double> mat = construct_reaction_matrix(reaction_data, flux);

    const int fe56 = data::ZAI_MAP.at(260560);
    const int fe57 = data::ZAI_MAP.at(260570);
    const int mn56 = data::ZAI_MAP.at(250560);
    const double ng = (2.0*1e12 + 3.0*4e12)*BARN_TO_CM2;
    const double np = 0.5*4e12*BARN_TO_CM2;
    REQUIRE_THAT( mat.coeff(fe56, fe57), WithinRel(ng, 1e-12) );
    REQUIRE_THAT( mat.coeff(fe56, mn56), WithinRel(np, 1e-12) );
    REQUIRE_THAT( mat.coeff(fe56, fe56), WithinRel(-(ng + np), 1e-12) );
}

TEST_CASE( "split_rate_by_state" ) {
    // no MF=10, all to the ground state
    auto split = split_rate_by_state(10.0, {});
    REQUIRE( split == std::map<int, double>{{0, 10.0}} );

    // isomer only, the rest goes to the ground state
    split = split_rate_by_state(10.0, {{1, 4.0}});
    REQUIRE( split == std::map<int, double>{{0, 6.0}, {1, 4.0}} );

    // all states given, MF=3 total is kept, MF=10 gives the fractions
    split = split_rate_by_state(10.0, {{0, 2.0}, {1, 3.0}});
    REQUIRE_THAT( split.at(0), WithinRel(4.0, 1e-12) );
    REQUIRE_THAT( split.at(1), WithinRel(6.0, 1e-12) );

    // isomer more than the total, scale it to the total
    split = split_rate_by_state(10.0, {{1, 20.0}});
    REQUIRE( split == std::map<int, double>{{1, 10.0}} );
}

TEST_CASE( "construct_reaction_matrix with isomer production" ) {
    // Ni62 (n,p) -> Co62 and Co62m, values like JEFF-4
    const ReactionData reaction_data = {
        {280620, 103, {1e6, 2e7}, {0.1, 0.1}},
        {280620, 103, {1e6, 2e7}, {0.03, 0.03}, 0},
        {280620, 103, {1e6, 2e7}, {0.07, 0.07}, 1},
    };
    const core::Spectrum<double> flux({1e6, 2e7}, {1e10});
    const Eigen::SparseMatrix<double> mat = construct_reaction_matrix(reaction_data, flux);

    const double rate = 0.1*1e10*BARN_TO_CM2;
    const int ni62 = data::ZAI_MAP.at(280620);
    REQUIRE_THAT( mat.coeff(ni62, ni62), WithinRel(-rate, 1e-12) );
    REQUIRE_THAT( mat.coeff(ni62, data::ZAI_MAP.at(270620)), WithinRel(0.3*rate, 1e-12) );
    REQUIRE_THAT( mat.coeff(ni62, data::ZAI_MAP.at(270621)), WithinRel(0.7*rate, 1e-12) );
    REQUIRE_THAT( mat.coeff(ni62, data::ZAI_MAP.at(PROTON)), WithinRel(rate, 1e-12) );
}

TEST_CASE( "construct_reaction_matrix isomer not in the nuclide list goes to the ground state" ) {
    // Co62 has no 5th isomer
    REQUIRE( data::ZAI_MAP.count(270625) == 0 );
    const ReactionData reaction_data = {
        {280620, 103, {1e6, 2e7}, {0.1, 0.1}},
        {280620, 103, {1e6, 2e7}, {0.04, 0.04}, 5},
    };
    const core::Spectrum<double> flux({1e6, 2e7}, {1e10});
    const Eigen::SparseMatrix<double> mat = construct_reaction_matrix(reaction_data, flux);

    const double rate = 0.1*1e10*BARN_TO_CM2;
    const int ni62 = data::ZAI_MAP.at(280620);
    REQUIRE_THAT( mat.coeff(ni62, data::ZAI_MAP.at(270620)), WithinRel(rate, 1e-12) );
    REQUIRE( mat.coeff(ni62, data::SINK_INDEX) == 0.0 );
}

TEST_CASE( "compute_reaction_products with a final state" ) {
    check_products(270590, 102, 270600, {});
    const ReactionProductT p = compute_reaction_products(270590, 102, 1);
    REQUIRE( std::get<0>(p) == 270601 );
}

namespace {

    // A small reaction tape: Fe56 at two temperatures and Ni62.
    // The second Fe56 temperature must not be used.
    const std::string test_reaction_tape_string =
        "nis test reaction tape                                               1 0  0     \n"
        " 2.605600+4 5.545400+1          0          0          0          02631 1451     \n"
        " 0.000000+0 0.000000+0          0          0          0          62631 1451     \n"
        " 1.000000+0 2.000000+7          0          0         10          82631 1451     \n"
        " 2.936000+2 0.000000+0          0          0          0          22631 1451     \n"
        "                                3        102          5          02631 1451     \n"
        "                                3        103          5          02631 1451     \n"
        "                                                                  2631 1  0     \n"
        "                                                                  2631 0  0     \n"
        " 2.605600+4 5.545400+1          0          0          0          02631 3102     \n"
        " 0.000000+0 0.000000+0          0          0          1          32631 3102     \n"
        "          3          1                                            2631 3102     \n"
        " 1.000000+0 2.000000+0 1.000000+1 3.000000+0 1.000000+2 3.000000+02631 3102     \n"
        "                                                                  2631 3  0     \n"
        " 2.605600+4 5.545400+1          0          0          0          02631 3103     \n"
        " 0.000000+0 0.000000+0          0          0          1          22631 3103     \n"
        "          2          1                                            2631 3103     \n"
        " 1.000000+1 5.000000-1 1.000000+2 5.000000-1                      2631 3103     \n"
        "                                                                  2631 3  0     \n"
        "                                                                  2631 0  0     \n"
        "                                                                     0 0  0     \n"
        " 2.605600+4 5.545400+1          0          0          0          02631 1451     \n"
        " 0.000000+0 0.000000+0          0          0          0          62631 1451     \n"
        " 1.000000+0 2.000000+7          0          0         10          82631 1451     \n"
        " 6.000000+2 0.000000+0          0          0          0          12631 1451     \n"
        "                                3        102          5          02631 1451     \n"
        "                                                                  2631 1  0     \n"
        "                                                                  2631 0  0     \n"
        " 2.605600+4 5.545400+1          0          0          0          02631 3102     \n"
        " 0.000000+0 0.000000+0          0          0          1          22631 3102     \n"
        "          2          1                                            2631 3102     \n"
        " 1.000000+0 1.000000+2 1.000000+2 1.000000+2                      2631 3102     \n"
        "                                                                  2631 3  0     \n"
        "                                                                  2631 0  0     \n"
        "                                                                     0 0  0     \n"
        " 2.806200+4 6.139600+1          0          0          0          02837 1451     \n"
        " 0.000000+0 0.000000+0          0          0          0          62837 1451     \n"
        " 1.000000+0 2.000000+7          0          0         10          82837 1451     \n"
        " 2.936000+2 0.000000+0          0          0          0          12837 1451     \n"
        "                                3        103          5          02837 1451     \n"
        "                                                                  2837 1  0     \n"
        "                                                                  2837 0  0     \n"
        " 2.806200+4 6.139600+1          0          0          0          02837 3103     \n"
        " 0.000000+0 0.000000+0          0          0          1          22837 3103     \n"
        "          2          1                                            2837 3103     \n"
        " 1.000000+6 1.000000-1 2.000000+7 1.000000-1                      2837 3103     \n"
        "                                                                  2837 3  0     \n"
        "                                                                  2837 0  0     \n"
        "                                                                     0 0  0     \n"
        "                                                                    -1 0  0     \n";

    void check_test_reaction_data(const ReactionData& data){
        REQUIRE( data.size() == 3 );
        REQUIRE( data[0].parent_zai == 260560 );
        REQUIRE( data[0].mt_number == 102 );
        REQUIRE( data[0].energies == std::vector<double>{1.0, 10.0, 100.0} );
        REQUIRE( data[0].values == std::vector<double>{2.0, 3.0, 3.0} );
        REQUIRE( data[1].parent_zai == 260560 );
        REQUIRE( data[1].mt_number == 103 );
        REQUIRE( data[2].parent_zai == 280620 );
        REQUIRE( data[2].mt_number == 103 );
        REQUIRE( data[2].values == std::vector<double>{0.1, 0.1} );
    }
}

TEST_CASE( "read_reaction_data_from_tape" ) {
    const njoy::ENDFtk::tree::Tape tape(test_reaction_tape_string);
    check_test_reaction_data(read_reaction_data_from_tape(tape));
}

TEST_CASE( "read_reaction_data_from_file reads one material at a time" ) {
    const auto filename = std::filesystem::temp_directory_path() / "nis_test_reaction_tape.endf";
    {
        std::ofstream out(filename);
        out << test_reaction_tape_string;
    }
    const ReactionData data = read_reaction_data_from_file(filename.string());
    std::filesystem::remove(filename);
    check_test_reaction_data(data);
}

TEST_CASE( "read_reaction_data_from_file with a missing file" ) {
    REQUIRE_THROWS_AS( read_reaction_data_from_file("no_such_file.endf"), std::runtime_error );
}
