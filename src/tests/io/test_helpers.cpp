#include <catch2/catch_test_macros.hpp>

#include <iostream>
#include <map>

#include <Eigen/Core>

#include "nis/nuclides.hpp"
#include "nis/util.hpp"
#include "nis/io/helpers.hpp"

using namespace nis;
using namespace nis::io;

TEST_CASE( "Vec to inventory", "[vector_from_inventory]" ) {
    const std::map<std::string, double> inv = {
        {"He3", 23.4},
        {"Mg27", 78.345e4},
        {"Xe123", 45.67},
        {"NotValid", 235.07},
        // empty string should not be a problem
        {"", 1.03},
    };
    const Eigen::VectorX<double> vec = helpers::vector_from_inventory<double>(inv);

    const int he3_index = data::ZAI_MAP.at(20030);
    const int mg27_index = data::ZAI_MAP.at(120270);
    const int Xe123_index = data::ZAI_MAP.at(541230);

    REQUIRE( vec.size() == data::ZAI_SIZE );
    REQUIRE( vec(he3_index) == 23.4 );
    REQUIRE( vec(mg27_index) == 78.345e4 );
    REQUIRE( vec(Xe123_index) == 45.67 );
    // REQUIRE( vec(sink_index) == 235.07 );
    for(int i=0;i<vec.size();++i){
        if((i != he3_index) && (i != mg27_index) && (i != Xe123_index)){
            REQUIRE( vec(i) == 0.0 );
        }
    }
}

TEST_CASE( "Inventory to vec", "[inventory_from_vector]" ) {
    const int he3_index = data::ZAI_MAP.at(20030);
    const int mg27_index = data::ZAI_MAP.at(120270);
    const int Xe123_index = data::ZAI_MAP.at(541230);
    Eigen::VectorX<double> vec = Eigen::VectorX<double>::Zero(data::ZAI_SIZE);
    vec(he3_index) = 23.4;
    vec(mg27_index) = 78.345e4;
    vec(Xe123_index) = 45.67;
    vec(data::SINK_INDEX) = 235.07;

    const std::map<std::string, double> inv = helpers::inventory_from_vector<double>(vec);

    REQUIRE( inv.size() == 4 );
    // empty string is unmapped
    REQUIRE( inv.at("") == 235.07 );
    REQUIRE( inv.at("He3") == 23.4 );
    REQUIRE( inv.at("Mg27") == 78.345e4 );
    REQUIRE( inv.at("Xe123")  == 45.67 );
}
TEST_CASE( "Input with spectrum and cooling", "[input]" ) {
    std::istringstream is(
        "dt 3600\n"
        "\n"
        "nuclides\n"
        "    Co59 1e+20\n"
        "\n"
        "spectrum\n"
        "    0 10 20\n"
        "    5e+15 1e+15\n"
        "\n"
        "cooling\n"
        "    60 120 3e+06\n");
    core::BasicModel<double> model;
    is >> model;

    REQUIRE( model.dt == 3600.0 );
    REQUIRE( model.nuclides.size() == 1 );
    REQUIRE( model.nuclides.at("Co59") == 1e20 );
    REQUIRE( model.spectrum != nullptr );
    REQUIRE( model.spectrum->x() == std::vector<double>{0.0, 10.0, 20.0} );
    REQUIRE( model.spectrum->y() == std::vector<double>{5e15, 1e15} );
    REQUIRE( model.cooling == std::vector<double>{60.0, 120.0, 3e6} );

    // write it and read it again
    std::ostringstream os;
    os << model;
    std::istringstream again(os.str());
    core::BasicModel<double> copy;
    again >> copy;
    REQUIRE( copy.dt == model.dt );
    REQUIRE( copy.nuclides == model.nuclides );
    REQUIRE( copy.cooling == model.cooling );
}

TEST_CASE( "Input without cooling", "[input]" ) {
    std::istringstream is("dt 10\n\nnuclides\n    H3 1e+10\n");
    core::BasicModel<double> model;
    is >> model;
    REQUIRE( model.cooling.empty() );
    REQUIRE( model.spectrum == nullptr );
}
