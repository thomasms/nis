#include <catch2/catch_test_macros.hpp>

#include "zai.hpp"

using namespace nis::zai;

TEST_CASE( "Z from ZAI", "[get_z_from_zai]" ) {
    REQUIRE( get_z_from_zai(0) == 0 );
    REQUIRE( get_z_from_zai(1) == 0 );
    REQUIRE( get_z_from_zai(100) == 0 );
    REQUIRE( get_z_from_zai(10010) == 1 );
    REQUIRE( get_z_from_zai(120320) == 12 );
    REQUIRE( get_z_from_zai(1234567) == 123 );
    REQUIRE( get_z_from_zai(91234567) == 9123 );
}

TEST_CASE( "A from ZAI", "[get_a_from_zai]" ) {
    REQUIRE( get_a_from_zai(0) == 0 );
    REQUIRE( get_a_from_zai(1) == 0 );
    REQUIRE( get_a_from_zai(10) == 1 );
    REQUIRE( get_a_from_zai(430) == 43 );
    REQUIRE( get_a_from_zai(100) == 10 );
    REQUIRE( get_a_from_zai(10010) == 1 );
    REQUIRE( get_a_from_zai(120320) == 32 );
    REQUIRE( get_a_from_zai(1234567) == 456 );
    REQUIRE( get_a_from_zai(91234567) == 456 );
}

TEST_CASE( "A from ZAI", "[get_i_from_zai]" ) {
    REQUIRE( get_i_from_zai(0) == 0 );
    REQUIRE( get_i_from_zai(1) == 1 );
    REQUIRE( get_i_from_zai(10) == 0 );
    REQUIRE( get_i_from_zai(435) == 5 );
    REQUIRE( get_i_from_zai(101) == 1 );
    REQUIRE( get_i_from_zai(10012) == 2 );
    REQUIRE( get_i_from_zai(120320) == 0 );
    REQUIRE( get_i_from_zai(120329) == 9 );
    REQUIRE( get_i_from_zai(1234567) == 7 );
    REQUIRE( get_i_from_zai(1234560) == 0 );
    REQUIRE( get_i_from_zai(91234560) == 0 );
    REQUIRE( get_i_from_zai(91234567) == 7 );
}