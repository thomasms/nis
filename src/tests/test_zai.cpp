#include <catch2/catch_test_macros.hpp>

#include "nis/zai.hpp"

using namespace nis::zai;

TEST_CASE( "Z from ZAI", "[get_z]" ) {
    REQUIRE( get_z(0) == 0 );
    REQUIRE( get_z(1) == 0 );
    REQUIRE( get_z(100) == 0 );
    REQUIRE( get_z(10010) == 1 );
    REQUIRE( get_z(120320) == 12 );
    REQUIRE( get_z(1234567) == 123 );
    REQUIRE( get_z(91234567) == 9123 );
}

TEST_CASE( "A from ZAI", "[get_a]" ) {
    REQUIRE( get_a(0) == 0 );
    REQUIRE( get_a(1) == 0 );
    REQUIRE( get_a(10) == 1 );
    REQUIRE( get_a(430) == 43 );
    REQUIRE( get_a(100) == 10 );
    REQUIRE( get_a(10010) == 1 );
    REQUIRE( get_a(120320) == 32 );
    REQUIRE( get_a(1234567) == 456 );
    REQUIRE( get_a(91234567) == 456 );
}

TEST_CASE( "A from ZAI", "[get_i]" ) {
    REQUIRE( get_i(0) == 0 );
    REQUIRE( get_i(1) == 1 );
    REQUIRE( get_i(10) == 0 );
    REQUIRE( get_i(435) == 5 );
    REQUIRE( get_i(101) == 1 );
    REQUIRE( get_i(10012) == 2 );
    REQUIRE( get_i(120320) == 0 );
    REQUIRE( get_i(120329) == 9 );
    REQUIRE( get_i(1234567) == 7 );
    REQUIRE( get_i(1234560) == 0 );
    REQUIRE( get_i(91234560) == 0 );
    REQUIRE( get_i(91234567) == 7 );
}

TEST_CASE( "ZAI from Z, A, I", "[make]" ) {
    REQUIRE( make(1, 1, 0) == 10010 );
    REQUIRE( make(92, 238, 1) == 922381 );
    REQUIRE( make(103, 258, 2) == 1032582 );
    REQUIRE( make(118, 295, 0) == 1182950 );
    REQUIRE( make(118, 295, 1) == 1182951 );
}

TEST_CASE( "ZAI from triple", "[make]" ) {
    REQUIRE( make(triple(1, 1, 0)) == 10010 );
    REQUIRE( make(triple(92, 238, 1)) == 922381 );
    REQUIRE( make(triple(103, 258, 2)) == 1032582 );
    REQUIRE( make(triple(118, 295, 0)) == 1182950 );
    REQUIRE( make(triple(118, 295, 1)) == 1182951 );
}

TEST_CASE( "Z, A, I from ZAI", "[decompose]" ) {
    REQUIRE( decompose(10010) == triple(1, 1, 0) );
    REQUIRE( decompose(922381) == triple(92, 238, 1) );
    REQUIRE( decompose(1032582) == triple(103, 258, 2) );
    REQUIRE( decompose(1182950) == triple(118, 295, 0) );
    REQUIRE( decompose(1182951) == triple(118, 295, 1) );
}

TEST_CASE( "ZAI from name", "[get_from_name]" ) {
    REQUIRE( get_from_name("h1") == 10010 );
    REQUIRE( get_from_name("H1") == 10010 );
    REQUIRE( get_from_name("H1m") == 10011 );
    REQUIRE( get_from_name("H1u") == 10019 );
    REQUIRE( get_from_name("He4") == 20040 );
    REQUIRE( get_from_name("u235s") == 922357 );
    REQUIRE( get_from_name("ra223") == 882230 );
    REQUIRE( get_from_name("Ra234o") == 882343 );
}

TEST_CASE( "name from ZAI", "[get_name_from]" ) {
    REQUIRE( get_name_from(10010) == "H1" );
    REQUIRE( get_name_from(10011) == "H1m" );
    REQUIRE( get_name_from(10019) == "H1u" );
    REQUIRE( get_name_from(20040) == "He4" );
    REQUIRE( get_name_from(922357) == "U235s" );
    REQUIRE( get_name_from(882230) == "Ra223" );
    REQUIRE( get_name_from(882343) == "Ra234o" );
}