#include <catch2/catch_test_macros.hpp>

#include "nis/util.hpp"

using namespace nis::util;

TEST_CASE( "Map from array", "[make_map_index_lookup]" ) {
    constexpr std::array<int, 3> arr = {9,3,7};
    const auto res = make_map_index_lookup(arr);

    REQUIRE( res.size() == 3 );
    REQUIRE( res.at(9) == 0 );
    REQUIRE( res.at(3) == 1 );
    REQUIRE( res.at(7) == 2 );
}

TEST_CASE( "Break up an integer into separate digits", "[get_digits]" ){

    REQUIRE( get_digits<int>(0) == std::list<int>{0} );
    REQUIRE( get_digits<int>(101) == std::list<int>{1,0,1} );
    REQUIRE( get_digits<int>(-5102) == std::list<int>{5,1,0,2} );
    REQUIRE( get_digits<int>(10102000) == std::list<int>{1,0,1,0,2,0,0,0} );
    REQUIRE( get_digits<long long>(12345678901239810) == std::list<int>{1,2,3,4,5,6,7,8,9,0,1,2,3,9,8,1,0} );
    REQUIRE( get_digits<unsigned int>(782111271) == std::list<int>{7,8,2,1,1,1,2,7,1} );
}