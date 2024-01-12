#include <catch2/catch_test_macros.hpp>

#include "util.hpp"

using namespace nis::util;

TEST_CASE( "Map from array", "[make_map_index_lookup]" ) {
    constexpr std::array<int, 3> arr = {9,3,7};
    const auto res = make_map_index_lookup(arr);

    REQUIRE( res.size() == 3 );
    REQUIRE( res.at(9) == 0 );
    REQUIRE( res.at(3) == 1 );
    REQUIRE( res.at(7) == 2 );
}
