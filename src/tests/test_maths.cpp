#include <catch2/catch_test_macros.hpp>

#include "nis/maths.hpp"

using namespace nis::maths;

TEST_CASE( "dot_product", "[dot_product]" ) {
    const auto res = dot_product<double>({9.2,3,7}, {1,-3,10.5});

    REQUIRE( res == 73.7 );
}
