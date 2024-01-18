#include <catch2/catch_test_macros.hpp>

#include "nis/endf/decay.hpp"

using namespace nis::endf;

TEST_CASE( "Expand the rtype from endf", "[expand_rtypes]" ) {
    REQUIRE( expand_rtypes(1.2) == std::vector<int>{1, 2} );
    REQUIRE( expand_rtypes(1.0) == std::vector<int>{1} );
    REQUIRE( expand_rtypes(0) == std::vector<int>{0} );
    REQUIRE( expand_rtypes(0.1) == std::vector<int>{0} );
    REQUIRE( expand_rtypes(4.5552) == std::vector<int>{4,5,5,5,2} );
    REQUIRE( expand_rtypes(1.30231) == std::vector<int>{1,3} );
    REQUIRE( expand_rtypes(10) == std::vector<int>{10} );
    // this is a bug
    REQUIRE( expand_rtypes(10.1) == std::vector<int>{10, 1} );
    // 11 is not mapped
    REQUIRE( expand_rtypes(11) == std::vector<int>{} );
    REQUIRE( expand_rtypes(11.2) == std::vector<int>{} );
    // 30 is not mapped
    REQUIRE( expand_rtypes(30) == std::vector<int>{} );
}