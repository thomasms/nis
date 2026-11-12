#include <catch2/catch_test_macros.hpp>

#include "nis/endf/decay.hpp"
#include "nis/zai.hpp"

using namespace nis;
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


TEST_CASE( "compute_products Li7 beta" ) {
    // Li7 beta decay
    DecayProductT p = compute_products(30080, 1.0, 0);
    const int primary = std::get<0>(p);
    const std::vector<int> secondaries = std::get<1>(p);
    // primary -> Z+1 = 40080
    REQUIRE ( primary == 40080);
    // secondaries - no secondaries for beta minus
    REQUIRE ( secondaries == std::vector<int>{} );
}

TEST_CASE( "compute_products Li7 beta then alpha" ) {
    // Li7 beta decay
    DecayProductT p = compute_products(30080, 1.4, 0);
    const int primary = std::get<0>(p);
    const std::vector<int> secondaries = std::get<1>(p);
    // primary -> Z+1 = 40080,
    // then Z-2, N-4 = 20040,
    REQUIRE ( primary == 20040);
    // secondaries
    REQUIRE ( secondaries == std::vector<int>{zai::ALPHA} );
}

TEST_CASE( "compute_products Li7 beta then 2 alphas" ) {
    // Li7 beta decay
    DecayProductT p = compute_products(30080, 1.44, 0);
    const int primary = std::get<0>(p);
    const std::vector<int> secondaries = std::get<1>(p);
    // primary -> Z+1 = 40080,
    // then Z-2, N-4 = 20040,
    // then another Z-2, N-4 = 0
    REQUIRE ( primary == 0);
    // secondaries
    REQUIRE ( secondaries == std::vector<int>{zai::ALPHA, zai::ALPHA} );
}

// TODO: how to best handle these unphysical processes?
TEST_CASE( "compute_products Li7 beta then 3 alphas - not physically possible" ) {
    // Li7 beta decay
    DecayProductT p = compute_products(30080, 1.444, 0);
    const int primary = std::get<0>(p);
    const std::vector<int> secondaries = std::get<1>(p);
    // primary -> Z+1 = 40080,
    // then Z-2, N-4 = 20040,
    // then another Z-2, N-4 = 0
    // then another Z-2, N-4 = -20040
    REQUIRE ( primary == -20040);
    // secondaries
    REQUIRE ( secondaries == std::vector<int>{zai::ALPHA, zai::ALPHA, zai::ALPHA} );
}
TEST_CASE( "compute_products proton emission" ) {
    // Li5 -> He4 + p: Z-1 and A-1, and the proton is an H1 atom
    DecayProductT p = compute_products(30050, 7.0, 0);
    REQUIRE ( std::get<0>(p) == 20040 );
    REQUIRE ( std::get<1>(p) == std::vector<int>{zai::PROTON} );

    // Fe46 beta+ then proton: Fe46 -> Mn46 -> Cr45 + p
    p = compute_products(260460, 2.7, 0);
    REQUIRE ( std::get<0>(p) == 240450 );
    REQUIRE ( std::get<1>(p) == std::vector<int>{zai::PROTON} );
}

TEST_CASE( "compute_products neutron emission" ) {
    // Li9 beta- then neutron: Li9 -> Be9 -> Be8 + n. Neutrons are not tracked, so there is no secondary.
    DecayProductT p = compute_products(30090, 1.5, 0);
    REQUIRE ( std::get<0>(p) == 40080 );
    REQUIRE ( std::get<1>(p).empty() );
}
