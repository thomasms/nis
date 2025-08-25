#include <catch2/catch_test_macros.hpp>

#include "nis/spectrum.hpp"

using namespace nis::core;


TEST_CASE( "get_overlap for no overlap" ) {
    const Edge<double> edge1(0.0, 1.0);
    const Edge<double> edge2(1.0, 2.0);
    std::tuple<bool, Edge<double>> overlap = get_overlap(edge1, edge2);
    const Edge<double> overlap_edge = std::get<1>(overlap);
    REQUIRE ( std::get<0>(overlap) == false);
    REQUIRE ( overlap_edge.lower() == 1.0);
    REQUIRE ( overlap_edge.upper() == 1.0);
}

TEST_CASE( "get_overlap for no overlap2" ) {
    const Edge<double> edge1(0.0, 0.5);
    const Edge<double> edge2(1.0, 2.0);
    std::tuple<bool, Edge<double>> overlap = get_overlap(edge1, edge2);
    const Edge<double> overlap_edge = std::get<1>(overlap);
    REQUIRE ( std::get<0>(overlap) == false);
    REQUIRE ( overlap_edge.lower() == 0.5);
    REQUIRE ( overlap_edge.upper() == 1.0);
}

TEST_CASE( "get_overlap for with overlap" ) {
    const Edge<double> edge1(0.0, 1.0);
    const Edge<double> edge2(0.9, 2.0);
    std::tuple<bool, Edge<double>> overlap = get_overlap(edge1, edge2);
    const Edge<double> overlap_edge = std::get<1>(overlap);
    REQUIRE ( std::get<0>(overlap) == true);
    REQUIRE ( overlap_edge.lower() == 0.9);
    REQUIRE ( overlap_edge.upper() == 1.0);
}

TEST_CASE( "get_overlap for with overlap 2" ) {
    const Edge<double> edge1(1.0, 2.0);
    const Edge<double> edge2(0.2, 1.1);
    std::tuple<bool, Edge<double>> overlap = get_overlap(edge1, edge2);
    const Edge<double> overlap_edge = std::get<1>(overlap);
    REQUIRE ( std::get<0>(overlap) == true);
    REQUIRE ( overlap_edge.lower() == 1.0);
    REQUIRE ( overlap_edge.upper() == 1.1);
}

TEST_CASE( "get_overlap for with overlap 3" ) {
    const Edge<double> edge1(1.0, 2.0);
    const Edge<double> edge2(1.2, 1.4);
    std::tuple<bool, Edge<double>> overlap = get_overlap(edge1, edge2);
    const Edge<double> overlap_edge = std::get<1>(overlap);
    REQUIRE ( std::get<0>(overlap) == true);
    REQUIRE ( overlap_edge.lower() == 1.2);
    REQUIRE ( overlap_edge.upper() == 1.4);
}
