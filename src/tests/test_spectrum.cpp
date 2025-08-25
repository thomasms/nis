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

using SpectrumD = Spectrum<double>;

TEST_CASE( "convert_by_energy" ) {
    SpectrumD spectrum(std::vector<double>{0.0, 1.0}, std::vector<double>{1.0});
    SpectrumD new_spec = convert_by_energy<double>(spectrum, std::vector<double>{0.5, 1.0});
    REQUIRE ( std::vector<double>{0.5} == new_spec.y());

    new_spec = convert_by_energy<double>(spectrum, std::vector<double>{0.0, 0.5});
    REQUIRE ( std::vector<double>{0.5} == new_spec.y());

    new_spec = convert_by_energy<double>(spectrum, std::vector<double>{1.0, 1.5});
    REQUIRE ( std::vector<double>{0.0} == new_spec.y());

    new_spec = convert_by_energy<double>(spectrum, std::vector<double>{-1.0, 0.0});
    REQUIRE ( std::vector<double>{0.0} == new_spec.y());

    SpectrumD spectrum2(std::vector<double>{0.0, 1.0, 2.0}, std::vector<double>{1.0, 0.7});
    new_spec = convert_by_energy<double>(spectrum2, std::vector<double>{0.0, 1.0, 2.0});
    REQUIRE ( std::vector<double>{1.0, 0.7} == new_spec.y());

    SpectrumD spectrum3(std::vector<double>{-1.0, 0.0, 1.0}, std::vector<double>{5.0, 8.0});
    new_spec = convert_by_energy<double>(spectrum3, std::vector<double>{0.0, 0.5, 0.75, 1.0});
    REQUIRE ( std::vector<double>{4.0, 2.0, 2.0} == new_spec.y());

    new_spec = convert_by_energy<double>(spectrum3, std::vector<double>{-10.0, 0.5, 0.75, 1.0});
    REQUIRE ( std::vector<double>{9.0, 2.0, 2.0} == new_spec.y());
}

TEST_CASE( "convert_by_lethargy" ) {
    SpectrumD spectrum(std::vector<double>{0.0, 1.0}, std::vector<double>{1.0});
    SpectrumD new_spec = convert_by_lethargy<double>(spectrum, std::vector<double>{0.5, 1.0});
    REQUIRE ( std::vector<double>{0.69314718055994529} == new_spec.y());

    new_spec = convert_by_lethargy<double>(spectrum, std::vector<double>{0.0, 0.5});
    REQUIRE ( std::vector<double>{1.0} == new_spec.y());

    new_spec = convert_by_lethargy<double>(spectrum, std::vector<double>{1.0, 1.5});
    REQUIRE ( std::vector<double>{0.0} == new_spec.y());

    new_spec = convert_by_lethargy<double>(spectrum, std::vector<double>{-1.0, 0.0});
    REQUIRE ( std::vector<double>{1.0} == new_spec.y());

    SpectrumD spectrum2(std::vector<double>{0.0, 1.0, 2.0}, std::vector<double>{1.0, 0.7});
    new_spec = convert_by_lethargy<double>(spectrum2, std::vector<double>{0.0, 1.0, 2.0});
    REQUIRE ( std::vector<double>{1.0, 0.7} == new_spec.y());

    SpectrumD spectrum3(std::vector<double>{-1.0, 0.0, 1.0}, std::vector<double>{5.0, 8.0});
    new_spec = convert_by_lethargy<double>(spectrum3, std::vector<double>{0.0, 0.5, 0.75, 1.0});
    REQUIRE ( std::vector<double>{8.0, 3.24372086486531508, 2.30145657961424677} == new_spec.y());

    new_spec = convert_by_lethargy<double>(spectrum3, std::vector<double>{-10.0, 0.5, 0.75, 1.0});
    REQUIRE ( std::vector<double>{13.0, 3.24372086486531508, 2.30145657961424677} == new_spec.y());
}