#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>

#include "nis/spectrum.hpp"

using namespace nis::core;

using Catch::Matchers::Approx;

void check_overlap(const Edge<double>& edge1, const Edge<double>& edge2, const Edge<double>& expected_edge, bool has_overlap){
    std::tuple<bool, Edge<double>> overlap = get_overlap(edge1, edge2);
    std::tuple<bool, Edge<double>> roverlap = get_overlap(edge2, edge1);
    const Edge<double> overlap_edge = std::get<1>(overlap);
    const Edge<double> roverlap_edge = std::get<1>(roverlap);
    REQUIRE ( std::get<0>(overlap) == has_overlap);
    REQUIRE ( overlap_edge.lower() == expected_edge.lower());
    REQUIRE ( overlap_edge.upper() == expected_edge.upper());
    REQUIRE ( std::get<0>(roverlap) == has_overlap);
    REQUIRE ( roverlap_edge.lower() == expected_edge.lower());
    REQUIRE ( roverlap_edge.upper() == expected_edge.upper());
}

TEST_CASE( "get_overlap for no overlap" ) {
    const Edge<double> e1(0.0, 1.0);
    const Edge<double> e2(1.0, 2.0);
    const Edge<double> expected(1.0, 1.0);
    check_overlap(e1, e2, expected, false);
}

TEST_CASE( "get_overlap for no overlap2" ) {
    const Edge<double> e1(0.0, 0.5);
    const Edge<double> e2(1.0, 2.0);
    const Edge<double> expected(0.5, 1.0);
    check_overlap(e1, e2, expected, false);
}

TEST_CASE( "get_overlap for with overlap" ) {
    const Edge<double> e1(0.0, 1.0);
    const Edge<double> e2(0.9, 2.0);
    const Edge<double> expected(0.9, 1.0);
    check_overlap(e1, e2, expected, true);
}

TEST_CASE( "get_overlap for with overlap 2" ) {
    const Edge<double> e1(1.0, 2.0);
    const Edge<double> e2(0.2, 1.1);
    const Edge<double> expected(1.0, 1.1);
    check_overlap(e1, e2, expected, true);
}

TEST_CASE( "get_overlap for with overlap 3" ) {
    const Edge<double> e1(1.0, 2.0);
    const Edge<double> e2(1.2, 1.4);
    const Edge<double> expected(1.2, 1.4);
    check_overlap(e1, e2, expected, true);
}

using SpectrumD = Spectrum<double>;
using VecD = std::vector<double>;

TEST_CASE( "convert_by_energy" ) {
    SpectrumD spectrum(VecD{0.0, 1.0}, VecD{1.0});
    SpectrumD new_spec = convert_by_energy<double>(spectrum, VecD{0.5, 1.0});
    REQUIRE ( VecD{0.5} == new_spec.y());

    new_spec = convert_by_energy<double>(spectrum, VecD{0.0, 0.5});
    REQUIRE ( VecD{0.5} == new_spec.y());

    new_spec = convert_by_energy<double>(spectrum, VecD{1.0, 1.5});
    REQUIRE ( VecD{0.0} == new_spec.y());

    new_spec = convert_by_energy<double>(spectrum, VecD{-1.0, 0.0});
    REQUIRE ( VecD{0.0} == new_spec.y());

    SpectrumD spectrum2(VecD{0.0, 1.0, 2.0}, VecD{1.0, 0.5});
    new_spec = convert_by_energy<double>(spectrum2, VecD{0.0, 1.0, 2.0});
    REQUIRE ( VecD{1.0, 0.5} == new_spec.y());
    new_spec = convert_by_energy<double>(spectrum2, VecD{1.0, 2.0});
    REQUIRE ( VecD{0.5} == new_spec.y());
    new_spec = convert_by_energy<double>(spectrum2, VecD{1.5, 2.0});
    REQUIRE ( VecD{0.25} == new_spec.y());
    new_spec = convert_by_energy<double>(spectrum2, VecD{2.0, 2.0});
    REQUIRE ( VecD{0.} == new_spec.y());
    new_spec = convert_by_energy<double>(spectrum2, VecD{2.0, 10000.0});
    REQUIRE ( VecD{0.} == new_spec.y());
    new_spec = convert_by_energy<double>(spectrum2, VecD{0.0, 10000.0});
    REQUIRE ( VecD{1.5} == new_spec.y());

    SpectrumD spectrum3(VecD{-1.0, 0.0, 1.0}, VecD{5.0, 8.0});
    new_spec = convert_by_energy<double>(spectrum3, VecD{0.0, 0.5, 0.75, 1.0});
    REQUIRE ( VecD{4.0, 2.0, 2.0} == new_spec.y());

    new_spec = convert_by_energy<double>(spectrum3, VecD{-10.0, 0.5, 0.75, 1.0});
    REQUIRE ( VecD{9.0, 2.0, 2.0} == new_spec.y());
}

TEST_CASE( "convert_by_lethargy" ) {
    SpectrumD spectrum(VecD{0.0, 1.0}, VecD{1.0});
    SpectrumD new_spec = convert_by_lethargy<double>(spectrum, VecD{0.5, 1.0});
    REQUIRE ( VecD{0.69314718055994529} == new_spec.y());

    new_spec = convert_by_lethargy<double>(spectrum, VecD{0.0, 0.5});
    REQUIRE ( VecD{1.0} == new_spec.y());

    new_spec = convert_by_lethargy<double>(spectrum, VecD{1.0, 1.5});
    REQUIRE ( VecD{0.0} == new_spec.y());

    new_spec = convert_by_lethargy<double>(spectrum, VecD{-1.0, 0.0});
    REQUIRE ( VecD{0.0} == new_spec.y());

    SpectrumD spectrum2(VecD{0.0, 1.0, 2.0}, VecD{1.0, 0.2});
    new_spec = convert_by_lethargy<double>(spectrum2, VecD{0.0, 1.0, 2.0});
    // when using Approx we must use REQUIRE_THAT,
    // not REQUIRE: https://android.googlesource.com/platform/external/catch2/+/2cc0c718/docs/matchers.md
    REQUIRE_THAT ( VecD({1.0, 0.2}), Approx(new_spec.y()) );
    new_spec = convert_by_lethargy<double>(spectrum2, VecD{1.0, 2.0});
    REQUIRE_THAT ( VecD({0.2}), Approx(new_spec.y()) );
    new_spec = convert_by_lethargy<double>(spectrum2, VecD{1.5, 2.0});
    // 0.2*log(2.0/1.5)/(log(2/1))
    REQUIRE ( VecD{0.08300749985576875} == new_spec.y());

    SpectrumD spectrum3(VecD{-1.0, 0.0, 1.0}, VecD{5.0, 8.0});
    new_spec = convert_by_lethargy<double>(spectrum3, VecD{0.0, 0.5, 0.75, 1.0});
    REQUIRE ( VecD{8.0, 3.24372086486531508, 2.30145657961424677} == new_spec.y());

    new_spec = convert_by_lethargy<double>(spectrum3, VecD{-10.0, 0.5, 0.75, 1.0});
    REQUIRE ( VecD{13.0, 3.24372086486531508, 2.30145657961424677} == new_spec.y());
}