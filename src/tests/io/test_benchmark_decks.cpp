#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "nis/model.hpp"
#include "nis/io/helpers.hpp"

using namespace nis;

namespace {
    // set by CMake to the benchmarks folder of the source tree
    const std::filesystem::path decks_dir = std::filesystem::path(NIS_BENCHMARKS_DIR) / "fns";

    std::vector<std::filesystem::path> all_decks(){
        std::vector<std::filesystem::path> decks;
        for(const auto& entry: std::filesystem::recursive_directory_iterator(decks_dir)){
            if(entry.path().extension() == ".inp"){
                decks.push_back(entry.path());
            }
        }
        return decks;
    }

    // the number of measured times in the csv of the deck (the lines after the header)
    std::size_t measured_times(const std::filesystem::path& deck){
        std::ifstream in(std::filesystem::path(deck).replace_extension(".csv"));
        std::size_t lines = 0;
        std::string line;
        while(std::getline(in, line)){
            if(!line.empty()){
                ++lines;
            }
        }
        return lines > 0 ? lines - 1 : 0;
    }
}

TEST_CASE( "all FNS benchmark decks can be read", "[input][benchmark]" ) {
    const auto decks = all_decks();
    REQUIRE( decks.size() == 132 );

    for(const auto& deck: decks){
        INFO( deck.string() );
        const auto model = io::helpers::from_file<double, core::BasicModel>(deck.string());

        REQUIRE_FALSE( model.nuclides.empty() );
        REQUIRE( model.spectrum != nullptr );
        // 709 groups, so 710 bounds
        REQUIRE( model.spectrum->x().size() == 710 );
        REQUIRE( model.spectrum->y().size() == 709 );

        // the irradiation with a flux, then one cooling step to each measured time
        REQUIRE( model.steps.size() == measured_times(deck) + 1 );
        REQUIRE( model.steps[0].flux.has_value() );
        REQUIRE( model.steps[0].flux.value() > 0.0 );
        REQUIRE( model.dt == model.steps[0].dt );
        for(std::size_t i = 1; i < model.steps.size(); ++i){
            REQUIRE( model.steps[i].flux == 0.0 );
            REQUIRE( model.steps[i].dt > 0.0 );
        }
    }
}
