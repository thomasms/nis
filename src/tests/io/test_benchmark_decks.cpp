#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "nis/model.hpp"
#include "nis/io/helpers.hpp"

using namespace nis;

namespace {
    // set by CMake to the benchmarks folder of the source tree
    const std::filesystem::path benchmarks_dir = std::filesystem::path(NIS_BENCHMARKS_DIR);

    // benchmark name -> its decks (<case>.inp)
    std::map<std::string, std::vector<std::filesystem::path>> all_decks(){
        std::map<std::string, std::vector<std::filesystem::path>> decks;
        for(const auto& entry: std::filesystem::directory_iterator(benchmarks_dir)){
            if(!std::filesystem::exists(entry.path() / "benchmark.json")){
                continue;
            }
            auto& list = decks[entry.path().filename().string()];
            for(const auto& file: std::filesystem::recursive_directory_iterator(entry.path())){
                if(file.path().extension() == ".inp"){
                    list.push_back(file.path());
                }
            }
        }
        return decks;
    }
}

TEST_CASE( "all benchmark decks can be read", "[input][benchmark]" ) {
    const auto decks = all_decks();
    REQUIRE( decks.count("fns") == 1 );
    REQUIRE( decks.at("fns").size() == 132 );

    for(const auto& [benchmark, list]: decks){
        for(const auto& deck: list){
            INFO( deck.string() );
            // each deck is made from the case file next to it
            REQUIRE( std::filesystem::exists(std::filesystem::path(deck).replace_extension(".json")) );

            const auto model = io::helpers::from_file<double, core::BasicModel>(deck.string());
            REQUIRE_FALSE( model.nuclides.empty() );
            REQUIRE_FALSE( model.steps.empty() );
            REQUIRE( model.dt == model.steps[0].dt );

            bool irradiated = false;
            for(const auto& step: model.steps){
                REQUIRE( step.dt > 0.0 );
                REQUIRE( step.flux.has_value() );
                REQUIRE( step.flux.value() >= 0.0 );
                irradiated = irradiated || step.flux.value() > 0.0;
            }
            // a flux needs a spectrum
            if(irradiated){
                REQUIRE( model.spectrum != nullptr );
                REQUIRE( model.spectrum->x().size() == model.spectrum->y().size() + 1 );
            }
        }
    }
}

TEST_CASE( "the FNS decks: irradiation, then cooling", "[input][benchmark]" ) {
    // keep the map, so the list is not destroyed before the loop
    const auto decks = all_decks();
    for(const auto& deck: decks.at("fns")){
        INFO( deck.string() );
        const auto model = io::helpers::from_file<double, core::BasicModel>(deck.string());
        // 709 groups, so 710 bounds
        REQUIRE( model.spectrum != nullptr );
        REQUIRE( model.spectrum->x().size() == 710 );
        REQUIRE( model.steps.size() >= 2 );
        REQUIRE( model.steps[0].flux.value() > 0.0 );
        for(std::size_t i = 1; i < model.steps.size(); ++i){
            REQUIRE( model.steps[i].flux.value() == 0.0 );
        }
    }
}
