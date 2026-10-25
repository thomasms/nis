#include <iostream>
#include <string>

#include <CLI/CLI.hpp>

#include "nis/endf/reaction.hpp"
#include "nis/endf/reaction_cache.hpp"

using namespace nis;

int main(int argc, char **argv){

    std::string reaction_data = "";
    std::string output = "";

    CLI::App app{"Make a binary cache of the reaction data, nis reads it much faster than the ENDF tape"};
    app.add_option("-r,--reaction-data", reaction_data, "An single ENDF file containing all reaction data")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", output, "The cache file to write, give it to nis with --reaction-data")
      ->required();
    CLI11_PARSE(app, argc, argv);

    try{
      const auto data = endf::read_reaction_data_from_file(reaction_data);
      endf::write_reaction_cache(data, output);
      std::cout << data.size() << " cross sections written to " << output << "\n";
    }
    catch(const std::exception& ex){
      std::cerr << ex.what();
      return 1;
    }

    return 0;
}
