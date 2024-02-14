#include <CLI/CLI.hpp>

#include "nis/io/input.hpp"

using namespace nis::io;

int main(int argc, char **argv){

    std::string filename = "example.inp";

    CLI::App app{"NIS"};
    app.add_option("-f,--file", filename, "The name of the output file")
      ->required();
    CLI11_PARSE(app, argc, argv);

    input::InputData<double> data;
    data.dt = 30000.0;
    data.nuclides = {
        {"u234", 1e+20},
        {"u235", 1e+20},
        {"u236", 1e+20},
        {"u237", 1e+20},
        {"u238", 1e+20},
        {"u239", 1e+20}
    };

    input::to_file(data, filename);
    return 0;
}
