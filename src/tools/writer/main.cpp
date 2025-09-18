#include <memory>
#include <vector>

#include <CLI/CLI.hpp>

#include "nis/io/helpers.hpp"
#include "nis/io/input.hpp"

using namespace nis::core;
using namespace nis::io;

using Model = BasicModel<double>;


int main(int argc, char **argv){

    std::string filename = "example.inp";

    CLI::App app{"NIS"};
    app.add_option("-f,--file", filename, "The name of the output file")
      ->required();
    CLI11_PARSE(app, argc, argv);

    Model data;
    data.dt = 30000.0;
    data.nuclides = {
        {"u234", 1e+20},
        {"u235", 1e+20},
        {"u236", 1e+20},
        {"u237", 1e+20},
        {"u238", 1e+20},
        {"u239", 1e+20}
    };
    data.spectrum = std::make_shared<Spectrum<double>>(std::move(std::vector<double>{0.0, 10, 20}),
                                                       std::move(std::vector<double>{0.5e16, 0.1e16}));

    helpers::to_file<double, BasicModel>(data, filename);
    return 0;
}
