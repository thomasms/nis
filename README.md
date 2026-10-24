# NIS


Nuclear Inventory Simulator


## Dependencies

- Eigen (For matrix manipulation and solvers)
- ENDFtk (For reading ENDF files)
- Catch2 (For unit testing)

ENDFtk has its own dependencies, which are all header only but need to be pulled locally too.

All dependencies are header only, so only include paths need to be set.


If you follow the script in `docker/setup_local_deps.sh`, running it in the root directory of this repository, it should setup all the correct dependencies as needed.


## Build

Once you have all the dependencies locally, you can then compile via:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```


### python bindings

If you want the python bindings (default=OFF), then follow the below steps.


First make a virtual env and install pybind via:

```bash
python3 -m venv venv
source venv/bin/activate
python3 -m pip install pybind11 numpy scipy
```

Then build as before but with slightly different cmake args, as below.
```bash
source venv/bin/activate
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_PY_API=ON -Dpybind11_DIR=`pybind11-config --cmakedir` ..
make -j4
```


## Example

To test your build, after compiling as above, try the following:

```bash
./build/bin/nis --decay-data ./data/jeff311.dec --reaction-data ./data/jeff4.reac --input examples/h3_halflife.inp
```

If successful, you should see roughly half remains as H3 and the other half has beta decayed to He3.

```bash
H3 5.00518e+19
He3 4.99482e+19
```


After the nuclides, the output gives the total activity (Bq) and the decay heat (W).
The decay heat comes from the average decay energies in the decay data (MF=8, MT=457), split into
light particles (beta, electrons), electromagnetic (gamma, x-rays) and heavy particles (alpha, recoil).

For neutron activation, add a `spectrum` to the input (see [examples/activate.inp](examples/activate.inp)) and give the reaction data too:

```bash
./build/bin/nis --decay-data ./data/jeff311.dec --reaction-data ./data/jeff4.reac --input examples/activate.inp
```

The spectrum is the neutron flux in n/cm^2/s per group. The first line gives the group bounds in eV (ascending), the second line gives one value per group.
The reaction data must be group-wise (INT=1), for example as processed by NJOY GROUPR.

To decay the inventory after the irradiation, add a `cooling` line with the length of each step in seconds:

```
cooling
    60 600 3600
```

The output then gives the nuclides, activity and decay heat after the irradiation, and again after each
cooling step, under a `cooling <time>` line (the time since the end of the irradiation, in seconds).

Note, you will need to have pulled the decay data via LFS in the [./data](./data) subdirectory.


## Tests

C++ unit tests:

```bash
./build/bin/nistests
```

Python integration tests and benchmarks (needs the python bindings and the data from LFS):

```bash
export PYTHONPATH=./build/bin/python3:./src/tests/python
pytest -v src/tests/python/testsuite.py
```

The reaction benchmarks in [src/tests/python/reactiontests.py](src/tests/python/reactiontests.py) use JEFF-4.0 and check:
- thermal (0.0253 eV) and 14.1 MeV cross sections against reference values,
- spectrum handling (linear in flux, bin splitting),
- the structure of the full reaction matrix,
- burn-up and Co60 activation against the analytic Bateman solution.

### FNS decay heat benchmark

The FNS (JAEA) experiments irradiate foils with 14 MeV neutrons and measure the decay heat.
[scripts/fns_benchmark.py](scripts/fns_benchmark.py) runs NIS for each case of the
FISPACT-II validation inputs, and compares it with the measured heat and with the FISPACT-II
(TENDL-2017) results in those files. FISPACT-II is not run.

```bash
export PYTHONPATH=./build/bin/python3
python scripts/fns_benchmark.py --fns-dir path/to/fns \
    --decay-data data/jeff311.dec --reaction-data data/jeff4.reac \
    --output fns.csv --nuclides-output fns_nuclides.csv
# C/E statistics for NIS and FISPACT-II
python scripts/fns_stats.py fns.csv
# the nuclides that cause a difference between NIS and FISPACT-II
python scripts/fns_heat_by_nuclide.py fns_nuclides.csv --decay-data data/jeff311.dec
# plots of NIS (smooth curve) and FISPACT-II against the experiment, one per case
python scripts/fns_plot.py --fns-dir path/to/fns \
    --decay-data data/jeff311.dec --reaction-data data/jeff4.reac --output-dir fns_plots
```

The cases run in parallel, one process per core by default (`--jobs`). The plots need `matplotlib`.

NIS input decks for all the cases are in [benchmarks/fns](benchmarks/fns), one `<material>/<case>.inp` for
each case (initial atoms, irradiation time, 709 group flux, a cooling step to each measured time), with the
measured and FISPACT-II decay heat in `<case>.csv`. They are made by
[scripts/fns_make_decks.py](scripts/fns_make_decks.py), and run with the `nis` program by
[scripts/fns_run_decks.py](scripts/fns_run_decks.py):

```bash
./build/bin/nis -d data/jeff311.dec -r data/jeff4.reac -i benchmarks/fns/Co/2000exp_5min.inp
python scripts/fns_run_decks.py --nis build/bin/nis --decks benchmarks/fns \
    --decay-data data/jeff311.dec --reaction-data data/jeff4.reac --output fns_decks.csv
python scripts/fns_stats.py fns_decks.csv
```

Each `nis` run reads the nuclear data again, so `fns_run_decks.py` is slower than `fns_benchmark.py`.
