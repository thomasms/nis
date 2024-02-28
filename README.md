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

Once you have all the dependencies locally and have compiled SparseSuite, you can then compile via:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```


## Example

To test your build, after compiling as above, try the following:

```bash
./build/bin/nis --data ./data/jeff311.dec --input examples/h3_halflife.inp
```

If successful, you should see roughly half remains as H3 and the other half has beta decayed to He3.

```bash
H3 5.00518e+19
He3 4.99482e+19
```


Note, you will need to have pulled the decay data via LFS in the [./data](./data) subdirectory.
