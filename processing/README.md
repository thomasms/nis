# Data processing

Make NIS reaction data from the raw ENDF neutron files of a library (for example TENDL), with
[NJOY2016](https://github.com/njoy/NJOY2016) and [PREPRO](https://github.com/IAEA-NDS/PREPRO).
The result is one file of 1102 group cross sections (INT=1), the same form as `data/jeff4.reac`, and the NIS
cache of it (`.nisreac`).

```
processing/
    process.py          the command
    nisdata/            the code (a python package, only the standard library)
        pipeline.py     the steps for one nuclide
        decks.py        the NJOY and PREPRO inputs
        runner.py       many nuclides at the same time
        endf.py         read ENDF headers, join the processed files into one
        tools.py        get and build NJOY and PREPRO
        fetch.py        download TENDL
        ukaea1102.txt   the 1102 group bounds (eV)
    tools/              NJOY and PREPRO, made by 'tools' (not in git)
    raw/                the raw ENDF files, from 'fetch' (not in git)
    out/                the processed files, one for each nuclide (not in git)
```

## Quick start

You need git, cmake, make and gfortran. Python 3.9 or newer, no other packages.

```bash
python3 processing/process.py tools                       # get and build NJOY and PREPRO (about 1 minute)
python3 processing/process.py tendl --year 2025           # download, process all, join: data/tendl2025.reac
./build/bin/nis -d data/jeff311.dec -r data/tendl2025.nisreac -i examples/activate.inp
```

`tendl` does the 3 steps below. You can also do them one at a time:

```bash
python3 processing/process.py fetch --year 2025           # processing/raw/TENDL-2025-n/
python3 processing/process.py run processing/raw/TENDL-2025-n --output processing/out/tendl2025
python3 processing/process.py combine processing/out/tendl2025 --output data/tendl2025.reac \
    --title "TENDL-2025 Incident Neutron File (1102 groups)"
```

Any library works with `run` and `combine`: give the folders or files of the ENDF neutron files.
`run` finds the incident neutron files (NSUB=10) in all the folders and ignores the other files
(thermal scattering, photo-atomic, ...). `--only Fe56 Co59` processes only some nuclides.
Each nuclide must be there only one time, because its output file has the name of the nuclide
(`Fe56.asc`, `Ag110m.asc`).

`combine` makes the NIS cache too, with `build/bin/nis_cache` (build NIS first, or give `--nis-cache`).

## The steps for one nuclide

These are the steps of the UKAEA processing for FISPACT-II (the `dpp` prototype), with the same inputs:

| Step | Program | What it does |
|---|---|---|
| 1 | NJOY: RECONR, BROADR, UNRESR, HEATR, GASPR | PENDF at 293.6 K, with the kerma and damage (MT 301-447) and the gas production (MT 203-207) |
| 2 | LINEAR, RECENT, SIGMA1 | linear cross sections, resonances, Doppler broadening at 293.6 K |
| 3 | SIXPAK, LINEAR | isomer production (MF9) from the MF6 data, if the file has it |
| 4 | MERGER, ACTIVATE | MF9 into MF10 (activation cross sections for each isomer) |
| 5 | MERGER (4 times) | add the NJOY heating, damage and gas data, keep MF1-3, 8, 10, 33 and 40 |
| 6 | FIXUP, DICTIN | make the data consistent (sums), make the index (MF1/MT451) again |
| 7 | GROUPIE | average in the UKAEA 1102 groups |

The program checks that the output has the material and the 1102 group total cross section.

## Speed

- **Many cores.** `run` processes `-j` nuclides at the same time (the default is the number of cores). The
  largest files start first, so that no long nuclide runs alone at the end.
- **NJOY at the same time as PREPRO.** For each nuclide, NJOY runs at the same time as steps 2 to 4, because
  they do not need each other. NJOY is the longest step.
- **Files in memory.** The work folders are in `/dev/shm` (a memory file system) if there is one, because
  PREPRO reads and writes many scratch files. Give another place with `--work-dir`.
- **PREPRO with -O2.** `tools` compiles PREPRO with `-O2` (the PREPRO makefile uses `-O0`). This is about 2 times
  faster, and the output is the same (checked byte for byte for the test nuclides).
- **Done is not done again.** A nuclide with an output file is skipped, so a stopped run continues where it
  stopped. `--force` processes all again.

For example, 28 files from the NJOY2016 tests (TENDL, ENDF/B-VIII.0, from H1 to Pu245) take 41 s on 4 cores.

## Failures

A nuclide can fail if PREPRO or NJOY cannot process its data. The run goes on with the other nuclides. At
the end, the command gives the failed nuclides and the first error message of each. The full output of all
the steps is in `<output>/logs/<nuclide>.log`, and `<output>/summary.json` has the result and the time of each
step for each nuclide. `--keep` keeps the work folders, so you can run a step again by hand.

## Notes

- `fetch` downloads `https://tendl.imperial.ac.uk/tendl_<year>/tar_files/TENDL-n.tgz` (TENDL moved there from
  tendl.web.psi.ch). If the site moves again, give the address with `--url`, or download and unpack the files yourself and use `run`.
- `tools` takes the latest NJOY2016 and PREPRO. Give a tag with `--njoy-ref` and `--prepro-ref` to fix the
  versions. `--njoy` and `--prepro-dir` (or `NJOY_EXE` and `PREPRO_EXE_DIR`) use programs that you already have.
- The tests are in [processingtests.py](../src/tests/python/processingtests.py). The tests that run NJOY and PREPRO
  need `tools` first (they use two files from the NJOY2016 tests).
