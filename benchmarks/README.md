# NIS benchmarks

All the benchmarks use one standard format, for all types of benchmark (fusion, fission and other). One set of
tools reads it: [scripts/bench.py](../scripts/bench.py) (the code is in [scripts/nisbench](../scripts/nisbench)).

| Benchmark | Type | Cases | What is measured |
|---|---|---|---|
| [fns](fns) | fusion | 132 | decay heat of foils after 14 MeV neutron irradiation (FNS, JAEA) |
| [fission_decay_heat](fission_decay_heat) | fission | 44 | beta, gamma and total decay heat after the fission of U233, U235, U238, Pu239, Pu241 and Th232: pulses and irradiations of 10 s to 1e5 s (IAEA CoNDERC) |

## Folders

```
benchmarks/
    group_structures/<name>.txt     energy group bounds (eV, ascending), e.g. CCFE-709
    <benchmark>/
        benchmark.json              what the benchmark is: title, type, references, where the files come from
        <case>.json                 one case: the source of truth (the case name can have folders, e.g. Co/2000exp_5min)
        <case>.inp                  the NIS input deck, made from <case>.json
        results/                    the NIS reference results and plots
```

**Rules**

1. The case file (`.json`) is the source of truth. Do not edit a deck (`.inp`): edit the case file, then make
   the decks again with `scripts/bench.py decks`. `scripts/bench.py check` (and a test) finds a deck that is not
   made from its case file.
2. A converter makes the case files from the original files of the benchmark (`scripts/bench.py convert`), so the
   conversion can be done again and checked. The original files are not in the repository. `benchmark.json` says
   what they are and where they come from.
3. A deck is a normal NIS input. You can run it with the `nis` program without the tools:
   `./build/bin/nis -d data/jeff311.dec -r data/jeff4.nisreac -i benchmarks/fns/Co/2000exp_5min.inp`

## benchmark.json

```json
{
  "format": "nis-benchmark",
  "version": 1,
  "name": "fns",
  "title": "FNS decay heat",
  "category": "fusion",
  "description": "What is irradiated and measured.",
  "references": ["The papers or reports of the experiment"],
  "source": "Where the original files come from, and how they were converted."
}
```

| Key | |
|---|---|
| `name` | the folder name |
| `category` | `fusion`, `fission` or `other` |
| `title`, `description`, `source` | text, must not be empty |
| `references` | a list of text, can be empty |

## Case file (`<case>.json`)

```json
{
  "format": "nis-benchmark-case",
  "version": 1,
  "benchmark": "fns",
  "name": "Co/2000exp_5min",
  "title": "Co, 5 minute irradiation, 2000 experiment",
  "sample": {
    "mass_kg": 0.001,
    "atoms": {"Co59": 1.02186e+22}
  },
  "steps": {
    "dt_s":        [300.0,     34.98, 15.0, 15.0],
    "flux_n_m2_s": [1.116e+14, 0.0,   0.0,  0.0]
  },
  "spectrum": {"group_structure": "CCFE-709", "values": [0.0, 0.0, "... 709 values"]},
  "measurements": [
    {
      "quantity": "decay_heat",
      "unit": "uW/g",
      "steps": [1, 2, 3],
      "values": [0.02729, 0.03019, 0.03046],
      "uncertainties": [0.001621, 0.001793, 0.001791],
      "calculations": {"FISPACT-II (TENDL-2017)": [null, 0.03285, 0.03277]}
    }
  ]
}
```

| Key | |
|---|---|
| `name` | the path of the file in the benchmark folder, without `.json` |
| `sample.atoms` | the atoms of each nuclide at the start (`Co59`, `Ag110m`, ...) |
| `sample.mass_kg` | the sample mass. Needed for a unit per mass (`uW/g`, `W/kg`, ...) |
| `fissions` | the number of fissions. Needed for a unit per fission (`MeV/fission`, ...). Not given with a `fission_source` |
| `fission_source` | optional, in place of a flux: `nuclide` (e.g. `U235`), `neutron_energy_ev` (the energy of the fission yields) and `fissions_at_start` (a pulse, default 0). Each fission adds the independent fission yields of the nuclide. With it, `sample.atoms` can be empty, and `fissions` is the fissions of the source |
| `steps.fission_rate_per_s` | optional, with a `fission_source`: the fission rate of each step (0 for cooling) |
| `steps` | the time steps in order: the length (s), and the total neutron flux (n/m²/s, 0 for cooling). The example shows only the first steps |
| `spectrum` | the flux shape: `values` (one for each group, any scale), with `group_structure` (a file in `group_structures/`) or `energy_bounds_ev` (ascending, one more than the values). `null` if no step has a flux. For each step, NIS scales the shape so that its total is the flux of the step |
| `measurements` | one or more measured quantities |
| `measurements[].steps` | the step (0 is the first) at the end of which each value is measured, going up |
| `measurements[].values`, `uncertainties` | one for each of `steps`, in `unit`. The uncertainty is one standard deviation |
| `measurements[].calculations` | optional: the results of other codes, one value (or `null`) for each of `steps` |

The **cooling time** of a value is the time from the end of the last step with a flux to the end of its step
(from the start, if no step has a flux).

**Quantities and units**

| Quantity | NIS result | Units |
|---|---|---|
| `decay_heat` | total decay heat | `W`, `MeV/s`, `W/kg`, `kW/kg`, `mW/g`, `uW/g`, `W/fission`, `MeV/s/fission`, `MeV/fission`, `MeV/(fission/s)` |
| `beta_heat` | decay heat of beta particles and electrons | the same |
| `gamma_heat` | decay heat of gamma and X rays | the same |
| `alpha_heat` | decay heat of alpha particles | the same |
| `activity` | total activity | `Bq`, `Bq/kg`, `Bq/g`, `Bq/fission` |

`MeV/fission` is the decay heat per fission times the cooling time, t·P(t), which is the usual way to give the decay
heat after a fission pulse. `MeV/(fission/s)` is the decay heat over the fission rate of the last step with a
fission rate, the usual way to give it after an irradiation of some length.

A **fission** case can be given in two ways. With a `fission_source` (as in `fission_decay_heat`): no cross sections
are needed, only the fission yields; a pulse is `fissions_at_start: 1`, an irradiation is a step with a fission rate.
Or as a real irradiation: the fissile atoms in `sample.atoms`, a spectrum, an irradiation step with a flux, then the
cooling steps, and `fissions` for the units per fission. Run both with fission yields (`--fission-yields`).
A case with a `fission_source` has no deck (the `nis` program has no fission source), so run it with `pynis`.

`scripts/bench.py pulse U235 --decay-data data/jeff311.dec --fission-yields data/jeff311.nfy` gives the decay heat
of an ideal pulse (the yields of one fission at t = 0, then decay) in MeV/s per fission, and the delayed beta and
gamma energy over all time. Use it to check a pulse case: for cooling times much longer than the irradiation, the
two must agree.

## Results file

`scripts/bench.py run` writes one CSV for all the benchmarks, one row for each measured value:

```
benchmark,case,quantity,unit,time_s,experiment,uncertainty,nis,FISPACT-II (TENDL-2017)
fns,Co/2000exp_5min,decay_heat,uW/g,34.98,0.02729,0.001621,2.832815e-02,
fns,Co/2000exp_5min,decay_heat,uW/g,49.98,0.03019,0.001793,2.829702e-02,3.285241e-02
```

`time_s` is the cooling time. The columns after `nis` are the `calculations` of the cases (empty where a case has
no value). `stats` and `plot` read this file.

## The commands

```bash
uv run python scripts/bench.py list -v                        # the benchmarks and their cases
uv run python scripts/bench.py check                          # the case files are correct, the decks are up to date
uv run python scripts/bench.py run fns --decay-data data/jeff311.dec --reaction-data data/jeff4.reac \
    --output results.csv                                      # all the cases (or --cases Co Fe/1996exp_7hour)
uv run python scripts/bench.py stats results.csv              # C/E for NIS and each calculation
uv run python scripts/bench.py plot results.csv --summary summary.png --cases-dir plots
uv run python scripts/bench.py decks                          # make all the decks again from the case files
uv run python scripts/bench.py convert fns --source path/to/originals   # make a benchmark from its original files
```

`run` uses `pynis` by default. `--engine nis --nis build/bin/nis` runs each deck with the `nis` program in place
of `pynis` (the results are the same). `--fission-yields` and `--sf-yields` give the fission yield files.

## Add a benchmark

1. Write a converter `scripts/nisbench/convert/<name>.py` with `NAME`, `METADATA` (the content of
   `benchmark.json`), `SOURCE_HELP` and `convert(source, only)`, which returns the cases.
   [convert/fns.py](../scripts/nisbench/convert/fns.py) is an example. If the spectrum uses a new group structure,
   add its bounds to `group_structures/`.
2. Make the files: `uv run python scripts/bench.py convert <name> --source path/to/originals`. The tools check each
   case file when they write it.
3. Run it, and keep the reference results and the summary plot in `<name>/results/`.
4. Add it to the table at the top of this file.
