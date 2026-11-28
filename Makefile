# NIS: get the data, process it and compare the benchmark results, in 4 steps.
#
#   make setup     get the C++ dependencies, build NIS (CMake) and pynis (uv)
#   make data      get the decay data, and make the reaction data of each library (NJOY and PREPRO, hours)
#   make bench     run the benchmark for each reaction library, all with the same decay data
#   make compare   the C/E statistics of all the runs next to each other, and a summary plot of each benchmark
#
# Or all of it: make all. Choose the data on the command line, for example:
#
#   make compare REACTIONS="jeff4 tendl2025" DECAY=jeff4
#   make compare REACTIONS=jeff4 DECAY=jeff311
#
# REACTIONS: jeff4, or tendl<year> (tendl2025). DECAY: jeff311, jeff4, endfbvii1, endfbviii0
# (python3 processing/process.py libraries). The results are in results/:
# <benchmark>-<reactions>+<decay>.csv for each run, and <benchmark>-stats.txt and <benchmark>-summary.png for the
# comparison (<benchmark>-<quantity>-summary.png for a benchmark with more than one quantity).
# A run is done again only when its data file changes. "make clean-results" removes the results.
#
# A reaction data file of an older version of the processing is made again (hours). ALLOW_OLD=1 uses it as it is:
#   make bench ALLOW_OLD=1
#
# BENCHMARK chooses the benchmark, FNS if not given (the names are the folders in benchmarks/):
#   make compare BENCHMARK=irdff2               the IRDFF-II cross sections
#   make compare BENCHMARK=exfor14              the EXFOR 14 MeV cross sections (its summary plot is large, 1162 cases)
#   make compare BENCHMARK=fission_decay_heat   the fission decay heat, with the fission yields of YIELDS (= DECAY)
#   make compare BENCHMARK=all                  all of them
#
# Other targets: make test-data (the files the tests need), make test (the Python and C++ tests).

REACTIONS ?= jeff4 tendl2025
DECAY ?= jeff4
BENCHMARK ?= fns
# the fission yields of the fission benchmark: jeff311 or jeff4 (python3 processing/process.py yields <library>)
YIELDS ?= $(DECAY)
RESULTS ?= results
# ALLOW_OLD=1: use a reaction data file of an older version of the processing as it is (see processing/README.md)
ALLOW_OLD ?=
JOBS ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

PYTHON ?= python3
UV ?= uv
PROCESS := $(PYTHON) processing/process.py
ALLOW_OLD_FLAG := $(if $(ALLOW_OLD),--allow-old,)
BENCH := $(UV) run python scripts/bench.py
NJOY := processing/tools/bin/rustjoy
NIS := build/bin/nis
NIS_CACHE := build/bin/nis_cache

# the benchmarks of BENCHMARK=all, and the quantities of each (one summary plot for each quantity)
ALL_BENCHMARKS := fns irdff2 exfor14 fission_decay_heat
QUANTITIES_fns := decay_heat
QUANTITIES_irdff2 := cross_section
QUANTITIES_exfor14 := cross_section
QUANTITIES_fission_decay_heat := decay_heat beta_heat gamma_heat
BENCHMARKS := $(if $(filter all,$(BENCHMARK)),$(ALL_BENCHMARKS),$(BENCHMARK))

DECAY_FILE := data/$(DECAY).dec
REACTION_FILES := $(foreach r,$(REACTIONS),data/$(r).reac)
# the results file of a run: $(call run_file,<benchmark>,<reactions>)
run_file = $(RESULTS)/$(1)-$(2)+$(DECAY).csv
RUNS := $(foreach b,$(BENCHMARKS),$(foreach r,$(REACTIONS),$(call run_file,$(b),$(r))))
# the fission benchmark also needs the fission yields
fission_yields = $(if $(filter fission_decay_heat,$(1)),data/$(YIELDS).nfy,)

.PHONY: all setup data bench compare test-data test clean-results FORCE

all: compare

# ---- 1. setup ----

setup: dependencies/ENDFtk $(NIS)
	$(UV) sync

dependencies/ENDFtk:
	bash docker/setup_local_deps.sh

$(NIS): FORCE dependencies/ENDFtk
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
	cmake --build build -j $(JOBS)

# ---- 2. data ----

data: $(DECAY_FILE) $(REACTION_FILES)

# process.py does nothing if the file is there (and makes it again if it is an old git lfs pointer file),
# so it is called each time, and a run is done again only when the file really changes
data/%.dec: FORCE
	$(PROCESS) decay $*

data/%.nfy data/%.sfy: FORCE
	$(PROCESS) yields $*

$(NJOY):
	$(PROCESS) tools -j $(JOBS)

data/tendl%.reac: FORCE | $(NJOY)
	$(PROCESS) tendl --year $* -j $(JOBS) $(ALLOW_OLD_FLAG)

data/%.reac: FORCE | $(NJOY)
	$(PROCESS) reactions $* -j $(JOBS) $(ALLOW_OLD_FLAG)

# ---- 3. bench ----

bench: $(RUNS)

# The reaction cache (data/<reactions>.nisreac) is made again, with nis_cache, if it is older than the .reac file
# or of another cache version. NIS is built first (if it is up to date, this is fast), so that nis_cache writes
# the cache version that the current code reads.
# the label is the data, so the runs can not be mixed up in the comparison. $(1): benchmark, $(2): reactions
define RUN_RULE
$(call run_file,$(1),$(2)): data/$(2).reac $(DECAY_FILE) $(call fission_yields,$(1)) | $(NIS)
	@mkdir -p $(RESULTS)
	$(BENCH) run $(1) --label $(2)+$(DECAY) --decay-data $(DECAY_FILE) --reaction-data data/$(2).reac \
	    --nis-cache $(NIS_CACHE) $(if $(call fission_yields,$(1)),--fission-yields $(call fission_yields,$(1)),) --jobs $(JOBS) --output $$@
endef
$(foreach b,$(BENCHMARKS),$(foreach r,$(REACTIONS),$(eval $(call RUN_RULE,$(b),$(r)))))

# ---- 4. compare ----

compare: $(foreach b,$(BENCHMARKS),compare-$(b))

# the summary plot of a benchmark and quantity: $(call summary_file,<benchmark>,<quantity>)
summary_file = $(RESULTS)/$(1)$(if $(word 2,$(QUANTITIES_$(1))),-$(2),)-summary.png

define COMPARE_RULE
.PHONY: compare-$(1)
compare-$(1): $(foreach r,$(REACTIONS),$(call run_file,$(1),$(r)))
	$(BENCH) stats $$^ --output $(RESULTS)/$(1)-stats.txt
	$(foreach q,$(QUANTITIES_$(1)),$(BENCH) plot $$^ --benchmark $(1) --quantity $(q) --summary $(call summary_file,$(1),$(q))
	)@echo "results: $(RESULTS)/$(1)-stats.txt and $(foreach q,$(QUANTITIES_$(1)),$(call summary_file,$(1),$(q)))"
endef
$(foreach b,$(BENCHMARKS),$(eval $(call COMPARE_RULE,$(b))))

# ---- tests ----

test-data: data/jeff311.dec data/endfbvii1.dec data/jeff311.nfy data/jeff311.sfy data/jeff4.reac

test: setup
	./build/bin/nistests
	$(UV) run pytest src/tests/python/testsuite.py

clean-results:
	rm -rf $(RESULTS)

FORCE:
