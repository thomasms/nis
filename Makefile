# NIS: get the data, process it and compare the benchmark results, in 4 steps.
#
#   make setup     get the C++ dependencies, build NIS (CMake) and pynis (uv)
#   make data      get the decay data, and make the reaction data of each library (NJOY and PREPRO, hours)
#   make bench     run the benchmark for each reaction library, all with the same decay data
#   make compare   the C/E statistics of all the runs next to each other, and a summary plot
#
# Or all of it: make all. Choose the data on the command line, for example:
#
#   make compare REACTIONS="jeff4 tendl2025" DECAY=jeff4
#   make compare REACTIONS=jeff4 DECAY=jeff311
#
# REACTIONS: jeff4, or tendl<year> (tendl2025). DECAY: jeff311, jeff4, endfbvii1, endfbviii0
# (python3 processing/process.py libraries). The results are in results/:
# <benchmark>-<reactions>+<decay>.csv for each run, and stats.txt and summary.png for the comparison.
# A run is done again only when its data file changes. "make clean-results" removes the results.
#
# A reaction data file of an older version of the processing is made again (hours). ALLOW_OLD=1 uses it as it is:
#   make bench ALLOW_OLD=1
#
# BENCHMARK=irdff2 runs the IRDFF-II cross sections in place of FNS: make compare BENCHMARK=irdff2
#
# Other targets: make test-data (the files the tests need), make test (the Python and C++ tests).

REACTIONS ?= jeff4 tendl2025
DECAY ?= jeff4
BENCHMARK ?= fns
RESULTS ?= results
# ALLOW_OLD=1: use a reaction data file of an older version of the processing as it is (see processing/README.md)
ALLOW_OLD ?=
JOBS ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

PYTHON ?= python3
UV ?= uv
PROCESS := $(PYTHON) processing/process.py
ALLOW_OLD_FLAG := $(if $(ALLOW_OLD),--allow-old,)
BENCH := $(UV) run python scripts/bench.py
NJOY := processing/tools/bin/njoy

DECAY_FILE := data/$(DECAY).dec
REACTION_FILES := $(foreach r,$(REACTIONS),data/$(r).reac)
RUNS := $(foreach r,$(REACTIONS),$(RESULTS)/$(BENCHMARK)-$(r)+$(DECAY).csv)

.PHONY: all setup data bench compare test-data test clean-results FORCE

all: compare

# ---- 1. setup ----

setup: dependencies/ENDFtk build/bin/nis
	$(UV) sync

dependencies/ENDFtk:
	bash docker/setup_local_deps.sh

build/bin/nis: FORCE dependencies/ENDFtk
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

# the label is the data, so the runs can not be mixed up in the comparison
define RUN_RULE
$(RESULTS)/$(BENCHMARK)-$(1)+$(DECAY).csv: data/$(1).reac $(DECAY_FILE)
	@mkdir -p $(RESULTS)
	$(BENCH) run $(BENCHMARK) --label $(1)+$(DECAY) --decay-data $(DECAY_FILE) --reaction-data data/$(1).reac \
	    --jobs $(JOBS) --output $$@
endef
$(foreach r,$(REACTIONS),$(eval $(call RUN_RULE,$(r))))

# ---- 4. compare ----

compare: $(RUNS)
	$(BENCH) stats $(RUNS) --output $(RESULTS)/stats.txt
	$(BENCH) plot $(RUNS) --benchmark $(BENCHMARK) --summary $(RESULTS)/summary.png
	@echo "results: $(RESULTS)/stats.txt and $(RESULTS)/summary.png"

# ---- tests ----

test-data: data/jeff311.dec data/endfbvii1.dec data/jeff311.nfy data/jeff311.sfy data/jeff4.reac

test: setup
	./build/bin/nistests
	$(UV) run pytest src/tests/python/testsuite.py

clean-results:
	rm -rf $(RESULTS)

FORCE:
