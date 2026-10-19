# Regenerate manuscript tables, figures and number checks from results/*.tsv
PY ?= python3
REPO := $(abspath ..)

.PHONY: all tables figures check results
all: tables figures check

tables:
	$(PY) tables/make_tables.py --repo $(REPO)

figures:
	cd $(REPO)/paper && $(PY) $(REPO)/paper-reproducibility/figures/generate_figures.py

check:
	$(PY) checks/check_numbers.py --repo $(REPO)

results:
	@echo "Run the pipeline folders under pipelines/ on their listed hosts; see README.md"
