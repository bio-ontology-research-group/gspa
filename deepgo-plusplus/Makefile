# DeepGO-PlusPlus — reproducible retraining pipeline
# ===================================================
# Rebuild the frozen integrator (and the CPU-buildable net/lit components) from
# a UniProt + STRING + CAFA release. Run `make help` for the target list.
#
# The heavy, GPU-bound components (diam/foldseek/clean/interpro/mlp/prostt5) are
# produced upstream by GSPA's predictors / PLM heads (see README §"Components")
# and dropped into $(COMPONENTS). This Makefile builds the CPU components
# (net, lit), assembles the ground truth, and trains + evaluates the integrator
# — the part that must be re-runnable at every UniProt/STRING release.
#
# Configure inputs by copying config.mk.example -> config.mk and editing paths.
# Every input release should be recorded in VERSIONS.md for provenance.

-include config.mk

# ---- tool / interpreter -----------------------------------------------------
# Default to `uv run` so the pinned environment in pyproject.toml is used.
PY        ?= uv run python
CAFAEVAL  ?= uv run cafaeval

# ---- inputs (override in config.mk) ----------------------------------------
UNIPROT_SPROT_DAT   ?= $(HOME)/Public/software/cafa6/uniprot_sprot.dat.gz
TESTSUPERSET_ANNOTS ?= testsuperset_exp_annots.tsv
TRAIN_TERMS         ?= $(HOME)/Public/software/cafa6/kaggle-official/train_terms.tsv
IA                  ?= $(HOME)/Public/software/cafa6/kaggle-official/IA.tsv
DAG                 ?= $(HOME)/Public/software/cafa6/go-dag.tsv
OBO                 ?= $(HOME)/Public/software/cafa6/go.obo
TAXON               ?= $(HOME)/Public/software/cafa6/kaggle-official/testsuperset-taxon-list.tsv
STRING_DIR          ?= string                # dir of {taxid}.protein.links.v12.0.txt.gz
T0                  ?= 20260202              # CAFA6 t0 (final submission deadline)

# ---- pipeline knobs --------------------------------------------------------
# Components fed to the integrator. The shipped default is the 6 GPU/homology
# components; append ,net (recommended) and/or ,lit once those are built.
COMPONENT_LIST ?= diam,foldseek,clean,interpro,mlp,prostt5
NET_MIN_CONF   ?= 400
NET_TOPK       ?= 50
LIT_TOPK       ?= 30

# ---- layout ----------------------------------------------------------------
PIPE       := pipeline
WORK       ?= build
COMPONENTS ?= $(WORK)/components
GT         := $(WORK)/gt
MODELS     := models
MODEL      ?= $(MODELS)/deepgo_plusplus_integrator.json

.PHONY: help all components integrator eval test clean realclean check-inputs

help:
	@echo "DeepGO-PlusPlus reproducible retraining"
	@echo ""
	@echo "  make index          build text_string_index.tsv + net_index.tsv (UniProt pass)"
	@echo "  make gt             build CAFA6 ground truth (gt_{no,limited,partial,all}.tsv)"
	@echo "  make net            build the Net-KNN component (STRING PPI)         -> $(COMPONENTS)/net.tsv.gz"
	@echo "  make lit            build the BM25 literature component              -> $(COMPONENTS)/lit.tsv.gz"
	@echo "  make integrator     train + freeze the per-aspect logreg stacker     -> $(MODEL)"
	@echo "  make eval           score the OOF predictions with official cafaeval"
	@echo "  make all            index -> gt -> net -> integrator"
	@echo "  make test           run the pytest regression suite (no GPU/network)"
	@echo "  make clean          remove \$(WORK) (keeps models/)"
	@echo ""
	@echo "Configure inputs in config.mk (see config.mk.example). COMPONENT_LIST=$(COMPONENT_LIST)"
	@echo "Upstream GPU components (diam/foldseek/clean/interpro/mlp/prostt5) must be"
	@echo "placed under $(COMPONENTS)/ as <component>.tsv[.gz] — see README §Components."

# --- one streaming UniProt pass: text + STRING xref + taxon index ------------
$(WORK)/text_string_index.tsv: $(PIPE)/build_text_string_index.py | $(WORK)
	$(PY) $(PIPE)/build_text_string_index.py $(UNIPROT_SPROT_DAT) $@

$(WORK)/net_index.tsv: $(WORK)/text_string_index.tsv
	cut -f1-3 $< > $@

index: $(WORK)/net_index.tsv

# --- CAFA6 ground truth from a dated GOA filter ------------------------------
$(GT)/gt_no.tsv: $(PIPE)/build_groundtruth.py $(TESTSUPERSET_ANNOTS) | $(WORK)
	$(PY) $(PIPE)/build_groundtruth.py --annots $(TESTSUPERSET_ANNOTS) --t0 $(T0) --outdir $(GT)

gt: $(GT)/gt_no.tsv

# --- Net-KNN component over STRING (CPU) -------------------------------------
$(COMPONENTS)/net.tsv.gz: $(PIPE)/build_net_component.py $(WORK)/net_index.tsv $(GT)/gt_no.tsv | $(COMPONENTS)
	$(PY) $(PIPE)/build_net_component.py \
	  --index $(WORK)/net_index.tsv --train-terms $(TRAIN_TERMS) \
	  --queries $(GT)/test_proteins.txt --string-dir $(STRING_DIR) \
	  --out $(COMPONENTS)/net.tsv --min-conf $(NET_MIN_CONF) --topk $(NET_TOPK)
	gzip -f $(COMPONENTS)/net.tsv

net: $(COMPONENTS)/net.tsv.gz

# --- BM25 literature component (CPU) ----------------------------------------
$(COMPONENTS)/lit.tsv.gz: $(PIPE)/build_lit_component.py $(WORK)/text_string_index.tsv $(GT)/gt_no.tsv | $(COMPONENTS)
	$(PY) $(PIPE)/build_lit_component.py \
	  --index $(WORK)/text_string_index.tsv --train-terms $(TRAIN_TERMS) \
	  --queries $(GT)/test_proteins.txt --out $(COMPONENTS)/lit.tsv --topk $(LIT_TOPK)
	gzip -f $(COMPONENTS)/lit.tsv

lit: $(COMPONENTS)/lit.tsv.gz

# --- train + freeze the integrator ------------------------------------------
# Requires every component named in COMPONENT_LIST to exist under $(COMPONENTS).
integrator: $(GT)/gt_no.tsv
	$(PY) $(PIPE)/train_integrator.py \
	  --components $(COMPONENTS) --gt $(GT)/gt_no.tsv \
	  --dag $(DAG) --ia $(IA) --train-terms $(TRAIN_TERMS) --taxon $(TAXON) \
	  --model logreg --features scores --component-list $(COMPONENT_LIST) \
	  --save-model $(MODEL) --out $(WORK)/ltr_run
	@echo "froze integrator -> $(MODEL)  (OOF preds: $(WORK)/ltr_run/preds/ltr.tsv)"

# --- evaluate the OOF predictions with the official cafaeval -----------------
eval: $(WORK)/ltr_run/preds/ltr.tsv
	$(CAFAEVAL) $(OBO) $(WORK)/ltr_run/preds $(GT)/gt_no.tsv \
	  -ia $(IA) -norm cafa -prop max -th_step 0.01 -out_dir $(WORK)/ltr_run/results
	@$(PY) -c "import sys; \
rows=open('$(WORK)/ltr_run/results/evaluation_best_f_w.tsv').read().splitlines(); \
h=rows[0].split('\t'); fi=h.index('f_w'); ni=h.index('ns'); \
m={r.split('\t')[ni]:float(r.split('\t')[fi]) for r in rows[1:]}; \
mf=m.get('molecular_function',0); bp=m.get('biological_process',0); cc=m.get('cellular_component',0); \
print(f'no-knowledge f_w  MF={mf:.3f} BP={bp:.3f} CC={cc:.3f}  mean={(mf+bp+cc)/3:.4f}')"

all: index gt net integrator

test:
	uv run pytest -q tests

$(WORK) $(COMPONENTS):
	mkdir -p $@

check-inputs:
	@for f in "$(UNIPROT_SPROT_DAT)" "$(TRAIN_TERMS)" "$(IA)" "$(DAG)" "$(OBO)" "$(TAXON)"; do \
	  [ -e "$$f" ] && echo "  ok   $$f" || echo "  MISS $$f"; done

clean:
	rm -rf $(WORK)

realclean: clean
	@echo "models/ retained on purpose; remove by hand if you really mean it"
