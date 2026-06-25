# CAFA6 winners → DeepGO-PlusPlus: what would actually move our number

Cross-cutting synthesis of the 1st–6th place writeups (+ figures). The question this
answers: **beyond what DG++/`cafa-baseline` already does, what work would improve the
classifier — and by roughly how much?**

## 1. Which signal appears in which solution

| signal | 1st | 2nd | 3rd | 4th | 5th | 6th | DG++ has it? |
|---|:--:|:--:|:--:|:--:|:--:|:--:|---|
| Homology kNN (BLAST/DIAMOND) | ✅ | — | — | — | ✅ | — | ✅ `diam` |
| PLM head (ESM2 / ProtT5) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ `mlp`/`prostt5`/`esm2_3b`/`esm2_head` |
| Structure (FoldSeek / 3Di) | ✅ | — | — | ✅ | ✅ | ✅ | ✅ `foldseek` (GPU-gated) |
| InterPro / domain LR | ✅ | — | — | — | — | — | ✅ `interpro`/`interpro_lr` |
| **PPI / Net-KNN (STRING)** | ✅ | — | — | — | ✅ BP-only | — | ✅ `net`/`net_union` (under review) |
| **Literature / PubMed text** | ✅ | ✅ | — | — | ✅ | ✅ | ⚠️ `lit` (BM25) — weaker form |
| **NEA: electronic GOA + evidence as features** | (via GOA) | ✅ 38-feat | — | — | ✅ top-1000 | ✅ encoder | ❌ **missing** |
| **Taxon encoding (one-hot / lineage)** | ✅ 21-dim | ✅ ~30 | ✅ | — | ✅ hier. | ✅ lineage | ❌ **missing in integrator** |
| Cross-ontology features (other-aspect labels) | — | ✅ SVD | — | — | (NEA) | — | ❌ missing |
| Conditional targets / MCM hierarchy | — | ✅ | ✅ | — | ✅ MCM | — | ✅ `esm2_head` MCM (C-HMCNN) |
| Learned integrator (LTR / GCN / weighted) | ✅ XGB-LTR | ✅ GCN | ✅ GCN | ⚠️ blend | ✅ wmean | ✅ NN | ✅ per-aspect logreg integrator |
| Taxon-constraint / consistency postproc | — | ✅ | — | — | ✅ | — | ✅ gspa-core **SAT4J** (not wired into DG++) |
| Test-distribution-matched validation | ✅ | — | — | — | ✅ group-CV | — | ✅ clean-A/B no-knowledge recon |

## 2. The unambiguous consensus
- **Integration is the lever, not new architectures.** 3rd/4th are CAFA5 redux; 1st's
  +0.12 over its best single component came from LTR. Our `cafa-baseline` (0.359→0.483
  on no-knowledge) is the same lever; the remaining gap to 0.524 is mostly **signals we
  don't have yet + a pre-t0 training population for the integrator**.
- **Three signals every strong team used that we lack: NEA, literature embeddings,
  taxon encoding.** These are the highest-EV additions.
- **Net/PPI is the weakest, anti-generalizing component** (GOAlpha Fig 2: 0.45 vs
  0.56–0.59; Fig 3: the *only* component that drops public→private, 0.29→0.26; 5th:
  BP-only). Our own clean-eval agrees (`diam` dominates `net`). → keep net, but **BP-only
  and low-weight**; do not expect it to be a win.

## 3. Prioritised backlog (highest expected value first)

### P1 — NEA: non-experimental GOA annotations + evidence codes as features ⭐
The single most consistent missing signal (2nd, 5th, 6th; root = its7171 CAFA5-3rd).
~30 % of electronic labels later become experimental. **CPU-only, no GPU, cheap.**
- Build: download `goa_uniprot_all.gaf` (pin a **pre-t0** release — 226/228), keep
  non-experimental codes (IEA/IBA/ISO/ISS/ND/…), propagate, take top-1000 frequent terms.
- Two consumers: (a) flat per-GO assignment vector → straight into the logreg integrator
  (cheap, do first); (b) optional term-evidence encoder for a NN head (later).
- **Leakage discipline**: features must be the *electronic* annotations available before
  t0; never the experimental labels we score. Validate on clean-B (novel) proteins.
- Expected: this is what carried 2nd over 3rd (≤0.04, "70–80 % from articles" — but NEA
  is the other recurring lever); realistically **+0.01–0.03 f_w** on limited/partial,
  smaller on pure no-knowledge (novel proteins have fewer electronic annotations too).

### P2 — Taxon encoding into the integrator ⭐
Every supervised solution used it; ≈ +0.02–0.03 on base models (2nd, CAFA5 est.).
- Add a taxon feature to `cafa-baseline`/Light integrator: one-hot of ~30 best-represented
  taxa + `[UNK]`, or a Phylum/Order-level bucket (5th). We already have taxon per protein
  in `net_index.tsv`. Trivial, CPU-only.
- Expected: **+0.01–0.02 f_w**, free.

### P3 — Upgrade `lit` from BM25 to a PubMed-embedding component
Literature was 2nd place's biggest new win and is in 1st/5th/6th. We have a BM25 `lit`;
the top teams use learned text embeddings (MedCPT, PubMedBERT) or TF-IDF⊕PLM.
- Cheapest upgrade: TF-IDF(title+abstract) like 2nd place (CPU). Better: MedCPT/PubMedBERT
  embedding (GPU once, then a cheap head).
- **Caveat from the figures**: literature **tops public but drops on private** (GOXML
  0.48→0.42) — it helps *known* proteins, hurts on *novel*. Keep it **down-weighted**, and
  strictly pre-t0 text only. Our clean-eval already shows `lit` is a top component on the
  current (limited-knowledge-ish) GT — re-check it on clean-B before trusting it.
- Expected: **+0.01–0.03** on limited/partial; ~0 or negative on pure no-knowledge.

### P4 — Wire gspa-core SAT4J taxon-constraint + consistency into DG++ postprocessing
Our genuine differentiator. 2nd (avg term/maxchild/minparent) and 5th (forbidden-GO drop,
drop-known-labels, max-propagation) do cheap versions; we have a real SAT consistency
checker. Apply as a post-pass: drop NEVER_IN/ONLY_IN taxon violations, enforce true-path.
- Expected: small but free precision gain (+0.005–0.01), and a paper-worthy distinctive.
- Also add the **cheap wins** unconditionally: drop (protein,GO) already in `train_terms`;
  top-N per aspect; score≥0.01. (Check `cafa-baseline` does these.)

### P5 — kNN scoring fix (IC-weighting + Σweight normalization)
Our `net`/`diam`/`foldseek` kNN votes weight by similarity and normalize by per-protein
**max**. 5th place weights each vote by term **IC** and normalizes by **Σ weight**
(weighted mean). Cheap to A/B-test on the existing components.
- Expected: **+0.005–0.015** on the kNN components, possibly more for `diam`.

### P6 — Cross-ontology features (limited/partial knowledge only)
2nd place: SVD-512 of known MF/CC one-hot labels → feature to predict BP. Only helps the
limited/partial test sets (where the protein already has some labels). Lower priority since
our headline target is no-knowledge.

### P7 — Integrator training population (already in the plan, not a writeup item)
Train the integrator on the **pre-t0 population** with component scores re-run on the train
FASTA, so a richer XGBoost-LTR (1st place) / GCN (2nd/3rd) can help without the
in-distribution leakage we saw. This is the structural gap to 0.524, independent of new signals.

## 4. What NOT to chase
- **Bigger/more PLMs** — diminishing returns confirmed both ways (our Phase-B: ESM2-3B
  redundant with mlp; only ProstT5 complementary; 2nd place: fine-tuning didn't help).
- **More training data** — 4th place: 150k experimental-only beat 550k/1.88M. Quality>quantity.
- **IA-weighted *loss* as a default** — contradictory across teams (4th hurt, 5th `1+IA`,
  6th `sqrt(IA)`); keep IA for the *metric*, treat loss-weighting as a per-aspect ablation.
- **Net/PPI as a headline signal** — weakest, anti-generalizing; BP-only, low weight.

## 5. Bottom line
The path to ~0.50→0.52 on our faithful no-knowledge recon is **add NEA + taxon to the
learned integrator (P1+P2, CPU, days), upgrade literature carefully (P3), then a pre-t0
LTR integrator (P7)** — plus our SAT consistency pass (P4) as the distinctive. New PLMs and
more data are *not* the lever. This matches GOAlpha's own conclusion: homology+PLM recover
most of the score; structure, **literature, and electronic-annotation evidence** add the
complementary signal that separates the medal band — through **integration**, not architecture.
