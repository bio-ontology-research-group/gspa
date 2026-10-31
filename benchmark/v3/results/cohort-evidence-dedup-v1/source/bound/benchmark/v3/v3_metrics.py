"""Reusable, dependency-free metrics for the GSPA v3 production benchmarks.

Everything here is a pure function of predictions and labels. Nothing in this
module produces predictions: it only scores the output of ``gspa integrate``.

Design: every metric is computed from an *additive sufficient-statistics
vector* per genome. Pooling genomes is summing vectors, and a genome-level
bootstrap is resampling the vectors, so confidence intervals cost
O(genomes x vector length) per replicate regardless of the number of
(protein, term) pairs.

Conventions (see docs/benchmarks/v3-protocol.md):

* Evaluation universe = proteins that carry at least one label of the kind
  being scored. Predictions for unlabelled proteins are ignored, because the
  absence of a label is not evidence of a negative.
* A prediction counts at threshold ``t`` when ``score >= t``.
* Fixed thresholds and the Fmax sweep use the grid 0.01 .. 0.99.
* ``fmax`` selects its threshold on the evaluated labels. It is a descriptive
  upper bound, not a held-out estimate; fixed-threshold metrics are the
  deployable ones.
"""
from __future__ import annotations

import math
import random
import re
from bisect import bisect_left
from collections import defaultdict
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

GRID: List[float] = [i / 100 for i in range(1, 100)]
N_BINS = 10
GO_ROOTS = frozenset({"GO:0008150", "GO:0003674", "GO:0005575"})

_COMPLETE_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

Pair = Tuple[str, str]


# --------------------------------------------------------------------- ids

def normalize_ec(value: str) -> str:
    v = value.strip()
    if v.upper().startswith("EC:"):
        v = v[3:]
    return v


def is_complete_ec(value: str) -> bool:
    """True for a four-level EC with numeric levels (no ``-`` wildcards, no ``n`` preliminary levels)."""
    return bool(_COMPLETE_EC.match(normalize_ec(value)))


def normalize_go(value: str) -> str:
    v = value.strip().replace("_", ":")
    return v


# ---------------------------------------------------------------- ontology

class Ontology:
    """Minimal GO view: is_a + part_of parents, alt ids, obsolete terms."""

    def __init__(self) -> None:
        self.parents: Dict[str, Set[str]] = defaultdict(set)
        self.alt: Dict[str, str] = {}
        self.obsolete: Set[str] = set()
        self._anc: Dict[str, frozenset] = {}

    def primary(self, term: str) -> str:
        return self.alt.get(term, term)

    def ancestors(self, term: str) -> frozenset:
        """The term itself plus every is_a / part_of ancestor."""
        term = self.primary(term)
        hit = self._anc.get(term)
        if hit is not None:
            return hit
        seen: Set[str] = set()
        stack = [term]
        while stack:
            t = stack.pop()
            if t in seen:
                continue
            seen.add(t)
            stack.extend(self.parents.get(t, ()))
        out = frozenset(seen)
        self._anc[term] = out
        return out

    def propagate_truth(self, pairs: Iterable[Pair]) -> Set[Pair]:
        out: Set[Pair] = set()
        for pid, term in pairs:
            term = self.primary(term)
            if term in self.obsolete:
                continue
            for a in self.ancestors(term):
                if a not in GO_ROOTS:
                    out.add((pid, a))
        return out

    def propagate_scores(self, pred: Mapping[Pair, float]) -> Dict[Pair, float]:
        """Ancestors inherit the maximum score of their descendants."""
        out: Dict[Pair, float] = {}
        for (pid, term), s in pred.items():
            term = self.primary(term)
            if term in self.obsolete:
                continue
            for a in self.ancestors(term):
                if a in GO_ROOTS:
                    continue
                k = (pid, a)
                if s > out.get(k, -1.0):
                    out[k] = s
        return out


def load_obo(path: str) -> Ontology:
    onto = Ontology()
    cur: Optional[str] = None
    in_term = False
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("["):
                in_term = line == "[Term]"
                cur = None
                continue
            if not in_term or ":" not in line:
                continue
            key, _, val = line.partition(":")
            val = val.split("!")[0].strip()
            if key == "id":
                cur = val
            elif cur is None:
                continue
            elif key == "alt_id":
                onto.alt[val] = cur
            elif key == "is_obsolete" and val == "true":
                onto.obsolete.add(cur)
            elif key == "is_a":
                onto.parents[cur].add(val.split()[0])
            elif key == "relationship":
                parts = val.split()
                if len(parts) >= 2 and parts[0] == "part_of":
                    onto.parents[cur].add(parts[1])
    return onto


# ---------------------------------------------------- sufficient statistics

# Vector layout of a pair-level GenomeStats vector.
_I_NPOS, _I_NPAIRS, _I_BRIER = 0, 1, 2
_I_TP = 3
_I_FP = _I_TP + len(GRID)
_I_BIN = _I_FP + len(GRID)          # N_BINS x (n, sum_p, sum_y)
PAIR_VEC_LEN = _I_BIN + 3 * N_BINS


def _grid_index(tau: float) -> int:
    i = round(tau * 100) - 1
    if not (0 <= i < len(GRID)) or abs(GRID[i] - tau) > 1e-9:
        raise ValueError(f"threshold {tau} is not on the 0.01..0.99 grid")
    return i


def candidate_universe(truth: Set[Pair], *preds: Mapping[Pair, float]) -> Set[Pair]:
    """Declared evaluation universe of one genome: every true pair plus every
    pair emitted by any compared arm, restricted to labelled proteins.

    Closed-world assumption: a pair in this universe without a label is
    treated as false. Labels are incomplete, so true-but-unlabelled pairs are
    counted as errors (this biases precision down, equally for all arms).
    """
    labelled = {p for p, _ in truth}
    out = set(truth)
    for pred in preds:
        out.update(k for k in pred if k[0] in labelled)
    return out


def pair_stats(pred: Mapping[Pair, float], truth: Set[Pair],
               universe: Optional[Set[Pair]] = None) -> List[float]:
    """Sufficient statistics of one genome for pair-level (protein, term) metrics.

    Calibration (Brier, ECE) is computed over the fixed ``universe`` with an
    unemitted pair scored 0, so omitting a prediction cannot flatter an arm and
    arms are compared on the same pairs. Pass the same ``universe`` for every
    arm of a comparison; the default is truth plus this arm's own pairs, which
    is only appropriate for a single-arm description.
    """
    if universe is None:
        universe = candidate_universe(truth, pred)
    labelled = {p for p, _ in truth}
    true_scores: List[float] = []
    false_scores: List[float] = []
    vec = [0.0] * PAIR_VEC_LEN
    vec[_I_NPOS] = float(len(truth))
    for key in universe:
        if key[0] not in labelled:
            continue
        s = pred.get(key, 0.0)
        if not (0.0 <= s <= 1.0):
            raise ValueError(f"score {s!r} for {key} is outside [0, 1]")
        y = 1.0 if key in truth else 0.0
        (true_scores if y else false_scores).append(s)
        b = min(N_BINS - 1, int(s * N_BINS))
        vec[_I_BIN + 3 * b] += 1
        vec[_I_BIN + 3 * b + 1] += s
        vec[_I_BIN + 3 * b + 2] += y
        vec[_I_BRIER] += (s - y) ** 2
        vec[_I_NPAIRS] += 1
    true_scores.sort()
    false_scores.sort()
    for i, t in enumerate(GRID):
        vec[_I_TP + i] = float(len(true_scores) - bisect_left(true_scores, t))
        vec[_I_FP + i] = float(len(false_scores) - bisect_left(false_scores, t))
    return vec


def sum_vectors(vectors: Iterable[Sequence[float]]) -> List[float]:
    total: Optional[List[float]] = None
    for v in vectors:
        if total is None:
            total = list(v)
        else:
            for i, x in enumerate(v):
                total[i] += x
    return total or []


def _prf(tp: float, fp: float, n_pos: float) -> Tuple[float, float, float]:
    p = tp / (tp + fp) if tp + fp > 0 else 0.0
    r = tp / n_pos if n_pos > 0 else 0.0
    f = 2 * p * r / (p + r) if p + r > 0 else 0.0
    return p, r, f


def precision_recall_at(vec: Sequence[float], tau: float) -> Dict[str, float]:
    """Fixed-threshold micro precision/recall/F1 (the deployable metrics)."""
    i = _grid_index(tau)
    p, r, f = _prf(vec[_I_TP + i], vec[_I_FP + i], vec[_I_NPOS])
    return {"threshold": tau, "precision": p, "recall": r, "f1": f,
            "n_predicted": vec[_I_TP + i] + vec[_I_FP + i]}


def fmax(vec: Sequence[float]) -> Dict[str, float]:
    """Global micro Fmax over the grid. Threshold is chosen on the evaluated
    labels, so this is a descriptive upper bound (``threshold_on_test``)."""
    best = {"fmax": 0.0, "threshold": float("nan"), "precision": 0.0, "recall": 0.0}
    for i, t in enumerate(GRID):
        if vec[_I_TP + i] + vec[_I_FP + i] == 0:
            continue
        p, r, f = _prf(vec[_I_TP + i], vec[_I_FP + i], vec[_I_NPOS])
        if f > best["fmax"]:
            best = {"fmax": f, "threshold": t, "precision": p, "recall": r}
    best["threshold_on_test"] = True
    return best


def calibration(vec: Sequence[float]) -> Dict[str, object]:
    """Expected calibration error (equal-width bins), Brier score, reliability table."""
    n = vec[_I_NPAIRS]
    bins = []
    ece = 0.0
    for b in range(N_BINS):
        nb, sp, sy = vec[_I_BIN + 3 * b: _I_BIN + 3 * b + 3]
        if nb:
            ece += nb / n * abs(sp / nb - sy / nb)
            bins.append({"bin": b, "n": nb, "mean_score": sp / nb, "frac_true": sy / nb})
    return {"ece": ece if n else float("nan"),
            "brier": vec[_I_BRIER] / n if n else float("nan"),
            "n_pairs": n, "reliability": bins}


# --------------------------------------------------------------- EC top-1

def ec_top1_outcomes(pred: Mapping[Pair, float], truth: Set[Pair]) -> Dict[str, Optional[bool]]:
    """Per labelled protein: was the highest-scoring *complete* EC a true EC?

    Value None = the protein received no complete-EC prediction (counts as
    not correct, reported as lacking coverage). Ties are broken by the
    lexicographically smaller EC so the result is deterministic.
    """
    true_by_p: Dict[str, Set[str]] = defaultdict(set)
    for pid, ec in truth:
        true_by_p[pid].add(normalize_ec(ec))
    best: Dict[str, Tuple[float, str]] = {}
    for (pid, ec), s in pred.items():
        if pid not in true_by_p or not is_complete_ec(ec):
            continue
        ec = normalize_ec(ec)
        cur = best.get(pid)
        if cur is None or s > cur[0] or (s == cur[0] and ec < cur[1]):
            best[pid] = (s, ec)
    return {pid: (best[pid][1] in true_by_p[pid] if pid in best else None) for pid in true_by_p}


def top1_vector(outcomes: Mapping[str, Optional[bool]]) -> List[float]:
    """[n_labelled, n_with_prediction, n_correct]"""
    return [float(len(outcomes)),
            float(sum(1 for v in outcomes.values() if v is not None)),
            float(sum(1 for v in outcomes.values() if v))]


def top1_metrics(vec: Sequence[float]) -> Dict[str, float]:
    n, cov, ok = vec
    return {"n_labelled": n, "coverage": cov / n if n else float("nan"),
            "top1_accuracy": ok / n if n else float("nan"),
            "top1_precision_when_called": ok / cov if cov else float("nan")}


# ---------------------------------------------------- correction and harm

def pair_correction_vector(base: Mapping[Pair, float], arm: Mapping[Pair, float],
                           truth: Set[Pair], tau: float) -> List[float]:
    """Pair-level corrections and harms of ``arm`` against a matched baseline.

    Universe = pairs called by either arm at ``tau`` plus all true pairs, on
    labelled proteins. A pair is *correct* for an arm when its call equals
    its label. corrected = base wrong and arm right; harmed = base right and
    arm wrong. Returns [corrected, harmed].
    """
    universe = {p for p, _ in truth}
    keys = {k for k, s in base.items() if s >= tau and k[0] in universe}
    keys |= {k for k, s in arm.items() if s >= tau and k[0] in universe}
    keys |= truth
    corrected = harmed = 0
    for k in keys:
        y = k in truth
        b_ok = (base.get(k, -1.0) >= tau) == y
        a_ok = (arm.get(k, -1.0) >= tau) == y
        if a_ok and not b_ok:
            corrected += 1
        elif b_ok and not a_ok:
            harmed += 1
    return [float(corrected), float(harmed)]


def top1_correction_vector(base: Mapping[str, Optional[bool]],
                           arm: Mapping[str, Optional[bool]]) -> List[float]:
    """[corrected, harmed] on the proteins labelled in both outcomes."""
    corrected = harmed = 0
    for pid in base.keys() & arm.keys():
        b, a = bool(base[pid]), bool(arm[pid])
        corrected += a and not b
        harmed += b and not a
    return [float(corrected), float(harmed)]


def mcnemar_exact(corrected: float, harmed: float) -> float:
    """Two-sided exact binomial test on discordant pairs."""
    n = int(corrected + harmed)
    if n == 0:
        return 1.0
    k = int(min(corrected, harmed))
    tail = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return min(1.0, 2 * tail)


# ------------------------------------------------------------- bootstrap

def genome_bootstrap(vectors: Mapping[str, Sequence[float]],
                     stat: Callable[[List[float]], float],
                     n_boot: int = 1000, seed: int = 13,
                     alpha: float = 0.05) -> Dict[str, float]:
    """Percentile CI from resampling *genomes* with replacement.

    ``stat`` maps the summed vector of a resample to a number. For paired
    comparisons, concatenate arm and baseline vectors per genome and let
    ``stat`` split them, so both arms see identical resamples.
    """
    ids = sorted(vectors)
    point = stat(sum_vectors(vectors[g] for g in ids)) if ids else float("nan")
    if len(ids) < 2 or n_boot <= 0:
        return {"estimate": point, "ci_low": float("nan"), "ci_high": float("nan"),
                "n_genomes": len(ids), "n_boot": 0}
    rng = random.Random(seed)
    draws = []
    for _ in range(n_boot):
        sample = [vectors[ids[rng.randrange(len(ids))]] for _ in ids]
        v = stat(sum_vectors(sample))
        if not math.isnan(v):
            draws.append(v)
    draws.sort()
    if not draws:
        return {"estimate": point, "ci_low": float("nan"), "ci_high": float("nan"),
                "n_genomes": len(ids), "n_boot": 0}
    lo = draws[int(math.floor(alpha / 2 * (len(draws) - 1)))]
    hi = draws[int(math.ceil((1 - alpha / 2) * (len(draws) - 1)))]
    return {"estimate": point, "ci_low": lo, "ci_high": hi,
            "n_genomes": len(ids), "n_boot": len(draws)}


def paired_stat(fn: Callable[[Sequence[float]], float], split: int) -> Callable[[List[float]], float]:
    """Difference statistic fn(arm) - fn(base) for vectors laid out as arm + base."""
    return lambda v: fn(v[:split]) - fn(v[split:])
