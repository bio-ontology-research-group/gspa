"""Unit tests for scripts/meteor_threshold_repair.py on a toy network.

Run with:  python3 -m pytest scripts/tests/test_meteor_threshold_repair.py -q
       or: python3 scripts/tests/test_meteor_threshold_repair.py
"""
import os
import sys

import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import meteor_threshold_repair as M  # noqa: E402


def toy():
    """A -> B -> C -> biomass, plus a side branch B -> D (D consumed by nothing) and an unused X.

    Reactions (columns): EX_A (uptake, reversed sign: +A), r1 A->B, r2 B->C, r3 B->D, bio C->, r4 X->Y
    Metabolites (rows): A, B, C, D, X, Y
    """
    rxns = ["EX_A_e", "r1_c", "r2_c", "r3_c", "biomass_GmNeg", "r4_c"]
    mets = ["A", "B", "C", "D", "X", "Y"]
    S = np.zeros((6, 6))
    S[0, 0] = 1.0            # EX_A, reversed sign: v > 0 is uptake
    S[0, 1], S[1, 1] = -1.0, 1.0
    S[1, 2], S[2, 2] = -1.0, 1.0
    S[1, 3], S[3, 3] = -1.0, 1.0
    S[2, 4] = -1.0
    S[4, 5], S[5, 5] = -1.0, 1.0
    lb = np.array([0.0, 0.0, 0.0, 0.0, 0.0, -10.0])
    ub = np.array([10.0, 10.0, 10.0, 10.0, 10.0, 10.0])
    return rxns, mets, csr_matrix(S), lb, ub


def test_excluded_reactions_dead_end_and_single_direction():
    rxns, mets, S, lb, ub = toy()
    excl = M.find_excluded_reactions(S, lb, ub, rxns, "biomass_GmNeg")
    # D is carried by one reaction only -> r3 excluded; A is carried by EX_A (+1) and r1 (-1): mixed signs, fine.
    # X/Y are carried by the reversible r4 only -> single carrier -> r4 excluded.
    assert 3 in excl and 5 in excl
    assert 1 not in excl and 2 not in excl and 4 not in excl


def test_maxbio_and_grow_support():
    rxns, mets, S, lb, ub = toy()
    oi = rxns.index("biomass_GmNeg")
    avail = np.ones(6, bool)
    core = np.zeros(6, bool)
    core[[0, 1]] = True                      # draft knows EX_A and r1 only
    assert M.maxbio_active(S, lb, ub, oi, core) == 0.0
    sup, v, obj = M.grow_support(S, lb, ub, oi, avail, gamma_min=0.1, core=core)
    assert sup is not None
    added = set(np.where(sup & ~core)[0])
    assert added == {2, oi}                  # r2 and the biomass reaction bridge the gap
    assert abs(obj - 0.2) < 1e-6             # minimal non-core flux: 0.1 through r2 + 0.1 through biomass
    repaired = core | sup
    assert abs(M.maxbio_active(S, lb, ub, oi, repaired) - 10.0) < 1e-6
    # infeasible when the bridge is unavailable
    avail2 = avail.copy()
    avail2[2] = False
    assert M.grow_support(S, lb, ub, oi, avail2, gamma_min=0.1, core=core)[0] is None


def test_threshold_draft_and_prediction_mapping():
    ecs = [np.array([0]), np.array([1, 2]), np.array([], dtype=np.int64)]
    pred = np.array([[0.9, 0.2, 0.1], [0.1, 0.6, 0.0]], dtype=np.float32)
    draft, maxscore = M.threshold_draft(pred, ecs, 0.5)
    assert draft.tolist() == [True, True, False]
    assert np.allclose(maxscore, [0.9, 0.6, 0.0])
    assert M.current_ec("1.1.1.1", {"1.1.1.1": "1.1.1.2", "1.1.1.2": "1.1.1.3"}) == "1.1.1.3"
    assert M.current_ec("1.1.1.1", {"1.1.1.1": "1.1.1.2", "1.1.1.2": "1.1.1.1"}) in {"1.1.1.1", "1.1.1.2"}  # cycle terminates


def test_default_medium():
    rxns = ["EX_cpd00027_e", "EX_cpd99999_e", "EX_cpd00001_e", "rxn1_c"]
    lb = np.array([-1000.0, -1000.0, -1000.0, 0.0])
    ub = np.array([1000.0, 1000.0, 1000.0, 1000.0])
    lb2, ub2, open_ex = M.apply_default_medium(rxns, lb, ub, {"default": ["cpd00027"]})
    assert open_ex == {0, 2}                 # named medium compound and the always-open H2O
    assert (lb2[0], ub2[0]) == (-100.0, 100.0)
    assert (lb2[1], ub2[1]) == (0.0, 1000.0)  # blocked uptake, secretion allowed
    assert (lb2[3], ub2[3]) == (0.0, 1000.0)  # internal reaction untouched


def test_formula_elements():
    assert M.formula_elements("C6H12O6") == {"C": 6, "H": 12, "O": 6}
    assert M.formula_elements("H2O") == {"H": 2, "O": 1}
    assert M.formula_elements("C10H12N5O13P3R") == {"C": 10, "H": 12, "N": 5, "O": 13, "P": 3, "R": 1}
    assert M.formula_elements("(C2H4)n") is None
    assert M.formula_elements(None) == {}
    assert M.formula_elements("C3H5.5O") == {"C": 3, "H": 5.5, "O": 1}


class _FakeU:
    """Minimal UniversalData stand-in for NetworkProfiler tests."""

    def __init__(self):
        self.allrxns = ["r1_c", "r2_c", "SNK_B_c", "EX_C_e", "biomass_GmNeg", "biomass_GmPos"]
        self.allmet = ["A_c", "B_c", "C_e", "P_c", "Q_c", "bio_c"]
        self.rxn_index = {r: i for i, r in enumerate(self.allrxns)}
        self.rxn_meta = {
            "r1_c": {"lb": 0.0, "ub": 1000.0, "stoich": {"A_c": -1.0, "B_c": 1.0}},
            "r2_c": {"lb": -1000.0, "ub": 1000.0, "stoich": {"A_c": -1.0, "P_c": 1.0}},
            "SNK_B_c": {"lb": 0.0, "ub": 0.0, "stoich": {"B_c": -1.0}},
            "EX_C_e": {"lb": 0.0, "ub": 1000.0, "stoich": {"C_e": -1.0}},
            "biomass_GmNeg": {"lb": 0.0, "ub": 1000.0, "stoich": {"P_c": -1.0, "bio_c": 1.0}},
            "biomass_GmPos": {"lb": 0.0, "ub": 1000.0, "stoich": {"Q_c": -1.0, "bio_c": 1.0}},
        }
        self.met_meta = {"A_c": {"formula": "C2"}, "B_c": {"formula": "C2"}, "C_e": {"formula": "C"},
                         "P_c": {"formula": "C3"}, "Q_c": {"formula": "C"}, "bio_c": {"formula": None}}

    def biomass_id(self, gram):
        return "biomass_GmPos" if gram == "positive" else "biomass_GmNeg"


def test_profiler_submodel_deadends_and_mass_balance():
    P = M.NetworkProfiler(_FakeU())
    assert P.boundary_rxns == ["SNK_B_c", "EX_C_e"]           # single-metabolite reactions only, biomass excluded
    assert P.balanced["r1_c"] is True
    assert P.balanced["r2_c"] is False                        # C2 -> C3
    assert P.balanced["biomass_GmNeg"] is False               # product without formula
    keep = np.array([True, False, False, False, False, False])  # r1 only
    rxns = P.submodel(keep, "biomass_GmNeg")
    assert rxns == ["r1_c", "SNK_B_c", "biomass_GmNeg"]       # sink of B pulled in, biomass appended
    lb = {r: P.U.rxn_meta[r]["lb"] for r in rxns}
    ub = {r: P.U.rxn_meta[r]["ub"] for r in rxns}
    # B: produced by r1 (irreversible) and consumed by the (closed) sink -> not a dead end;
    # bio_c: produced by biomass only -> dead end; A: only consumed -> not a dead end
    assert P.deadends(rxns, lb, ub) == 1
    keep2 = np.array([True, False, True, False, True, False])
    rxns2 = P.submodel(keep2, "biomass_GmNeg")
    assert P.deadends(rxns2, {r: P.U.rxn_meta[r]["lb"] for r in rxns2}, {r: P.U.rxn_meta[r]["ub"] for r in rxns2}) == 1


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
