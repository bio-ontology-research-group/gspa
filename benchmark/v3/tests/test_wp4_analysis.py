"""Tests for the WP4 sealed recurrence analysis (design §5, §8, §10 row 8).

Synthetic production-shaped outputs: ``<genome>.integrated.tsv``, ``batch_stats.tsv``,
``verified_run.json`` and ``features/<genome>/candidate_features.jsonl``.
"""
import csv
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import analyze_wp4_recurrence as M  # noqa: E402

TARGET = "GCA_000007065.1"           # test split (§2)
TRAIN = "GCA_000007005.1"            # train split
BLOCK = "Methanosarcinaceae"
F1, A1 = "EC:1.1.1.1", "EC:2.2.2.2"  # p1: toward (+ complementary)
F2, A2 = "EC:3.3.3.3", "EC:4.4.4.4"  # p2: away
F3, A3 = "EC:5.5.5.5", "EC:6.6.6.6"  # p3: unchanged
F4 = "EC:7.7.7.7"                    # p4: no alternative
F5 = "EC:8.8.8.8"                    # p5: no candidate
F6, A6 = "EC:9.9.9.9", "EC:1.2.3.4"  # p6: F unavailable in G0
U1, U2 = "EC:1.1.1.10", "EC:1.1.1.11"  # p7: unlabelled multi-candidate

# protein -> {function: (q, complementarity)} for R0 and R2 UNIVERSAL.
R0_Q = {"p1": {F1: (0.2, 0.0), A1: (0.2, 0.0)}, "p2": {F2: (0.5, 0.0), A2: (0.1, 0.0)},
        "p3": {F3: (0.3, 0.0), A3: (0.3, 0.0)}, "p4": {F4: (0.4, 0.0)},
        "p6": {F6: None, A6: (0.3, 0.0)}, "p7": {U1: (0.3, 0.0), U2: (0.3, 0.0)}}
R2_Q = {"p1": {F1: (0.6, 0.3), A1: (0.2, 0.0)}, "p2": {F2: (0.1, 0.0), A2: (0.5, 0.0)},
        "p3": {F3: (0.3, 0.0), A3: (0.3, 0.0)}, "p4": {F4: (0.4, 0.0)},
        "p6": {F6: (0.5, 0.0), A6: (0.3, 0.0)}, "p7": {U1: (0.9, 0.0), U2: (0.1, 0.0)}}
# The S1 shuffle turns p2 toward its label and leaves p7 flat: the null is not zero.
S1_Q = {"p1": {F1: (0.2, 0.0), A1: (0.2, 0.0)}, "p2": {F2: (0.9, 0.0), A2: (0.1, 0.0)},
        "p3": {F3: (0.3, 0.0), A3: (0.3, 0.0)}, "p4": {F4: (0.4, 0.0)},
        "p6": {F6: (0.5, 0.0), A6: (0.3, 0.0)}, "p7": {U1: (0.3, 0.0), U2: (0.3, 0.0)}}
LABELS = [(TARGET, "p1", F1, "ECO:0000269", "experimental"),
          (TARGET, "p2", F2, "ECO:0000269", "experimental"),
          (TARGET, "p3", F3, "ECO:0000269", "experimental"),
          (TARGET, "p4", F4, "ECO:0000269", "experimental"),
          (TARGET, "p5", F5, "ECO:0000269", "experimental"),
          (TARGET, "p6", F6, "ECO:0000269", "experimental"),
          (TARGET, "p7", U1, "ECO:0000305", "curated_inferential"),  # excluded tier
          (TARGET, "p3", "EC:5.5.5.-", "ECO:0000269", "experimental"),  # incomplete EC
          (TARGET, "p1", "GO:0003824", "IDA", "experimental")]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def feature_row(genome, protein, fid, value, model=None, local=True):
    available = value is not None
    q, comp = value if available else (0.0, 0.0)
    features = {
        "cohort.w.available": 1.0 if available else 0.0,
        "cohort.w.complementary_available": 1.0 if available else 0.0,
        "cohort.w.joint_coverage": q, "cohort.w.complementarity": comp,
        "cohort.w.complementary_unannotated": 0.0,
        "cohort.w.support_fraction": q if available else 0.0,
        "cohort.w.local_available": 1.0 if local else 0.0,
        "cohort.w.local_window_truncated": 0.0, "cohort.w.local_window_unknown": 0.0,
        "cohort.w.family_known": 1.0,
    }
    row = {"schema": "gspa.context-feature-row/1", "genome_id": genome, "protein_id": protein,
           "function_type": "EC", "function_id": fid, "role": "CATALYTIC",
           "candidate_key": f"{genome}|{protein}|{fid}|CATALYTIC", "features": features,
           "available_contexts": {"sequence_evidence": True, "cohort": available}}
    if model is not None:
        probability, above = model
        row["model"] = {"id": "m", "sha256": "0" * 64, "kind": "logistic",
                        "score_status": "scored" if probability is not None else "unsupported",
                        "margin": None, "probability": probability, "display_threshold": 0.5,
                        "above_display_threshold": above}
    return row


class Fixture(unittest.TestCase):
    def setUp(self):
        base = Path(os.environ.get("TMPDIR", str(Path.home() / "tmp")))
        base.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=base)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.runs = self.root / "runs"
        self.arms_rows = []

    def write_arm(self, regime, ctx, kind, table, variant="base", models=None, receipt=True,
                  genome=TARGET):
        arm_dir = self.runs / BLOCK / regime / ctx / kind / variant.replace(":", "_")
        inference = arm_dir / "inference"
        inference.mkdir(parents=True, exist_ok=True)
        integrated = inference / f"{genome}.integrated.tsv"
        lines = ["protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\tlikelihood_logodds"]
        features = []
        for protein, functions in sorted(table.items()):
            for fid, value in sorted(functions.items()):
                lines.append(f"{protein}\tEC\t{fid}\t\t{0.1 if value is None else value[0]}\t0")
                model = (models or {}).get((protein, fid))
                features.append(feature_row(genome, protein, fid, value, model))
        lines.append("p1\tGO\tGO:0003824\tF\t0.4\t0")
        integrated.write_text("\n".join(lines) + "\n")
        (inference / "batch_stats.tsv").write_text(
            f"genome_id\tout_path\tn_claims\tn_integrated\n{genome}\t{integrated}\t1\t{len(lines) - 1}\n")
        feature_path = arm_dir / "features" / genome / "candidate_features.jsonl"
        feature_path.parent.mkdir(parents=True, exist_ok=True)
        feature_path.write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in features))
        if receipt:
            (inference / "verified_run.json").write_text(json.dumps(
                {"returncode": 0, "spec": {}, "artifacts": {str(integrated): sha(integrated)}}))
        self.arms_rows.append({"regime": regime, "ctx": ctx, "kind": kind, "variant": variant,
                               "inference_dir": str(inference), "run_dir": str(arm_dir),
                               "holdout": "CLADE" if regime == "R0" else "STRAIN_CLUSTER",
                               "block": BLOCK, "target": TARGET, "model": "-",
                               "manifest": "", "sealed_genomes": genome})
        return arm_dir

    def write_arms(self, name="arms.tsv"):
        path = self.root / "plans" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        fields = list(self.arms_rows[0])
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(self.arms_rows)
        return path

    def write_labels(self, rows=LABELS):
        path = self.root / "labels" / "labels.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("genome_id\tprotein_id\tfunction_id\tevidence\ttier\tswissprot_accession\t"
                        "match_basis\n" + "".join(f"{g}\t{p}\t{f}\t{e}\t{t}\tP0\tsequence\n"
                                                  for g, p, f, e, t in rows))
        return path

    def write_claims(self):
        claims = self.root / "claims"
        claims.mkdir(parents=True, exist_ok=True)
        lines = []
        for protein, functions in R0_Q.items():
            for fid in functions:
                lines.append({"protein_id": protein, "function_type": "EC", "function_id": fid[3:]})
        (claims / f"{TARGET}.claims.jsonl").write_text("".join(json.dumps(r) + "\n" for r in lines))
        return claims

    def build_standard(self):
        r0_models = {("p1", F1): (0.4, False), ("p2", F2): (0.8, True), ("p3", F3): (0.5, False),
                     ("p3", A3): (0.5, False), ("p1", A1): (0.3, False), ("p6", F6): (None, None)}
        r2_models = {("p1", F1): (0.9, True), ("p2", F2): (0.2, False), ("p3", F3): (0.5, False),
                     ("p3", A3): (0.5, False), ("p1", A1): (0.7, True), ("p6", F6): (None, None)}
        self.write_arm("R0", "UNIVERSAL", "none", R0_Q)
        self.write_arm("R2", "UNIVERSAL", "none", R2_Q)
        self.write_arm("R2", "UNIVERSAL", "none", S1_Q, variant="S1:1")
        self.write_arm("R0", "UNIVERSAL", "logistic", R0_Q, models=r0_models)
        self.write_arm("R2", "UNIVERSAL", "logistic", R2_Q, models=r2_models)
        self.write_arm("R0", "SEQONLY", "logistic", R0_Q, models=r0_models)
        self.write_arm("R0-SC", "UNIVERSAL", "none", R0_Q)
        return self.write_arms()

    def seal_and_score(self):
        arms = self.build_standard()
        claims = self.write_claims()
        seal_dir = self.root / "seal"
        M.seal([arms], self.runs, seal_dir, claims_dir=claims)
        metrics = M.score(self.write_labels(), seal_dir, self.root / "score")
        return seal_dir, metrics


class SealAndScore(Fixture):
    def test_seal_writes_label_blind_table_and_hash(self):
        arms = self.build_standard()
        receipt = M.seal([arms], self.runs, self.root / "seal", claims_dir=self.write_claims())
        table = self.root / "seal" / M.TABLE_FILE
        self.assertEqual(receipt["table"]["sha256"], sha(table))
        self.assertTrue(receipt["label_blind"])
        with table.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(list(rows[0]), M.SEAL_COLUMNS)
        self.assertEqual({r["function_type"] for r in rows}, {"EC"})
        unavailable = next(r for r in rows if r["arm_id"].endswith("R0/UNIVERSAL/none/base")
                           and r["protein_id"] == "p6" and r["function_id"] == F6)
        self.assertEqual(unavailable["joint_coverage"], "NA")  # null, never 0 (§5)
        self.assertEqual({r["split"] for r in rows}, {"test"})
        for column in M.SEAL_COLUMNS:
            self.assertNotIn("label", column)
        with self.assertRaises(M.AnalysisError):  # never overwrites a seal
            M.seal([arms], self.runs, self.root / "seal")

    def test_seal_refuses_label_files_and_missing_or_wrong_receipts(self):
        self.write_arm("R0", "UNIVERSAL", "none", R0_Q)
        arms = self.write_arms("labels.tsv")
        with self.assertRaises(M.AnalysisError):
            M.seal([arms], self.runs, self.root / "s1")
        arms = self.write_arms()
        receipt = Path(self.arms_rows[0]["inference_dir"]) / "verified_run.json"
        receipt.write_text(json.dumps({"returncode": 0, "artifacts": {"x/" + TARGET + ".integrated.tsv": "0" * 64}}))
        with self.assertRaises(M.AnalysisError):
            M.seal([arms], self.runs, self.root / "s2")
        receipt.unlink()
        with self.assertRaises(M.AnalysisError):
            M.seal([arms], self.runs, self.root / "s3")
        self.assertEqual(M.seal([arms], self.runs, self.root / "s4",
                                allow_missing_receipts=True)["receipts_required"], False)

    def test_score_refuses_a_changed_table_or_seal(self):
        seal_dir, _ = self.seal_and_score()
        labels = self.root / "labels" / "labels.tsv"
        frozen = sha(seal_dir / M.SEAL_FILE)
        M.score(labels, seal_dir, self.root / "again", expect_seal_sha256=frozen)
        with self.assertRaises(M.AnalysisError):
            M.score(labels, seal_dir, self.root / "x", expect_seal_sha256="0" * 64)
        table = seal_dir / M.TABLE_FILE
        table.write_text(table.read_text().replace("0.6", "0.7", 1))
        with self.assertRaises(M.AnalysisError):
            M.score(labels, seal_dir, self.root / "y")
        self.assertEqual(M.main(["score", "--labels", str(labels), "--seal", str(seal_dir),
                                 "--out", str(self.root / "z")]), 2)

    def test_score_refuses_a_changed_analysis_plan(self):
        seal_dir, _ = self.seal_and_score()
        receipt = json.loads((seal_dir / M.SEAL_FILE).read_text())
        receipt["analysis_plan"]["delta"] = 0.01
        (seal_dir / M.SEAL_FILE).write_text(json.dumps(receipt))
        with self.assertRaises(M.AnalysisError):
            M.score(self.root / "labels" / "labels.tsv", seal_dir, self.root / "x")

    def test_t2_partitions_the_complete_denominator(self):
        _, metrics = self.seal_and_score()
        self.assertEqual(metrics["labels"]["strict_ec_assertions"], 6)
        pooled = metrics["T2"]["UNIVERSAL"]["pooled"]
        self.assertEqual(pooled["denominator"], 6)
        self.assertEqual({c: pooled[c] for c in M.T2_CATEGORIES},
                         {"toward": 1, "away": 1, "unchanged": 1, "f_unavailable": 1,
                          "no_alternative": 1, "no_candidate": 1, "arm_missing": 0})
        for ctx in M.RELATIONS:
            pooled = metrics["T2"][ctx]["pooled"]
            self.assertEqual(sum(pooled[c] for c in M.T2_CATEGORIES), pooled["denominator"])
        self.assertEqual(metrics["T2"]["PATHWAY"]["pooled"]["arm_missing"], 6)
        comp = metrics["T2"]["UNIVERSAL"]["complementary_of_toward"]
        self.assertEqual((comp["denominator"], comp["complementary"]), (1, 1))

    def test_t1_ledger_sums(self):
        _, metrics = self.seal_and_score()
        row = metrics["T1"]["rows"][0]
        self.assertEqual(row["assertions"], 6)
        self.assertEqual(row["no_candidate"] + row["f_candidate_only"] + row["multi_candidate"], 6)
        self.assertEqual((row["no_candidate"], row["f_candidate_only"], row["multi_candidate"]), (1, 1, 4))
        multi = row["multi_UNIVERSAL"]
        self.assertEqual(multi["denominator"], 4)
        self.assertEqual((multi["f_unavailable_g0"], multi["f_available_both"]), (1, 3))

    def test_null_is_computed_and_not_zero(self):
        _, metrics = self.seal_and_score()
        t3 = metrics["T3"]["UNIVERSAL"]
        self.assertEqual(t3["null_status"], "computed")
        self.assertEqual(t3["null_by_seed"][0]["null_toward"], 1)
        self.assertFalse(t3["null_all_zero"])
        self.assertFalse(t3["toward_exceeds_all_seeds"])  # 1 real vs 1 null: not above
        self.assertEqual(metrics["T3"]["PATHWAY"]["null_status"], "null_missing")
        t4 = metrics["T4"]["UNIVERSAL"]
        self.assertEqual(t4["denominator"], 5)
        self.assertEqual(sum(t4[c] for c in ("narrowed", "broadened", "unchanged", "h_undefined",
                                             "arm_missing")), 5)
        self.assertEqual(t4["null_by_seed"][0]["null_narrowed"], 1)  # p2 concentrates under S1

    def test_power_gate_reports_not_testable(self):
        _, metrics = self.seal_and_score()
        gate = metrics["T2"]["UNIVERSAL"]["power_gate"]
        self.assertEqual((gate["assertions_available_both"], gate["minimum"]), (3, 20))
        self.assertFalse(gate["testable"])
        self.assertEqual(metrics["T2"]["UNIVERSAL"]["verdict"], "not_testable")

    def test_t6_learned_paired_gains_harms_ties_and_unconfirmed(self):
        _, metrics = self.seal_and_score()
        t6 = metrics["T6"]
        r2 = next(r for r in t6["per_arm"] if r["regime"] == "R2" and r["ctx"] == "UNIVERSAL")
        self.assertEqual(r2["denominator"], 6)
        self.assertEqual(sum(r2[c] for c in M.RECOVERY_CATEGORIES), 6)
        # p1 recovered; p2, p3 below; p4 has no model score; p5 absent; p6 abstains.
        self.assertEqual({c: r2[c] for c in M.RECOVERY_CATEGORIES},
                         {"recovered": 1, "below_threshold": 2, "abstained": 1, "not_scored": 1,
                          "candidate_absent": 1})
        # p3: F and its competitor tie at 0.5 -> pessimistic rank 2, not top-1.
        self.assertEqual(r2["strict_top1"], 2)  # p1 and p2 (p2's only scored candidate)
        self.assertAlmostEqual(r2["mrr"], (1.0 + 1.0 + 0.5) / 6)
        # p1's competitor above threshold is unconfirmed, never a false positive.
        self.assertEqual(r2["candidates"]["calls_unconfirmed"], 1)
        self.assertEqual(r2["candidates"]["calls_confirmed"], 1)
        pair = next(r for r in t6["paired"] if r["contrast"] == "R2_vs_R0" and r["ctx"] == "UNIVERSAL")
        self.assertEqual((pair["gains"], pair["harms"]), (1, 1))
        self.assertEqual(pair["mcnemar_p"], 1.0)
        text = json.dumps({k: v for k, v in metrics.items() if k != "unconfirmed_policy"}).lower()
        for forbidden in ("false_positive", "precision", "brier", "ece\""):
            self.assertNotIn(forbidden, text)

    def test_t7_controls_expect_zero_change(self):
        _, metrics = self.seal_and_score()
        rows = [r for r in metrics["T7"]["rows"] if r["control"] == "R0-SC"]
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["changed"], rows[0]["defect"]), (0, False))
        self.assertEqual(metrics["T7"]["defects"], 0)

    def test_t7_flags_a_changed_invariance_control(self):
        changed = {k: dict(v) for k, v in R0_Q.items()}
        changed["p1"] = {F1: (0.2 + 1e-9, 0.0), A1: (0.2, 0.0)}
        self.write_arm("R0", "UNIVERSAL", "none", R0_Q)
        self.write_arm("R0-SC", "UNIVERSAL", "none", changed)
        M.seal([self.write_arms()], self.runs, self.root / "seal", claims_dir=self.write_claims())
        metrics = M.score(self.write_labels(), self.root / "seal", self.root / "score")
        self.assertEqual(metrics["T7"]["defects"], 1)
        self.assertIn(M.NEGATIVE["T7_defect"], M.negative_results(metrics))

    def test_report_writes_tables_and_negative_wording(self):
        _, _ = self.seal_and_score()
        sentences = M.report(self.root / "score", self.root / "report")
        tables = (self.root / "report" / "v3-wp4-recurrence-tables.md").read_text()
        for heading in ("## T1 ledger", "## T2 narrowing", "## T3 complementarity",
                        "## T4 label-free", "## T5 add-one", "## T6 learned", "## T7 controls",
                        "## T8 true removal", "## T9 fragmentation", "## Power gate"):
            self.assertIn(heading, tables)
        self.assertIn(M.NEGATIVE["T2_not_testable"].format(ctx="UNIVERSAL"), sentences)
        self.assertIn(M.NEGATIVE["T6_no_gain"].format(kind="logistic", ctx="UNIVERSAL",
                                                      contrast="R2 vs R0"), sentences)
        negative = (self.root / "report" / "v3-wp4-negative-results.md").read_text()
        self.assertIn("not testable", negative)
        # The metrics file is bound by the score receipt.
        metrics = self.root / "score" / M.METRICS_FILE
        metrics.write_text(metrics.read_text().replace('"T1"', '"T1" ', 1))
        with self.assertRaises(M.AnalysisError):
            M.report(self.root / "score", self.root / "report2")


class InMemory(unittest.TestCase):
    """Direct §5/§8 rules on hand-built candidates (no files)."""

    @staticmethod
    def cand(q, comp=0.0, claims=True, support=None, prob=None, above=None):
        c = M.Cand()
        c.q, c.comp, c.in_claims = q, comp, claims
        c.support = q if support is None else support
        if prob is not None:
            c.status, c.prob, c.above = "scored", prob, above
        return c

    def analysis(self, proteins):
        """proteins: {protein: (F, {fid: (q0, qf)})} with labelled F."""
        meta = {"block": BLOCK, "target": TARGET, "split": "test", "tier": "",
                "arm_id": "a", "holdout": "CLADE"}
        r0, r2 = {}, {}
        labels = {}
        for protein, (F, functions) in proteins.items():
            labels[(TARGET, protein)] = {F}
            for fid, (q0, qf) in functions.items():
                r0[(protein, "EC", fid)] = self.cand(q0)
                r2[(protein, "EC", fid)] = self.cand(qf)
        arms = {(TARGET, "R0", "UNIVERSAL", "none", "base"): dict(meta, arm_id="r0"),
                (TARGET, "R2", "UNIVERSAL", "none", "base"): dict(meta, arm_id="r2")}
        cands = {(TARGET, "R0", "UNIVERSAL", "none", "base"): r0,
                 (TARGET, "R2", "UNIVERSAL", "none", "base"): r2}
        an = M.Analysis(arms, cands, {}, set(), labels)
        an.c_peers = {}
        return an

    def test_delta_is_strict_and_fixed(self):
        self.assertEqual(M.DELTA, 0.05)
        # s0(F) = 0.5; sF(F) = 0.54 -> unchanged; 0.56 -> toward; 0.44 -> away; 0.45 -> unchanged.
        cases = {"u1": (0.54, "unchanged"), "t1": (0.56, "toward"), "a1": (0.44, "away"),
                 "u2": (0.45 + 1e-9, "unchanged")}
        an = self.analysis({p: ("EC:1.1.1.1", {"EC:1.1.1.1": (0.5, s), "EC:2.2.2.2": (0.5, 1 - s)})
                            for p, (s, _) in cases.items()})
        for protein, (_, expected) in cases.items():
            state = an.classify(TARGET, protein, "EC:1.1.1.1", "UNIVERSAL")
            self.assertEqual(state["category"], expected, protein)

    def test_zero_total_is_unavailable_not_zero_share(self):
        an = self.analysis({"p": ("EC:1.1.1.1", {"EC:1.1.1.1": (0.0, 0.5), "EC:2.2.2.2": (0.0, 0.5)})})
        self.assertEqual(an.classify(TARGET, "p", "EC:1.1.1.1", "UNIVERSAL")["category"], "f_unavailable")

    def test_power_gate_threshold(self):
        for n, testable in ((19, False), (20, True)):
            proteins = {f"p{i}": ("EC:1.1.1.1", {"EC:1.1.1.1": (0.5, 0.8), "EC:2.2.2.2": (0.5, 0.2)})
                        for i in range(n)}
            t2 = M.t2_narrowing(self.analysis(proteins))["UNIVERSAL"]
            self.assertEqual(t2["power_gate"]["testable"], testable, n)
            self.assertEqual(t2["pooled"]["toward"], n)
            self.assertEqual(t2["verdict"], "narrowing" if testable else "not_testable")
            self.assertLess(t2["sign_test_one_sided_p"], 0.05)

    def test_balanced_toward_and_away_is_no_narrowing(self):
        proteins = {}
        for i in range(12):
            up = i % 2 == 0
            proteins[f"p{i}"] = ("EC:1.1.1.1", {"EC:1.1.1.1": (0.5, 0.8 if up else 0.2),
                                               "EC:2.2.2.2": (0.5, 0.2 if up else 0.8)})
        for i in range(10):
            proteins[f"q{i}"] = ("EC:1.1.1.1", {"EC:1.1.1.1": (0.5, 0.5), "EC:2.2.2.2": (0.5, 0.5)})
        t2 = M.t2_narrowing(self.analysis(proteins))["UNIVERSAL"]
        self.assertTrue(t2["power_gate"]["testable"])
        self.assertEqual((t2["pooled"]["toward"], t2["pooled"]["away"]), (6, 6))
        self.assertEqual(t2["verdict"], "no_narrowing")
        self.assertEqual(sum(t2["secondary_rank"][c] for c in ("rank_improved", "rank_worsened",
                                                               "rank_same", "rank_unavailable")),
                         t2["secondary_rank"]["denominator"])

    def test_pessimistic_ties(self):
        self.assertEqual(M.Analysis.rank_pessimistic({"F": 0.5, "A": 0.5, "B": 0.1}, "F"), 2)
        self.assertEqual(M.Analysis.rank_pessimistic({"F": 0.5, "A": 0.4}, "F"), 1)
        self.assertIsNone(M.Analysis.rank_pessimistic({"F": None, "A": 0.4}, "F"))
        an = self.analysis({"p": ("EC:1.1.1.1", {"EC:1.1.1.1": (0.5, 0.5), "EC:2.2.2.2": (0.5, 0.5)})})
        arm = {("p", "EC", "EC:1.1.1.1"): self.cand(0.5, prob=0.7, above=True),
               ("p", "EC", "EC:2.2.2.2"): self.cand(0.5, prob=0.7, above=True)}
        state = M.recovery_states(an, TARGET, arm)[("p", "EC:1.1.1.1")]
        self.assertEqual((state["state"], state["top1"], state["rr"]), ("recovered", False, 0.5))

    def test_unconfirmed_calls_are_never_false_positives(self):
        an = self.analysis({"p": ("EC:1.1.1.1", {"EC:1.1.1.1": (0.5, 0.5), "EC:2.2.2.2": (0.5, 0.5)})})
        arm = {("p", "EC", "EC:2.2.2.2"): self.cand(0.5, prob=0.9, above=True),
               ("x", "EC", "EC:3.3.3.3"): self.cand(0.5, prob=0.9, above=True),
               ("p", "EC", "EC:1.1.1.1"): self.cand(0.5, prob=0.9, above=True)}
        counts = M.call_counts(an, TARGET, arm)
        self.assertEqual((counts["calls_confirmed"], counts["calls_unconfirmed"]), (1, 2))
        self.assertEqual(set(counts), {"scored", "calls_confirmed", "calls_unconfirmed"})

    def test_t5_add_one_and_t8_removal_partition_their_denominators(self):
        an = self.analysis({"p": ("EC:1.1.1.1", {"EC:1.1.1.1": (0.2, 0.6), "EC:2.2.2.2": (0.2, 0.2)}),
                            "q": ("EC:3.3.3.3", {"EC:3.3.3.3": (0.2, 0.2), "EC:4.4.4.4": (0.2, 0.2)})})
        meta = {"block": BLOCK, "target": TARGET, "split": "test", "holdout": "CLADE"}
        add = (TARGET, "ADD1", "UNIVERSAL", "none", "addone:GCF_1")
        an.arms[add] = dict(meta, arm_id="add", tier="C")
        an.cands[add] = {("p", "EC", "EC:1.1.1.1"): self.cand(0.5), ("p", "EC", "EC:2.2.2.2"): self.cand(0.2),
                         ("q", "EC", "EC:3.3.3.3"): self.cand(0.3), ("q", "EC", "EC:4.4.4.4"): self.cand(0.3)}
        t5 = M.t5_add_one(an)["UNIVERSAL"]["C"]
        self.assertEqual(t5["denominator"], 2)
        # q: F and competitor gain equally -> tie is not "above" (pessimistic).
        self.assertEqual((t5["f_gain_above_competitor"], t5["f_gain_not_above_competitor"]), (1, 1))
        removals = []
        for copy, cls, in_gap in (("c1", "labelled", True), ("c2", "random", False), ("c3", "random", None)):
            key = (copy, "R0", "UNIVERSAL", "none", "removal")
            an.arms[key] = dict(meta, arm_id=f"rm-{copy}", tier="")
            an.cands[key] = {("other", "EC", "EC:1.1.1.1"): self.cand(0.5, prob=0.9, above=True)}
            if in_gap is not None:
                an.gap_available.add((f"rm-{copy}", copy))
                an.gaps[(f"rm-{copy}", copy)] = {"EC:1.1.1.1"} if in_gap else set()
            removals.append({"copy_genome": copy, "target_genome": TARGET, "removed_protein": "p",
                             "class": cls, "function_id": "EC:1.1.1.1"})
        t8 = M.t8_removal(an, removals)
        for row in t8["rows"]:
            self.assertEqual(row["f_in_observed_gaps"] + row["f_not_in_observed_gaps"]
                             + row["gaps_unavailable"], row["denominator"])
        comparison = t8["labelled_vs_random"][0]
        self.assertEqual((comparison["labelled_rate"], comparison["random_base_rate"]), (1.0, 0.0))
        self.assertTrue(comparison["labelled_above_base"])
        self.assertEqual(sum(r["f_reassigned_above_threshold_unconfirmed"] for r in t8["rows"]), 3)

    def test_partition_refuses_an_incomplete_denominator(self):
        with self.assertRaises(M.AnalysisError):
            M.partition({"a": 1, "b": 1}, ("a",), 2)

    def test_statistics(self):
        self.assertAlmostEqual(M.sign_test_one_sided(5, 0), 1 / 32)
        self.assertEqual(M.sign_test_one_sided(0, 0), None)
        self.assertEqual(M.mcnemar_exact(3, 3), 1.0)
        adjusted = M.holm({"a": 0.01, "b": 0.04, "c": None})
        self.assertAlmostEqual(adjusted["a"], 0.02)
        self.assertAlmostEqual(adjusted["b"], 0.04)
        self.assertIsNone(adjusted["c"])
        self.assertAlmostEqual(M.normalized_entropy({"a": 0.5, "b": 0.5}), 1.0)
        self.assertIsNone(M.normalized_entropy({"a": 1.0}))

    def test_labels_strict_policy(self):
        base = Path(os.environ.get("TMPDIR", str(Path.home() / "tmp")))
        with tempfile.TemporaryDirectory(dir=base) as tmp:
            path = Path(tmp) / "labels.tsv"
            path.write_text("genome_id\tprotein_id\tfunction_id\tevidence\ttier\n"
                            + "".join(f"{g}\t{p}\t{f}\t{e}\t{t}\n" for g, p, f, e, t in LABELS)
                            + f"OTHER\tq\t{F1}\tECO:0000269\texperimental\n")
            ec, go = M.load_labels(path, {TARGET})
            self.assertEqual(sum(len(v) for v in ec.values()), 6)
            self.assertNotIn(U1, ec.get((TARGET, "p7"), set()))
            self.assertEqual(go[(TARGET, "p1")], {"GO:0003824"})


if __name__ == "__main__":
    unittest.main()
