"""Regression tests for the follow-up review: leakage coverage, strict inputs, stale outputs,
path safety, engine identity and comparable calibration universes."""
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import v3_runner as R  # noqa: E402
from test_runner import Base, claim  # noqa: E402


class ClaimsScan(Base):
    def test_label_after_row_2000_is_found(self):
        lines = [claim(f"p{i}", "EC", "1.1.1.1", 0.5) for i in range(2500)]
        lines.append(json.dumps({"protein_id": "z", "function_type": "EC", "function_id": "1.1.1.1",
                                 "source": "diamond", "raw_score": 0.5, "is_true": True}))
        _, probs = R.scan_claims(self.write("big.jsonl", "\n".join(lines) + "\n"))
        self.assertTrue(any("label-like field is_true" in p and "line 2501" in p for p in probs), probs)

    def test_label_in_metadata_is_found(self):
        _, probs = R.scan_claims(self.write("m.jsonl", claim("a", "EC", "1.1.1.1", 0.5, truth=True) + "\n"))
        self.assertTrue(any("label-like" in p for p in probs))

    def test_malformed_and_nonfinite(self):
        text = "\n".join([
            "{not json",
            "[1, 2]",
            json.dumps({"protein_id": "a", "function_type": "EC", "function_id": "1.1.1.1", "source": "diamond", "raw_score": float("nan")}),
            json.dumps({"protein_id": "b", "function_type": "EC", "function_id": "1.1.1.1", "source": "diamond", "raw_score": "0.5"}),
            json.dumps({"protein_id": "c", "function_type": "EC", "function_id": "1.1.1.1", "source": "diamond", "raw_score": 0.5, "calibrated_prob": 1.2}),
            json.dumps({"protein_id": "d", "function_type": "EC", "function_id": "1.1.1.1", "raw_score": 0.5}),
        ]) + "\n"
        _, probs = R.scan_claims(self.write("bad.jsonl", text))
        blob = " ".join(probs)
        for frag in ("malformed JSON", "is not a JSON object","non-finite or non-numeric raw_score",
                     "calibrated_prob outside [0, 1]", "missing required field(s) source"):
            self.assertIn(frag, blob)

    def test_duplicate_source_claims(self):
        c = claim("a", "EC", "1.1.1.1", 0.5)
        _, probs = R.scan_claims(self.write("d.jsonl", c + "\n" + c + "\n"))
        self.assertTrue(any("duplicate" in p for p in probs))

    def test_bad_claims_abort_run(self):
        m = self.manifest()
        self.write("g1.claims.jsonl", "{broken\n")
        res = R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"))
        self.assertEqual(res["status"], "aborted: leakage")
        self.assertTrue(any("malformed JSON" in p for p in res["leakage_problems"]))


class SplitControls(Base):
    def with_groups(self, train="t1", calib="c1"):
        m = self.manifest()
        m["genomes"].append({"genome_id": "gt", "split": "train", "claims": self.write("t.jsonl", claim(train, "EC", "1.1.1.1", 0.9) + "\n")})
        m["genomes"].append({"genome_id": "gc", "split": "calibration", "claims": self.write("c.jsonl", claim(calib, "EC", "1.1.1.1", 0.9) + "\n")})
        return m

    def test_calibration_must_be_disjoint_from_test_and_train(self):
        self.assertEqual(R.leakage_problems(self.with_groups()), [])
        self.assertTrue(any("calibration and test" in p for p in R.leakage_problems(self.with_groups(calib="a1"))))
        self.assertTrue(any("train and calibration" in p for p in R.leakage_problems(self.with_groups(calib="t1"))))

    def test_homology_map_must_cover_all_proteins(self):
        m = self.with_groups()
        m["homology_clusters"] = self.write("h.tsv", "a1\tC1\nt1\tC2\n")  # a2, b1, c1 unmapped
        rep = R.leakage_report(m)
        self.assertTrue(any("unmapped" in p for p in rep["problems"]), rep)
        self.assertTrue(rep["controls"]["homology_clusters"].startswith("not_assessed"))

    def test_complete_homology_map_verifies_and_detects_spanning(self):
        m = self.with_groups()
        ids = ["a1", "a2", "b1", "t1", "c1"]
        m["homology_clusters"] = self.write("h.tsv", "".join(f"{p}\tC{i}\n" for i, p in enumerate(ids)))
        rep = R.leakage_report(m)
        self.assertEqual(rep["problems"], [])
        self.assertEqual(rep["controls"]["homology_clusters"], "verified")
        m["homology_clusters"] = self.write("h2.tsv", "a1\tX\na2\tY\nb1\tZ\nt1\tX\nc1\tW\n")
        self.assertTrue(any("homology clusters span train and test" in p for p in R.leakage_problems(m)))

    def test_calibration_vs_test_cluster_span(self):
        m = self.with_groups()
        m["homology_clusters"] = self.write("h.tsv", "a1\tX\na2\tY\nb1\tZ\nt1\tT\nc1\tX\n")
        self.assertTrue(any("span calibration and test" in p for p in R.leakage_problems(m)))

    def test_controls_never_claim_more_than_checked(self):
        rep = R.leakage_report(self.manifest())  # test genomes only
        c = rep["controls"]
        self.assertTrue(c["protein_id_disjoint"].startswith("not_assessed"))
        self.assertEqual(c["homology_clusters"], "not_declared")
        self.assertEqual(c["species_disjoint"], "not_declared")
        self.assertEqual(c["pretrained_model_leakage"], "not_assessed")
        m = self.manifest()
        m["declared_controls"] = {"diamond": "reference with panel proteins removed"}
        self.assertTrue(R.leakage_report(m)["controls"]["declared:diamond"].startswith("declared_unverified"))

    def test_controls_recorded_in_report(self):
        m = self.manifest()
        res = R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"))
        self.assertIn("pretrained_model_leakage", res["leakage_controls"])
        self.assertTrue(any("closed world" in a.lower() for a in res["evaluation"]["assumptions"]))


class PathAndIdSafety(Base):
    def test_unsafe_ids_rejected(self):
        for bad in ("../x", "a/b", "a\tb", "", ".hidden", "x" * 80):
            m = self.manifest()
            m["genomes"][0]["genome_id"] = bad
            self.assertTrue(R.manifest_problems(m), bad)
            m = self.manifest()
            m["arms"][1]["name"] = bad
            self.assertTrue(R.manifest_problems(m), bad)

    def test_tab_in_path_rejected(self):
        m = self.manifest()
        m["genomes"][0]["operons"] = "/x/a\tb"
        self.assertTrue(any("tabs/newlines" in p for p in R.manifest_problems(m)))

    def test_run_experiment_refuses_unsafe_manifest(self):
        m = self.manifest()
        m["arms"][0]["name"] = "../../escape"
        with self.assertRaises(R.ManifestError):
            R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "escape")))

    def test_shared_claims_file_is_ambiguous(self):
        m = self.manifest()
        m["genomes"][1]["claims"] = m["genomes"][0]["claims"]
        self.assertTrue(any("ambiguous" in p for p in R.manifest_problems(m)))

    def test_declared_aux_file_missing_is_a_skip_not_a_silent_run(self):
        m = self.manifest()
        m["genomes"][0]["operons"] = os.path.join(self.tmp, "absent_operons.tsv")
        res = R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"))
        self.assertEqual(res["execution_status"]["full"]["g1"], "skipped: declared input missing: operons")
        self.assertEqual(res["execution_status"]["full"]["g2"], "executed")


class StaleOutputs(Base):
    def run_it(self, m, **kw):
        return R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"), **kw)

    def out(self, arm, gid="g1"):
        return os.path.join(self.tmp, "out", "arms", arm, f"{gid}.integrated.tsv")

    def test_failed_rerun_does_not_inherit_previous_success(self):
        self.assertEqual(self.run_it(self.manifest())["n_not_executed"], 0)
        self.assertTrue(os.path.isfile(self.out("full")))
        bad = self.write("bad.py", "import sys\nsys.exit(3)\n")
        res = self.run_it(self.manifest([sys.executable, bad]))
        self.assertFalse(os.path.exists(self.out("full")))
        self.assertTrue(res["execution_status"]["full"]["g1"].startswith("failed"))
        self.assertEqual(res["evaluation"]["kinds"]["ec"]["arms"], {})

    def test_output_written_before_nonzero_exit_is_not_scored(self):
        m = self.manifest()
        m["arms"][1]["args"] = ["--fail-after-write"]
        res = self.run_it(m)
        self.assertEqual(res["runs"]["full"]["returncode"], 2)
        self.assertFalse(os.path.exists(self.out("full")))
        self.assertIn("full", res["execution_status"])
        self.assertTrue(all(s.startswith("failed") for s in res["execution_status"]["full"].values()))
        self.assertNotIn("full", res["evaluation"]["kinds"]["ec"]["arms"])

    def test_evaluate_only_reuses_verified_outputs(self):
        m = self.manifest()
        self.run_it(m)
        res = self.run_it(m, phases=("evaluate",))
        self.assertEqual(res["n_not_executed"], 0)
        self.assertIn("full", res["evaluation"]["kinds"]["ec"]["arms"])

    def test_evaluate_only_rejects_modified_output_changed_claims_and_unrun_arms(self):
        m = self.manifest()
        self.run_it(m)
        with open(self.out("full"), "a") as fh:
            fh.write("a1\tEC\t9.9.9.9\t\t0.99\n")
        self.write("g2.claims.jsonl", claim("b1", "EC", "3.1.1.1", 0.1) + "\n")
        res = self.run_it(m, phases=("evaluate",))
        self.assertIn("modified", res["execution_status"]["full"]["g1"])
        self.assertIn("claims changed", res["execution_status"]["full"]["g2"])
        fresh = os.path.join(self.tmp, "fresh")
        res = R.run_experiment(m, self.write("m2.json", json.dumps(m)), fresh, phases=("evaluate",))
        self.assertTrue(all("no run record" in s for s in res["execution_status"]["full"].values()))

    def test_evaluate_only_rejects_recorded_failure(self):
        m = self.manifest()
        m["arms"][1]["args"] = ["--fail-after-write"]
        self.run_it(m)
        res = self.run_it(m, phases=("evaluate",))
        self.assertTrue(all(s.startswith("failed") for s in res["execution_status"]["full"].values()))


class PredictionParsing(Base):
    def run_bad(self, flag):
        m = self.manifest()
        m["arms"][1]["args"] = [flag]
        return R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"))

    def test_nan_and_out_of_range_scores_fail_the_cell(self):
        for flag in ("--bad-score", "--score-above-one"):
            res = self.run_bad(flag)
            self.assertTrue(all("unparseable output" in s for s in res["execution_status"]["full"].values()), flag)
            self.assertNotIn("full", res["evaluation"]["kinds"]["ec"]["arms"])
            self.assertIn("stage0", res["evaluation"]["kinds"]["ec"]["arms"])
            self.assertTrue(res["evaluation"]["eval_failures"])

    def test_truncated_and_headerless_outputs(self):
        with self.assertRaises(R.PredictionFormatError):
            R.parse_integrated(self.write("a.tsv", "protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\np\tEC\n"), "ec")
        with self.assertRaises(R.PredictionFormatError):
            R.parse_integrated(self.write("b.tsv", "foo\tbar\n"), "ec")
        with self.assertRaises(R.PredictionFormatError):
            R.parse_integrated(self.write("c.tsv", "protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\np\tEC\t1.1.1.1\t\tabc\n"), "ec")


class ComparableCalibration(Base):
    def test_arms_share_one_universe(self):
        m = self.manifest()
        m["arms"].append({"name": "dropped", "args": ["--drop-decoys"]})
        res = R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"))
        arms = res["evaluation"]["kinds"]["ec"]["arms"]
        n = {a: arms[a]["calibration"]["n_pairs"] for a in arms}
        self.assertEqual(len(set(n.values())), 1, n)
        self.assertEqual(n["stage0"], 6)  # 3 true pairs + 3 decoys; dropped arm scores its omitted decoys 0
        u = res["evaluation"]["kinds"]["ec"]["universe"]
        self.assertEqual(sum(v["pairs"] for v in u.values()), 6)
        # omitting decoys legitimately removes false positives, but the universe is unchanged
        self.assertAlmostEqual(arms["dropped"]["fixed_threshold"]["precision"]["estimate"], 1.0)


class EngineIdentity(Base):
    def test_jar_is_hashed_and_binary_mismatch_rejected(self):
        jar = self.write("engine.jar", "jar-bytes")
        other = self.write("other.jar", "x")
        m = self.manifest(["java", "-jar", jar])
        ident = R.engine_identity(m, self.tmp)
        self.assertEqual(ident["jar"]["sha256"], R.sha256_file(jar))
        self.assertEqual(ident["command"], ["java", "-jar", jar])
        self.assertNotIn("warning", ident)
        m["gspa"]["binary"] = other
        self.assertTrue(any("differs from the -jar" in p for p in R.manifest_problems(m)))

    def test_unidentified_engine_is_flagged_and_arm_launch_args_recorded(self):
        m = self.manifest()
        res = R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"))
        self.assertIn("warning", res["engine"])  # a python stub identifies no engine content
        cmd = res["runs"]["full"]["command"]
        self.assertIn("--batch-manifest", cmd)
        self.assertTrue(res["engine"]["launcher"]["sha256"])
        with open(os.path.join(self.tmp, "out", "arms", "full", R.RUN_RECORD)) as fh:
            self.assertEqual(json.load(fh)["command"], cmd)


if __name__ == "__main__":
    unittest.main()
