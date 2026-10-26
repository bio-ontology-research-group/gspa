"""Runner tests: manifest validation, leakage controls, orchestration, scoring, optional real-engine smoke.

The orchestration tests use a *stub* executable that writes integrated TSVs in the
production column layout. It exercises the runner's plumbing (batch manifest,
status accounting, hashing, parsing, scoring). It says nothing about the engine;
``RealEngineSmoke`` (opt-in, ``GSPA_V3_ENGINE_SMOKE=1``) runs ``gspa integrate`` itself.
"""
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import v3_runner as R  # noqa: E402

STUB = '''#!/usr/bin/env python3
import sys
args = sys.argv[1:]
bm = args[args.index("--batch-manifest") + 1]
priors_none = "none" in args
rows = [l.rstrip("\\n").split("\\t") for l in open(bm)][1:]
for r in rows:
    with open(r[2], "w") as fh:
        fh.write("protein_id\\ttype\\tfunction_id\\tgo_aspect\\tposterior_prob\\n")
        for l in open(r[1]):
            import json
            c = json.loads(l)
            p = c["raw_score"] if priors_none else min(1.0, c["raw_score"] * (0.5 if c.get("metadata", {}).get("decoy") else 1.0))
            fh.write("\\t".join([c["protein_id"], c["function_type"], c["function_id"], "", "%.6f" % p]) + "\\n")
'''


def claim(pid, ftype, fid, score, **meta):
    return json.dumps({"protein_id": pid, "function_type": ftype, "function_id": fid,
                       "source": "diamond", "raw_score": score, "metadata": meta})


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.stub = self.write("stub_gspa.py", STUB)

    def write(self, name, text):
        p = os.path.join(self.tmp, name)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def manifest(self, command=None):
        # Two test genomes; "decoy" claims are halved by the stub unless priors are off.
        c1 = self.write("g1.claims.jsonl", "\n".join([
            claim("a1", "EC", "1.1.1.1", 0.9), claim("a1", "EC", "1.1.1.2", 0.8, decoy=True),
            claim("a2", "EC", "2.7.1.1", 0.7), claim("a2", "EC", "2.7.1.2", 0.6, decoy=True)]) + "\n")
        c2 = self.write("g2.claims.jsonl", "\n".join([
            claim("b1", "EC", "3.1.1.1", 0.95), claim("b1", "EC", "3.1.1.2", 0.85, decoy=True)]) + "\n")
        l1 = self.write("g1.ec.tsv", "protein_id\tec\na1\t1.1.1.1\na2\t2.7.1.1\n")
        l2 = self.write("g2.ec.tsv", "protein_id\tec\nb1\t3.1.1.1\n")
        return {
            "schema_version": 1, "experiment_id": "unit", "dataset_status": "synthetic",
            "gspa": {"command": command or [sys.executable, self.stub]},
            "settings": {"threshold": 0.5, "bootstrap": 50, "seed": 3},
            "arms": [{"name": "stage0", "args": ["--enable-priors", "none"]}, {"name": "full", "args": []}],
            "baseline_arm": "stage0",
            "genomes": [
                {"genome_id": "g1", "split": "test", "claims": c1, "labels": {"ec": l1}},
                {"genome_id": "g2", "split": "test", "claims": c2, "labels": {"ec": l2}},
            ],
        }


class ManifestValidation(Base):
    def test_valid(self):
        self.assertEqual(R.manifest_problems(self.manifest()), [])

    def test_rejects_bad_fields(self):
        m = self.manifest()
        m["baseline_arm"] = "nope"
        m["genomes"][0]["split"] = "dev"
        m["genomes"][1]["genome_id"] = "g1"
        m["settings"]["threshold"] = 0.505
        text = "; ".join(R.manifest_problems(m))
        for frag in ("baseline_arm", "split must", "unique", "grid"):
            self.assertIn(frag, text)

    def test_shipped_manifests_load(self):
        env = {"GSPA_JAR": "/x/gspa.jar", "GSPA_DATA": "/x/d", "GSPA_REF": "/x/r"}
        mdir = os.path.join(HERE, "..", "manifests")
        names = [n for n in os.listdir(mdir) if n.endswith(".json")]
        self.assertTrue(names)
        for n in names:
            man = R.load_manifest(os.path.join(mdir, n), env)
            self.assertEqual(man["dataset_status"], "real")
            self.assertEqual(sum(g["split"] == "test" for g in man["genomes"]), 13)
            self.assertNotIn("${", json.dumps(man))

    def test_env_expansion(self):
        p = self.write("m.json", json.dumps({"x": "${DATA}/a"}))
        with open(p) as fh:
            self.assertEqual(R._expand(json.load(fh), {"DATA": "/d"}), {"x": "/d/a"})


class Leakage(Base):
    def with_train(self, shared=True):
        m = self.manifest()
        tr = self.write("train.claims.jsonl", claim("a1" if shared else "t1", "EC", "1.1.1.1", 0.9) + "\n")
        m["genomes"].append({"genome_id": "gt", "split": "train", "claims": tr})
        return m

    def test_clean_split_passes(self):
        self.assertEqual(R.leakage_problems(self.with_train(shared=False)), [])

    def test_protein_overlap_detected(self):
        probs = R.leakage_problems(self.with_train(shared=True))
        self.assertTrue(any("also in training" in p for p in probs))

    def test_declared_train_protein_list(self):
        m = self.manifest()
        m["train_proteins"] = self.write("ids.txt", "b1\nzz\n")
        self.assertTrue(any("g2" in p for p in R.leakage_problems(m)))

    def test_species_and_cluster_overlap(self):
        m = self.with_train(shared=False)
        m["genomes"][-1]["species"] = "E. coli"
        m["genomes"][0]["species"] = "E. coli"
        self.assertTrue(any("species overlap" in p for p in R.leakage_problems(m)))

    def test_homology_cluster_spans_splits(self):
        m = self.with_train(shared=False)
        m["homology_clusters"] = self.write("clusters.tsv", "t1\tC1\na1\tC1\n")
        self.assertTrue(any("homology clusters" in p for p in R.leakage_problems(m)))

    def test_label_separation(self):
        m = self.manifest()
        m["genomes"][0]["claims"] = m["genomes"][0]["labels"]["ec"]
        self.assertTrue(any("also a label file" in p for p in R.leakage_problems(m)))
        m = self.manifest()
        m["arms"][1]["args"] = ["--x", m["genomes"][0]["labels"]["ec"]]
        self.assertTrue(any("appears in gspa arguments" in p for p in R.leakage_problems(m)))
        m = self.manifest()
        self.write("g1.claims.jsonl", json.dumps({"protein_id": "a1", "label": 1}) + "\n")
        self.assertTrue(any("label-like" in p for p in R.leakage_problems(m)))

    def test_leakage_aborts_run(self):
        m = self.with_train(shared=True)
        res = R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"))
        self.assertEqual(res["status"], "aborted: leakage")
        self.assertNotIn("runs", res)


class Orchestration(Base):
    def run_it(self, m, **kw):
        mp = self.write("m.json", json.dumps(m))
        return R.run_experiment(m, mp, os.path.join(self.tmp, "out"), **kw)

    def test_end_to_end_with_stub(self):
        res = self.run_it(self.manifest())
        self.assertEqual(res["status"], "evaluated")
        self.assertEqual(res["n_not_executed"], 0)
        self.assertEqual(res["execution_status"]["full"], {"g1": "executed", "g2": "executed"})
        # hashed artefacts
        self.assertEqual(len(res["runs"]["full"]["genomes"]["g1"]["sha256"]), 64)
        self.assertEqual(len(res["provenance"]["inputs"]["g1"]["claims"]["sha256"]), 64)
        self.assertEqual(set(res["provenance"]["code"]), {"v3_runner.py", "v3_metrics.py"})
        # stage0 keeps decoys at their raw scores; full halves them -> top-1 unchanged, pairs corrected
        ec = res["evaluation"]["kinds"]["ec"]
        self.assertAlmostEqual(ec["arms"]["stage0"]["ec_top1"]["top1_accuracy"], 1.0)
        cmp_ = ec["vs_baseline"]["full"]
        # at tau .5 decoys (.8,.6,.85) are FP for stage0; halved (.4,.3,.425) are not called by full
        self.assertEqual(cmp_["pair_correction"]["corrected"], 3.0)
        self.assertEqual(cmp_["pair_correction"]["harmed"], 0.0)
        self.assertAlmostEqual(ec["arms"]["stage0"]["fixed_threshold"]["precision"]["estimate"], 3 / 6)
        self.assertAlmostEqual(ec["arms"]["full"]["fixed_threshold"]["precision"]["estimate"], 1.0)
        self.assertTrue(ec["arms"]["full"]["global_fmax"]["threshold_on_test"])
        self.assertTrue(os.path.isfile(os.path.join(self.tmp, "out", "report.json")))

    def test_noop_control_baseline_against_itself(self):
        m = self.manifest()
        m["arms"] = [{"name": "stage0", "args": ["--enable-priors", "none"]},
                     {"name": "stage0_again", "args": ["--enable-priors", "none"]}]
        cmp_ = self.run_it(m)["evaluation"]["kinds"]["ec"]["vs_baseline"]["stage0_again"]
        self.assertEqual((cmp_["pair_correction"]["corrected"], cmp_["pair_correction"]["harmed"]), (0.0, 0.0))
        self.assertEqual(cmp_["pair_correction"]["mcnemar_p"], 1.0)
        for k in ("delta_f1_at_threshold", "delta_global_fmax", "delta_ece"):
            self.assertEqual((cmp_[k]["estimate"], cmp_[k]["ci_low"], cmp_[k]["ci_high"]), (0.0, 0.0, 0.0))

    def test_failed_engine_is_recorded_not_raised(self):
        bad = self.write("bad.py", "import sys\nsys.exit(3)\n")
        res = self.run_it(self.manifest([sys.executable, bad]))
        self.assertEqual(res["runs"]["full"]["status"], "failed")
        self.assertTrue(res["execution_status"]["full"]["g1"].startswith("failed"))
        self.assertEqual(res["evaluation"]["kinds"]["ec"]["arms"], {})

    def test_missing_input_is_skipped(self):
        m = self.manifest()
        os.remove(m["genomes"][1]["claims"])
        res = self.run_it(m)
        self.assertEqual(res["execution_status"]["full"]["g2"], "skipped: claims file missing")
        self.assertEqual(res["execution_status"]["full"]["g1"], "executed")
        self.assertEqual(res["n_not_executed"], 2)  # g2 in each of the two arms
        self.assertEqual(res["evaluation"]["kinds"]["ec"]["vs_baseline"]["full"]["matched_genomes"], ["g1"])

    def test_dry_run_writes_commands_only(self):
        res = self.run_it(self.manifest(), dry_run=True)
        self.assertEqual(res["status"], "dry_run")
        self.assertIn("--batch-manifest", res["runs"]["full"]["command"])
        self.assertNotIn("evaluation", res)

    def test_gradle_args_style_quotes(self):
        m = self.manifest(["./gradlew", ":gspa-cli:run"])
        m["gspa"]["args_style"] = "gradle"
        cmd = R.build_gspa_command(m, ["integrate", "--taxon", "a b"])
        self.assertEqual(cmd[-1], "--args=integrate --taxon 'a b'")

    def test_batch_manifest_layout(self):
        m = self.manifest()
        path, outs = R.write_batch_manifest(m["arms"][0], m["genomes"], os.path.join(self.tmp, "bm"))
        with open(path) as fh:
            rows = [l.split("\t") for l in fh.read().splitlines()]
        self.assertEqual(rows[0], R.BATCH_HEADER)
        self.assertEqual(len(rows[1]), 11)
        self.assertEqual(rows[1][2], outs["g1"])

    def test_parse_integrated_and_labels(self):
        p = self.write("i.tsv", "protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\n"
                       "p\tGO\tGO:0000001\tBP\t0.7\np\tEC\tEC:1.1.1.1\t\t0.6\np\tEC\t1.1.1.1\t\t0.8\n")
        self.assertEqual(R.parse_integrated(p, "ec"), {("p", "1.1.1.1"): 0.8})
        self.assertEqual(R.parse_integrated(p, "go"), {("p", "GO:0000001"): 0.7})
        g = self.write("t.tsv", "accession\taspect\tgo_term\nP1\tMF\tGO:0000002\n")
        self.assertEqual(R.parse_labels(g, "go"), {("P1", "GO:0000002")})


@unittest.skipUnless(os.environ.get("GSPA_V3_ENGINE_SMOKE") == "1", "set GSPA_V3_ENGINE_SMOKE=1 to run the real gspa CLI")
class RealEngineSmoke(Base):
    """Invokes the production ``gspa integrate`` through Gradle on a 2-protein fixture."""

    def test_real_integrate(self):
        repo = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
        m = self.manifest([os.path.join(repo, "gradlew"), "-q", "--console=plain", ":gspa-cli:run"])
        m["gspa"].update({"args_style": "gradle", "cwd": repo, "timeout_sec": 900})
        m["arms"] = [{"name": "stage0", "args": ["--enable-priors", "none", "--lite"]}]
        res = R.run_experiment(m, self.write("m.json", json.dumps(m)), os.path.join(self.tmp, "out"))
        self.assertEqual(res["runs"]["stage0"]["returncode"], 0, res["runs"]["stage0"].get("failed"))
        self.assertEqual(res["execution_status"]["stage0"], {"g1": "executed", "g2": "executed"})
        ec = res["evaluation"]["kinds"]["ec"]["arms"]["stage0"]["ec_top1"]
        self.assertGreater(ec["n_labelled"], 0)


if __name__ == "__main__":
    unittest.main()
