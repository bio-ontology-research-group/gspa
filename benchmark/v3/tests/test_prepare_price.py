"""Tests for prepare_price.py on a small synthetic artifact that mimics the real compact schema.

The artifact is synthetic *only as a test fixture* (it is built here, never exported as data).
One test (skipped when the local artifact is absent) runs the exporter on the real file.
"""
import copy
import gzip
import inspect
import json
import os
import pickle
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import prepare_price as P  # noqa: E402
import v3_runner as R  # noqa: E402

REAL_ARTIFACT = "/home/leechuck/Public/software/gspa/benchmark/manifests/reference_data/ibex_extracted/price149_all6_compact.pkl.gz"
REAL_MANIFEST = "/home/leechuck/Public/software/gspa/benchmark/manifests/suite6_price149_manifest.tsv"

MANIFEST_COLS = ["price_entry", "host_assembly", "proteome_protein_id", "ground_truth_ecs", "in_seed_universal_model",
                 "seed_reaction_ids", "metabolic_stratum"]
ROWS = [  # entry, assembly, proteome id, truth
    ("T1", "GCF_A.1", "P1.1", "1.1.1.1;1.1.1.2"),     # multi-EC
    ("T2", "GCF_A.1", "P2.1", "1.1.1.-;2.7.1.1"),     # partial EC dropped, complete kept
    ("T3", "GCF_A.1", "P3.1", "3.3.3.3"),             # predictors disagree on the protein
    ("T4", "GCF_A.1", "P4.1", "3.1.1.1"),             # tie at the top
    ("T6", "GCF_A.1", "NC_000001.1", "1.1.1.1"),      # replicon accession, absent from the artifact
    ("T7", "GCF_A.1", "P7.1", "4.-.-.-"),             # partial only
    ("T5", "GCA_B.1", "B5.1", "5.5.5.5"),             # (collide=True gives it T1's protein id)
]


def make_artifact(with_b=True, collide=False):
    """``collide`` gives assembly B a protein id that also occurs in assembly A (cohort namespacing)."""
    art = {"_meta": {"extraction_version": "test", "created": "x", "slurm_job_id": "0",
                     "manifest": {"sha256": "x"}, "predictors": {}, "ec_vocab": ["HUGE"]}}
    for p in P.PREDICTORS:
        art["_meta"]["predictors"][p] = {"label": p, "files": {"GCF_A.1": {"path": f"/raw/{p}/A.pkl", "sha256": "ab", "size": 1}}}
        bare = p == "dpz"  # DeepProZyme drops version suffixes
        v = (lambda s: s.split(".")[0]) if bare else (lambda s: s)
        scale = 6e-5 if p == "clean" else 1.0
        ta = {
            "T1": {"1.1.1.1": 0.9 * scale, "1.1.1.2": 0.4 * scale, "9.9.9.9": 1e-9},
            "T2": {"EC:2.7.1.1": 0.5 * scale, "2.7.1.1": 0.3 * scale, "1.1.1.-": 0.2 * scale},
            "T3": {"3.3.3.3": 0.7 * scale},
            "T4": {"3.1.1.1": 0.3, "3.1.1.2": 0.3, "3.1.1.3": 0.1} if p == "clean" else {"3.1.1.1": 0.6},
            "T7": {"4.1.1.1": 0.5 * scale},
        }
        stats = {t: {"matched_id": v(f"P{t[1]}.1"), "n_nonzero_vocab": 5, "n_nonzero_dropped_cols": 0, "max": 1.0}
                 for t in ta}
        if p == "enzbert":
            stats["T3"]["matched_id"] = "Q3.1"
        rec = {"n_proteins": 1000, "target_ids": sorted(ta), "target_scores": ta,
               "ec_top_scores": {"1.1.1.1": [(v("P1.1"), 0.9), (v("N9.1"), 0.05)]},
               "neighbor_ecs": {"T1": [(-1, v("N9.1"), [("7.6.2.10", 0.994)]), (2, v("N8.1"), [("1.8.3.7", 0.8956)])],
                                "T2": [(1, v("N9.1"), [("7.6.2.10", 0.994)])]},
               "meteor_active_ecs": [], "meteor_muted_ecs": [], "meteor_source": f"/raw/meteor/{p}.pkl",
               "raw_target_stats": stats}
        art[p] = {"GCF_A.1": rec}
        if with_b:
            art[p]["GCA_B.1"] = {"n_proteins": 500, "target_ids": ["T5"], "target_scores": {"T5": {"5.5.5.5": 0.8}},
                                 "ec_top_scores": {}, "neighbor_ecs": {}, "meteor_active_ecs": [], "meteor_muted_ecs": [],
                                 "meteor_source": "m", "raw_target_stats": {"T5": {"matched_id": "P1.1" if collide else "B5.1", "max": 0.8}}}
    return art


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def artifact_file(self, art, name="art.pkl.gz"):
        path = os.path.join(self.tmp, name)
        with gzip.open(path, "wb") as fh:
            pickle.dump(art, fh)
        return path

    def manifest_file(self, rows=ROWS, name="manifest.tsv", **overrides):
        path = os.path.join(self.tmp, name)
        with open(path, "w") as fh:
            fh.write("\t".join(MANIFEST_COLS) + "\n")
            for entry, asm, pid, truth in rows:
                cells = {"price_entry": entry, "host_assembly": asm, "proteome_protein_id": pid, "ground_truth_ecs": truth,
                         "in_seed_universal_model": "true", "seed_reaction_ids": "rxn1", "metabolic_stratum": "Stratum_A"}
                cells.update({k: v(entry) for k, v in overrides.items()})
                fh.write("\t".join(cells[c] for c in MANIFEST_COLS) + "\n")
        return path

    def do_export(self, out="out", art=None, **kw):
        a = self.artifact_file(art or make_artifact())
        m = kw.pop("manifest", None) or self.manifest_file()
        return P.export(a, m, os.path.join(self.tmp, out), **kw)

    @staticmethod
    def claims(path):
        with open(path) as fh:
            return [json.loads(l) for l in fh]


class Aliases(Base):
    def test_matched_id_is_canonical_and_versions_are_normalised(self):
        aliases, conflicts, idmaps = P.alias_map(make_artifact())
        self.assertEqual(aliases["GCF_A.1"]["T1"], "P1.1")           # dpz said 'P1', others 'P1.1'
        self.assertEqual(idmaps["GCF_A.1"]["P1"], "P1.1")
        self.assertEqual(conflicts["GCF_A.1"], ["T3"])               # enzbert says Q3.1
        self.assertNotIn("T3", aliases["GCF_A.1"])

    def test_claims_use_canonical_ids_for_every_predictor(self):
        rep = self.do_export()
        for v in ("clean", "dpz", "enzbert"):
            ids = {c["protein_id"] for c in self.claims(rep["claims"][v]["GCF_A.1"]["path"])}
            self.assertIn("P1.1", ids)
            self.assertNotIn("P1", ids)
            self.assertIn("N9.1", ids)                              # neighbour ids normalised too
            self.assertNotIn("N9", ids)
            self.assertNotIn("Q3.1", ids)
        self.assertEqual(rep["counts"]["alias_version_normalised"]["dpz"]["targets_rewritten"], 4)

    def test_conflicting_target_is_not_exported_or_labelled_under_a_guess(self):
        rep = self.do_export()
        ids = {c["protein_id"] for c in self.claims(rep["claims"]["clean"]["GCF_A.1"]["path"])}
        self.assertNotIn("P3.1", ids)
        reasons = {e["price_entry"]: e["reason"] for e in rep["exclusions"]}
        self.assertIn("disagree", reasons["T3"])


class Exclusions(Base):
    def test_replicon_row_is_excluded_with_reason(self):
        rep = self.do_export()
        reasons = {e["price_entry"]: e["reason"] for e in rep["exclusions"]}
        self.assertIn("replicon accession", reasons["T6"])
        self.assertEqual(rep["counts"]["excluded_rows"], 3)          # T3, T6, T7

    def test_partial_ec_is_dropped_from_labels_and_recorded(self):
        rep = self.do_export()
        with open(rep["labels"]["GCF_A.1"]["path"]) as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh][1:]
        self.assertIn(["P2.1", "2.7.1.1"], rows)
        self.assertFalse([r for r in rows if "-" in r[1]])
        self.assertEqual({"price_entry": "T2", "dropped": "1.1.1.-"}, rep["partial_ec_drops"][0] | {})
        self.assertIn(["P1.1", "1.1.1.2"], rows)                     # multi-EC keeps both complete ECs
        reasons = {e["price_entry"]: e["reason"] for e in rep["exclusions"]}
        self.assertIn("partial ECs only", reasons["T7"])
        self.assertNotIn("P7.1", {r[0] for r in rows})

    def test_alias_agreement_with_manifest_is_audited(self):
        rep = self.do_export()
        self.assertEqual(rep["counts"]["alias_disagrees_with_manifest"], [])
        m = self.manifest_file(rows=[("T1", "GCF_A.1", "WRONG.1", "1.1.1.1")], name="m2.tsv")
        rep = self.do_export("out2", manifest=m)
        self.assertEqual(rep["counts"]["alias_disagrees_with_manifest"], ["T1"])

    def test_duplicate_manifest_rows_rejected(self):
        m = self.manifest_file(rows=ROWS + [ROWS[0]], name="dup.tsv")
        with self.assertRaises(P.ExportError):
            P.read_manifest(m)


class ClaimsContent(Base):
    def setUp(self):
        super().setUp()
        self.rep = self.do_export()

    def clean(self):
        return self.claims(self.rep["claims"]["clean"]["GCF_A.1"]["path"])

    def test_raw_scores_and_weak_candidates_retained_unmodified(self):
        by = {(c["protein_id"], c["function_id"]): c for c in self.clean()}
        self.assertEqual(by[("P1.1", "1.1.1.1")]["raw_score"], 0.9 * 6e-5)   # tiny CLEAN score kept as raw
        self.assertEqual(by[("P1.1", "9.9.9.9")]["raw_score"], 1e-9)           # weak candidate not dropped
        self.assertTrue(all("calibrated_prob" not in c for c in self.clean()))
        self.assertTrue(all(c["evidence_type"] == "SEQUENCE_DEEPLEARNING" for c in self.clean()))

    def test_ec_prefix_duplicates_merge_to_the_maximum_once(self):
        two = [c for c in self.clean() if c["protein_id"] == "P2.1" and c["function_id"] == "2.7.1.1"]
        self.assertEqual(len(two), 1)
        self.assertEqual(two[0]["raw_score"], 0.5 * 6e-5)

    def test_ties_between_different_ecs_are_both_kept(self):
        t = {c["function_id"]: c["raw_score"] for c in self.clean() if c["protein_id"] == "P4.1"}
        self.assertEqual(t["3.1.1.1"], t["3.1.1.2"])
        self.assertEqual(self.rep["truncation"]["per_predictor_genome"]["clean"]["GCF_A.1"]["n_targets_with_ties_at_top"], 1)
        again = self.do_export("again")
        self.assertEqual(self.rep["claims"]["clean"]["GCF_A.1"]["sha256"], again["claims"]["clean"]["GCF_A.1"]["sha256"])

    def test_non_target_context_is_flagged_and_precedence_is_target_first(self):
        origins = {c["metadata"]["origin"] for c in self.clean()}
        self.assertEqual(origins, {"target_scores", "neighbor_ecs"})
        ctx = [c for c in self.clean() if c["metadata"]["origin"] == "neighbor_ecs"]
        self.assertTrue(ctx and all(c["protein_id"].startswith("N") for c in ctx))

    def test_ec_top_context_is_opt_in(self):
        rep = self.do_export("withtop", include_ec_top=True)
        origins = {c["metadata"]["origin"] for c in self.claims(rep["claims"]["clean"]["GCF_A.1"]["path"])}
        self.assertIn("ec_top_scores", origins)

    def test_truncation_is_disclosed_and_no_coordinates_exported(self):
        t = self.rep["truncation"]
        self.assertFalse(t["proteome_complete"])
        self.assertFalse(t["context_complete"])
        self.assertIn("NOT the whole proteome", t["statement"])
        self.assertIn("No gene coordinates", t["statement"])
        for c in self.clean():
            self.assertFalse({"start", "end", "strand", "contig", "offset"} & (set(c) | set(c["metadata"])))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "out", "layout")))

    def test_exact_raw_sources_listed(self):
        src = self.rep["raw_sources_needed_for_full_context"]
        self.assertEqual(src["clean"]["files"]["GCF_A.1"]["path"], "/raw/clean/A.pkl")
        self.assertIn("meteor_sources", src)


class NoTruthInInference(Base):
    def test_build_claims_has_no_truth_inputs(self):
        params = set(inspect.signature(P.build_claims).parameters)
        self.assertEqual(params, {"artifact", "predictor", "asm", "aliases", "namespaced", "include_ec_top", "idmap"})

    def test_claims_are_byte_identical_when_truth_columns_change(self):
        a = self.artifact_file(make_artifact())
        m1 = self.manifest_file(name="m1.tsv")
        rows2 = [(e, asm, pid, "1.2.3.4") for e, asm, pid, _ in ROWS]  # entirely different truth
        m2 = self.manifest_file(rows=rows2, name="m2.tsv", metabolic_stratum=lambda e: "Stratum_Z",
                                seed_reaction_ids=lambda e: "rxn999;rxn998")
        r1 = P.export(a, m1, os.path.join(self.tmp, "o1"))
        r2 = P.export(a, m2, os.path.join(self.tmp, "o2"))
        for v in r1["claims"]:
            for g in r1["claims"][v]:
                self.assertEqual(r1["claims"][v][g]["sha256"], r2["claims"][v][g]["sha256"], (v, g))
        self.assertNotEqual(r1["labels"]["GCF_A.1"]["sha256"], r2["labels"]["GCF_A.1"]["sha256"])

    def test_no_truth_derived_fields_in_claims_or_manifests(self):
        rep = self.do_export()
        banned = ("seed", "stratum", "truth", "ground", "label", "rxn", "in_seed")
        for v, files in rep["claims"].items():
            for f in files.values():
                with open(f["path"]) as fh:
                    text = fh.read().lower()
                self.assertFalse([b for b in banned if b in text], (v, f["path"]))
        with open(rep["manifests"]["clean"]) as fh:
            man = json.load(fh)
        self.assertEqual(R.scan_claims(man["genomes"][0]["claims"])[1], [])
        arg_blob = json.dumps(man["gspa"]) + json.dumps(man["arms"])
        self.assertNotIn(os.path.dirname(man["genomes"][0]["labels"]["ec"]), arg_blob)


class Namespace(Base):
    def test_auto_namespaces_only_on_collision(self):
        rep = self.do_export(art=make_artifact(collide=True))        # P1.1 occurs in both assemblies
        self.assertTrue(rep["parameters"]["namespaced"])
        self.assertFalse(self.do_export("plain")["parameters"]["namespaced"])
        ids = {c["protein_id"] for c in self.claims(rep["claims"]["clean"]["GCA_B.1"]["path"])}
        self.assertEqual(ids, {"GCA_B.1__P1.1"})
        with open(rep["labels"]["GCA_B.1"]["path"]) as fh:
            self.assertIn("GCA_B.1__P1.1\t5.5.5.5", fh.read())
        rep2 = self.do_export("o2", art=make_artifact(with_b=False),
                              manifest=self.manifest_file(rows=ROWS[:-1], name="m3.tsv"))
        self.assertFalse(rep2["parameters"]["namespaced"])

    def test_namespace_modes(self):
        a, _, _ = P.alias_map(make_artifact(with_b=False))
        self.assertFalse(P.namespace_decision(a, "auto"))
        self.assertTrue(P.namespace_decision(a, "always"))
        with self.assertRaises(P.ExportError):
            P.namespace_decision(a, "sometimes")


class VariantsAndArms(Base):
    def setUp(self):
        super().setUp()
        self.rep = self.do_export(variants=("clean", "dpz", "ensemble_pool4", "ensemble_all6"))

    def test_ensemble_is_exactly_the_union_of_single_predictor_claims(self):
        g = "GCF_A.1"
        singles = []
        for p in P.POOL4:
            single = self.do_export(f"s_{p}", variants=(p,))
            singles += self.claims(single["claims"][p][g]["path"])
        pool = self.claims(self.rep["claims"]["ensemble_pool4"][g]["path"])
        key = lambda c: json.dumps(c, sort_keys=True)  # noqa: E731
        self.assertEqual(sorted(map(key, pool)), sorted(map(key, singles)))
        keys = [(c["protein_id"], c["function_id"], c["source"]) for c in pool]
        self.assertEqual(len(keys), len(set(keys)))                  # no baseline added twice
        self.assertEqual({c["source"] for c in pool}, set(P.POOL4))

    def test_baseline_twice_is_rejected(self):
        c, _ = P.build_claims(make_artifact(), "clean", "GCF_A.1", {"T1": "P1.1"})
        with self.assertRaises(P.ExportError):
            P.merge_variant([c, c])

    def test_manifest_arms_share_one_claims_file_and_load_in_runner(self):
        env = {"GSPA_JAR": "/x.jar", "GSPA_REF": "/ref"}
        for v, path in self.rep["manifests"].items():
            man = R.load_manifest(path, env)
            self.assertEqual(R.manifest_problems(man), [])
            self.assertEqual(R.leakage_problems(man), [])
            self.assertEqual({a["name"] for a in man["arms"]}, {"stage0", "stage1"})
            self.assertEqual(man["baseline_arm"], "stage0")
            self.assertEqual(man["dataset_status"], "real")
            self.assertTrue(all(g["split"] == "test" for g in man["genomes"]))
            self.assertIn("TRUNCATED", man["description"])
            # one claims path per genome, independent of the arm: matched evidence across ablation arms
            self.assertEqual(len({g["claims"] for g in man["genomes"]}), len(man["genomes"]))
        single = R.load_manifest(self.rep["manifests"]["clean"], env)
        self.assertEqual({g["genome_id"] for g in single["genomes"]}, {"GCF_A.1", "GCA_B.1"})

    def test_one_predictor_per_manifest(self):
        for v in ("clean", "dpz"):
            for f in self.rep["claims"][v].values():
                self.assertEqual({c["source"] for c in self.claims(f["path"])}, {v})

    def test_report_hashes_match_files(self):
        for files in self.rep["claims"].values():
            for f in files.values():
                self.assertEqual(f["sha256"], P.sha256_file(f["path"]))
        self.assertEqual(self.rep["inputs"]["manifest"]["rows"], len(ROWS))


class Determinism(Base):
    def test_input_order_does_not_change_claims(self):
        art = make_artifact()
        shuffled = copy.deepcopy(art)
        for p in P.PREDICTORS:
            rec = shuffled[p]["GCF_A.1"]
            rec["target_ids"] = list(reversed(rec["target_ids"]))
            rec["target_scores"] = dict(reversed(list(rec["target_scores"].items())))
        r1 = self.do_export("a", art=art)
        r2 = self.do_export("b", art=shuffled)
        for v in r1["claims"]:
            for g in r1["claims"][v]:
                self.assertEqual(r1["claims"][v][g]["sha256"], r2["claims"][v][g]["sha256"])


class FullMatrixSchema(Base):
    def test_matrix_claims_apply_the_artifact_retention_rule(self):
        import numpy as np
        vals = np.array([[0.9, 0.02, 0.005, 0.0], [0.001, 0.0005, 0.0002, 0.0001]])
        cl = P.matrix_claims("clean", "GCF_A.1", ["P1", "P9.1"], ["EC:1.1.1.1", "1.1.1.2", "1.1.1.3", "1.1.1.4"], vals,
                             top_k=2, idmap={"P1": "P1.1"})
        got = {(c["protein_id"], c["function_id"]): c["raw_score"] for c in cl}
        self.assertEqual(got[("P1.1", "1.1.1.1")], 0.9)               # EC: prefix stripped, id canonicalised
        self.assertIn(("P1.1", "1.1.1.2"), got)                       # >= 0.01
        self.assertNotIn(("P1.1", "1.1.1.3"), got)                    # below 0.01 and outside the top 2
        self.assertEqual(len([k for k in got if k[0] == "P9.1"]), 2)  # weak protein keeps its top-k
        self.assertTrue(all(c["metadata"]["origin"] == "full_matrix" for c in cl))

    def test_npz_and_pickle_loading_rules(self):
        import numpy as np
        path = os.path.join(self.tmp, "m.npz")
        np.savez(path, index=np.array(["a", "b"]), columns=np.array(["EC:1.1.1.1"]), values=np.array([[0.5], [0.6]]))
        idx, cols, vals = P.load_matrix(path)
        self.assertEqual((idx, cols, vals.shape), (["a", "b"], ["EC:1.1.1.1"], (2, 1)))
        with self.assertRaises(P.MatrixLoadError):
            P.load_matrix(os.path.join(self.tmp, "x.pkl"))            # untrusted pickle refused
        with self.assertRaises(P.MatrixLoadError):
            P.load_matrix(os.path.join(self.tmp, "x.csv"))

    def test_full_matrix_export_is_whole_proteome_and_labels_still_align(self):
        import numpy as np
        rmm = {}
        for p in ("clean",):
            path = os.path.join(self.tmp, f"{p}.npz")
            np.savez(path, index=np.array(["P1", "P2.1", "P4.1", "ZZ.1"]),
                     columns=np.array(["EC:1.1.1.1", "2.7.1.1", "3.1.1.1"]),
                     values=np.array([[0.9, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.7], [0.3, 0.3, 0.3]]))
            rmm[p] = {"GCF_A.1": path}
        a = self.artifact_file(make_artifact(with_b=False))
        m = self.manifest_file(rows=ROWS[:-1], name="m4.tsv")
        rep = P.export(a, m, os.path.join(self.tmp, "full"), variants=("clean",), raw_matrix_map=rmm)
        self.assertTrue(rep["truncation"]["proteome_complete"])
        cl = self.claims(rep["claims"]["clean"]["GCF_A.1"]["path"])
        self.assertIn("ZZ.1", {c["protein_id"] for c in cl})          # a non-target protein is present
        self.assertIn("P1.1", {c["protein_id"] for c in cl})
        st = rep["truncation"]["per_predictor_genome"]["clean"]["GCF_A.1"]
        self.assertIn("P7.1", st["targets_absent_from_matrix"])       # honest report of missing targets
        with open(rep["manifests"]["clean"]) as fh:
            self.assertNotIn("TRUNCATED", json.load(fh)["description"])


@unittest.skipUnless(os.path.isfile(REAL_ARTIFACT) and os.path.isfile(REAL_MANIFEST), "local Price artifact not present")
class RealArtifact(Base):
    def test_real_export_counts(self):
        rep = P.export(REAL_ARTIFACT, REAL_MANIFEST, os.path.join(self.tmp, "real"), variants=("clean", "topec"))
        c = rep["counts"]
        self.assertEqual((c["manifest_rows"], c["labelled_targets"], c["excluded_rows"], c["genomes"]), (149, 148, 1, 22))
        self.assertEqual([e["price_entry"] for e in rep["exclusions"]], ["NP_384885"])
        self.assertEqual(c["alias_conflicts"], {})
        self.assertEqual(c["alias_version_normalised"]["dpz"]["targets_rewritten"], 19)
        self.assertEqual(c["alias_disagrees_with_manifest"], [])
        self.assertEqual(rep["inputs"]["artifact"]["extraction_version"], "v2-2026-10-19")
        self.assertEqual(rep["inputs"]["manifest"]["sha256"], P.sha256_file(REAL_MANIFEST))


@unittest.skipUnless(os.environ.get("GSPA_V3_ENGINE_SMOKE") == "1", "set GSPA_V3_ENGINE_SMOKE=1 to run the real gspa CLI")
class EngineAcceptsExport(Base):
    def test_gspa_integrate_reads_exported_claims_for_every_predictor(self):
        repo = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
        rep = self.do_export(variants=tuple(P.PREDICTORS) + ("ensemble_pool4",))
        with open(rep["manifests"]["ensemble_pool4"]) as fh:
            man = json.load(fh)
        man["gspa"].update({"command": [os.path.join(repo, "gradlew"), "-q", "--console=plain", ":gspa-cli:run"],
                            "args_style": "gradle", "cwd": repo, "timeout_sec": 900, "shared_args": []})
        man["gspa"].pop("binary")
        man["arms"] = [{"name": "stage0", "args": ["--enable-priors", "none", "--lite"]}]
        mp = os.path.join(self.tmp, "run.json")
        with open(mp, "w") as fh:
            json.dump(man, fh)
        res = R.run_experiment(man, mp, os.path.join(self.tmp, "run"))
        self.assertEqual(res["runs"]["stage0"]["returncode"], 0, res["runs"]["stage0"].get("failed"))
        self.assertEqual(res["n_not_executed"], 0)
        # all four pooled sources produced integrated rows (evidence_type resolved, none dropped)
        top1 = res["evaluation"]["kinds"]["ec"]["arms"]["stage0"]["ec_top1"]
        self.assertGreater(top1["coverage"], 0.5)


if __name__ == "__main__":
    unittest.main()
