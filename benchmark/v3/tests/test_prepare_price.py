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


VOCAB = ["1", "1.1", "1.1.1", "1.1.1.1", "1.1.1.2", "1.1.1.3", "1.1.1.n11", "2.7.1.1", "3.1.1.1", "3.1.1.2",
         "3.1.1.3", "3.3.3.3", "4.1.1.1", "5.5.5.5", "9.9.9.9"]
EXTRA_WEAK = {"T1": {"1.1.1.3": 5e-7, "1.1.1.n11": 3e-8},        # nonzero in target_full, absent from target_scores
              "T2": {"1.1": 0.01, "9.9.9.9": 2e-9}}              # '1.1' is a hierarchical class, never a candidate


def f32(x):
    import numpy as np
    return float(np.float32(x))


def full_matrix(target_ids, ta, extra=None):
    import numpy as np
    m = np.zeros((len(target_ids), len(VOCAB)), dtype=np.float32)
    for i, t in enumerate(target_ids):
        for ec, s in {**ta[t], **((extra or {}).get(t) or {})}.items():
            m[i, VOCAB.index(ec)] = s
    return m


def make_artifact(with_b=True, collide=False, full=True):
    """``collide`` gives assembly B a protein id that also occurs in assembly A (cohort namespacing).
    ``full=False`` mimics an artifact without target_full and keeps the messy truncated lists
    (EC: prefix duplicates, a hierarchical partial) used to test normalisation."""
    art = {"_meta": {"extraction_version": "test", "created": "x", "slurm_job_id": "0",
                     "manifest": {"sha256": "x"}, "predictors": {}, "ec_vocab": list(VOCAB) if full else ["HUGE"]}}
    for p in P.PREDICTORS:
        art["_meta"]["predictors"][p] = {"label": p, "files": {"GCF_A.1": {"path": f"/raw/{p}/A.pkl", "sha256": "ab", "size": 1,
                                                                            "mtime": "2026-01-01T00:00:00"}}}
        bare = p == "dpz"  # DeepProZyme drops version suffixes
        v = (lambda s: s.split(".")[0]) if bare else (lambda s: s)
        scale = 6e-5 if p == "clean" else 1.0
        ta = {
            "T1": {"1.1.1.1": 0.9 * scale, "1.1.1.2": 0.4 * scale, "9.9.9.9": 1e-9},
            "T2": ({"2.7.1.1": 0.5 * scale} if full else
                   {"EC:2.7.1.1": 0.5 * scale, "2.7.1.1": 0.3 * scale, "1.1.1.-": 0.2 * scale}),
            "T3": {"3.3.3.3": 0.7 * scale},
            "T4": {"3.1.1.1": 0.3, "3.1.1.2": 0.3, "3.1.1.3": 0.1} if p == "clean" else {"3.1.1.1": 0.6},
            "T7": {"4.1.1.1": 0.5 * scale},
        }
        if full:
            ta = {t: {ec: f32(s) for ec, s in d.items()} for t, d in ta.items()}
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
        if full:
            rec["target_full"] = full_matrix(rec["target_ids"], ta, EXTRA_WEAK)
        art[p] = {"GCF_A.1": rec}
        if with_b:
            b = {"n_proteins": 500, "target_ids": ["T5"], "target_scores": {"T5": {"5.5.5.5": f32(0.8)}},
                 "ec_top_scores": {}, "neighbor_ecs": {}, "meteor_active_ecs": [], "meteor_muted_ecs": [],
                 "meteor_source": "m", "raw_target_stats": {"T5": {"matched_id": "P1.1" if collide else "B5.1", "max": 0.8}}}
            if full:
                b["target_full"] = full_matrix(["T5"], b["target_scores"])
            art[p]["GCA_B.1"] = b
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
        self.assertEqual(by[("P1.1", "1.1.1.1")]["raw_score"], f32(0.9 * 6e-5))   # tiny CLEAN score kept as raw
        self.assertEqual(by[("P1.1", "9.9.9.9")]["raw_score"], f32(1e-9))          # weak candidate not dropped
        self.assertTrue(all("calibrated_prob" not in c for c in self.clean()))
        self.assertTrue(all(c["evidence_type"] == "SEQUENCE_DEEPLEARNING" for c in self.clean()))

    def test_all_nonzero_full_candidates_are_retained_beyond_target_scores(self):
        by = {(c["protein_id"], c["function_id"]): c for c in self.clean()}
        self.assertEqual(by[("P1.1", "1.1.1.3")]["raw_score"], f32(5e-7))         # absent from target_scores
        self.assertEqual(by[("P1.1", "1.1.1.3")]["metadata"]["origin"], "target_full")
        self.assertIn(("P1.1", "1.1.1.n11"), by)                                   # preliminary EC is a real output
        self.assertEqual(by[("P2.1", "9.9.9.9")]["raw_score"], f32(2e-9))
        self.assertFalse([k for k in by if k[1] in ("1", "1.1", "1.1.1")])         # hierarchical classes never claimed
        st = self.rep["truncation"]["per_predictor_genome"]["clean"]["GCF_A.1"]
        tc = st["target_candidates"]
        self.assertEqual((tc["source"], tc["truncated"], tc["fallback_reason"]), ("target_full", False, None))
        self.assertEqual(tc["lost_to_artifact_truncation"], 3)                     # T1: 1.1.1.3, 1.1.1.n11; T2: 9.9.9.9
        self.assertEqual(st["dropped_hierarchical_partial_nonzero"], 1)
        self.assertTrue(self.rep["candidate_coverage"]["clean"]["full_candidate_coverage"])
        self.assertEqual(self.rep["candidate_coverage"]["clean"]["candidates_a_truncated_artifact_would_lose"], 3)

    def test_ec_prefix_duplicates_merge_to_the_maximum_once(self):
        rep = self.do_export("legacy", art=make_artifact(full=False))              # no target_full: truncated fallback
        two = [c for c in self.claims(rep["claims"]["clean"]["GCF_A.1"]["path"])
               if c["protein_id"] == "P2.1" and c["function_id"] == "2.7.1.1"]
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
        self.assertEqual(origins, {"target_full", "neighbor_ecs"})
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

    def test_raw_retrieval_manifest_is_concise_and_exact(self):
        r = self.rep["raw_retrieval_manifest"]
        self.assertEqual(r["files"], len(P.PREDICTORS))
        with open(r["path"]) as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh]
        self.assertEqual(rows[0], ["predictor", "assembly", "path", "sha256", "size_bytes", "mtime"])
        self.assertIn(["clean", "GCF_A.1", "/raw/clean/A.pkl", "ab", "1", "2026-01-01T00:00:00"], rows)
        with open(self.rep["report_path"]) as fh:
            text = fh.read()
        self.assertNotIn("HUGE", text)                                             # ec_vocab is never dumped
        self.assertLess(len(text), 200_000)
        self.assertIn("GFF", r["also_needed"])

    def test_provisional_calibrator_and_evidence_types_are_documented(self):
        self.assertIn("not the identity and not learned", self.rep["calibration"])
        self.assertIn("Provisional", self.rep["evidence_type_provenance"])
        self.assertIn("unknown", self.rep["evidence_type_provenance"])
        self.assertNotIn("before any", self.rep["evidence_type_provenance"].lower())
        over = self.do_export("ov", evidence_types={"graphec": "SEQUENCE_DEEPLEARNING"}, variants=("graphec",))
        self.assertEqual({c["evidence_type"] for c in self.claims(over["claims"]["graphec"]["GCF_A.1"]["path"])},
                         {"SEQUENCE_DEEPLEARNING"})
        self.assertEqual(self.rep["evidence_types"]["graphec"], "STRUCTURE_DEEPLEARNING")
        with self.assertRaises(P.ExportError):
            self.do_export("bad", evidence_types={"nope": "X"})


class NoTruthInInference(Base):
    def test_build_claims_has_no_truth_inputs(self):
        params = set(inspect.signature(P.build_claims).parameters)
        self.assertEqual(params, {"artifact", "predictor", "asm", "aliases", "namespaced", "include_ec_top", "idmap",
                                  "use_target_full", "evidence_type"})

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
            rec["target_full"] = rec["target_full"][::-1].copy()                  # rows follow target_ids
            rec["target_scores"] = dict(reversed(list(rec["target_scores"].items())))
        r1 = self.do_export("a", art=art)
        r2 = self.do_export("b", art=shuffled)
        for v in r1["claims"]:
            for g in r1["claims"][v]:
                self.assertEqual(r1["claims"][v][g]["sha256"], r2["claims"][v][g]["sha256"])


class FullMatrixSchema(Base):
    VALS = None

    def matrix(self, **kw):
        import numpy as np
        vals = np.array([[0.9, 0.02, 0.005, 0.0], [0.001, 0.0005, 0.0002, 0.0001]])
        return P.matrix_claims("clean", "GCF_A.1", ["P1", "P9.1"], ["EC:1.1.1.1", "1.1.1.2", "1.1.1.3", "1.1.1.4"], vals,
                               idmap={"P1": "P1.1"}, **kw)

    def test_default_retention_keeps_every_nonzero_score_and_is_not_pruned(self):
        cl, info = self.matrix()
        got = {(c["protein_id"], c["function_id"]): c["raw_score"] for c in cl}
        self.assertEqual(got[("P1.1", "1.1.1.1")], 0.9)               # EC: prefix stripped, id canonicalised
        self.assertEqual(got[("P1.1", "1.1.1.3")], 0.005)             # below the old 0.01 rule, still kept
        self.assertEqual(len([k for k in got if k[0] == "P9.1"]), 4)
        self.assertEqual((info["policy"], info["n_nonzero"], info["n_retained"], info["pruned"]), ("all_nonzero", 7, 7, False))
        self.assertTrue(all(c["metadata"]["origin"] == "full_matrix" for c in cl))

    def test_lossy_retention_is_explicit_and_marks_pruning(self):
        cl, info = self.matrix(retention="threshold_topk", top_k=2)
        got = {(c["protein_id"], c["function_id"]) for c in cl}
        self.assertIn(("P1.1", "1.1.1.2"), got)                       # >= 0.01
        self.assertNotIn(("P1.1", "1.1.1.3"), got)                    # below 0.01 and outside the top 2
        self.assertEqual(len([k for k in got if k[0] == "P9.1"]), 2)  # weak protein keeps its top-k
        self.assertTrue(info["pruned"])
        self.assertEqual((info["min_score"], info["top_k"]), (0.01, 2))

    def test_matrix_validation(self):
        import numpy as np
        for bad in (np.nan, np.inf, -0.1, 1.5):
            with self.assertRaises(P.ExportError):
                P.matrix_claims("clean", "A", ["p"], ["1.1.1.1"], np.array([[bad]]))
        with self.assertRaises(P.ExportError):
            P.matrix_claims("clean", "A", ["p", "q"], ["1.1.1.1"], np.array([[0.5]]))
        with self.assertRaises(P.ExportError):
            P.matrix_claims("clean", "A", ["p"], ["1.1.1.1"], np.array([[0.5]]), retention="some")
        cl, info = P.matrix_claims("clean", "A", ["p"], ["1.1", "1.1.1.n3", "1.1.1.1"], np.array([[0.3, 0.2, 0.1]]))
        self.assertEqual({c["function_id"] for c in cl}, {"1.1.1.n3", "1.1.1.1"})   # class '1.1' dropped, n-EC kept
        self.assertEqual(info["dropped_hierarchical_partial_nonzero"], 1)

    def test_pruned_matrix_export_does_not_claim_full_coverage(self):
        import numpy as np
        path = os.path.join(self.tmp, "c.npz")
        np.savez(path, index=np.array(["P1", "P2.1"]), columns=np.array(["1.1.1.1", "2.7.1.1", "3.1.1.1"]),
                 values=np.array([[0.9, 0.002, 0.001], [0.0, 0.5, 0.0]]))
        a = self.artifact_file(make_artifact(with_b=False))
        m = self.manifest_file(rows=ROWS[:-1], name="m5.tsv")
        full = P.export(a, m, os.path.join(self.tmp, "f1"), variants=("clean",), raw_matrix_map={"clean": {"GCF_A.1": path}})
        self.assertTrue(full["candidate_coverage"]["clean"]["full_candidate_coverage"])
        self.assertEqual(full["parameters"]["retention"], "all_nonzero")
        lossy = P.export(a, m, os.path.join(self.tmp, "f2"), variants=("clean",), raw_matrix_map={"clean": {"GCF_A.1": path}},
                         retention="threshold_topk", top_k=1)
        self.assertFalse(lossy["candidate_coverage"]["clean"]["full_candidate_coverage"])
        with open(lossy["manifests"]["clean"]) as fh:
            self.assertIn("PRUNED", json.load(fh)["description"])
        with open(full["manifests"]["clean"]) as fh:
            self.assertIn("all nonzero", json.load(fh)["description"])

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


class FullTargetValidation(Base):
    def claims_of(self, art):
        aliases, _, idmaps = P.alias_map(art)
        return P.build_claims(art, "clean", "GCF_A.1", aliases["GCF_A.1"], False, False, idmaps["GCF_A.1"])

    def test_missing_target_full_falls_back_with_explicit_provenance(self):
        art = make_artifact()
        del art["clean"]["GCF_A.1"]["target_full"]
        claims, st = self.claims_of(art)
        tc = st["target_candidates"]
        self.assertEqual((tc["source"], tc["truncated"], tc["fallback_reason"]), ("target_scores", True, "target_full absent"))
        self.assertIsNone(tc["nonzero_in_target_full"])
        self.assertEqual({c["metadata"]["origin"] for c in claims if c["protein_id"] == "P1.1"}, {"target_scores"})
        self.assertNotIn(("P1.1", "1.1.1.3"), {(c["protein_id"], c["function_id"]) for c in claims})
        rep = self.do_export("fb", art=art)
        self.assertFalse(rep["candidate_coverage"]["clean"]["full_candidate_coverage"])
        self.assertTrue(rep["candidate_coverage"]["clean"]["truncation_reasons"]["target_full absent"])
        self.assertTrue(rep["candidate_coverage"]["mapred"]["full_candidate_coverage"])
        with open(rep["manifests"]["clean"]) as fh:
            self.assertIn("PRUNED or truncated", json.load(fh)["description"])

    def test_structural_problems_fall_back_not_crash(self):
        import numpy as np
        cases = {"shape": lambda a: a["clean"]["GCF_A.1"].update(target_full=np.zeros((2, len(VOCAB)), dtype=np.float32)),
                 "dup vocab": lambda a: a["_meta"].update(ec_vocab=VOCAB[:-1] + [VOCAB[0]]),
                 "no vocab": lambda a: a["_meta"].pop("ec_vocab"),
                 "1-D": lambda a: a["clean"]["GCF_A.1"].update(target_full=np.zeros(len(VOCAB)))}
        for name, mutate in cases.items():
            art = make_artifact()
            mutate(art)
            _, st = self.claims_of(art)
            self.assertTrue(st["target_candidates"]["truncated"], name)
            self.assertTrue(st["target_candidates"]["fallback_reason"], name)

    def test_inconsistent_values_are_rejected_not_masked(self):
        import numpy as np
        for name, mutate in {
            "altered value": lambda a: a["clean"]["GCF_A.1"]["target_full"].__setitem__((0, VOCAB.index("1.1.1.1")), 0.123),
            "row order": lambda a: a["clean"]["GCF_A.1"].update(target_full=a["clean"]["GCF_A.1"]["target_full"][::-1].copy()),
            "wrong vocabulary": lambda a: a["_meta"].update(ec_vocab=list(reversed(VOCAB))),
            "zeroed candidate": lambda a: a["clean"]["GCF_A.1"]["target_full"].__setitem__((0, VOCAB.index("1.1.1.1")), 0.0),
            "ec outside vocab": lambda a: a["clean"]["GCF_A.1"]["target_scores"]["T1"].update({"8.8.8.8": 0.5}),
        }.items():
            art = make_artifact()
            mutate(art)
            with self.assertRaises(P.ExportError, msg=name):
                self.claims_of(art)

    def test_nan_inf_negative_and_above_one_are_rejected(self):
        import numpy as np
        for bad in (np.nan, np.inf, -0.5, 1.5):
            art = make_artifact()
            art["clean"]["GCF_A.1"]["target_full"][0, VOCAB.index("1.1.1.3")] = bad
            with self.assertRaises(P.ExportError, msg=str(bad)):
                self.claims_of(art)

    def test_float32_rounding_is_tolerated(self):
        import numpy as np
        art = make_artifact()
        ts = art["clean"]["GCF_A.1"]["target_scores"]["T1"]
        ts["1.1.1.1"] = ts["1.1.1.1"] * (1 + 5e-7)                                 # within float32 rounding
        claims, _ = self.claims_of(art)
        self.assertTrue(claims)
        ts["1.1.1.1"] = ts["1.1.1.1"] * (1 + 1e-3)                                 # a real disagreement
        with self.assertRaises(P.ExportError):
            self.claims_of(art)

    def test_partial_vocab_entries_are_not_confused_with_preliminary_ecs(self):
        self.assertTrue(P.is_hierarchical_partial("1.1"))
        self.assertTrue(P.is_hierarchical_partial("1.1.1.-"))
        self.assertFalse(P.is_hierarchical_partial("1.1.1.n11"))
        self.assertFalse(P.is_hierarchical_partial("EC:1.1.1.1"))
        self.assertTrue(P.is_preliminary_ec("3.5.1.n2"))


class PrecisionRoundTrip(Base):
    def test_weak_distinct_scores_survive_serialisation(self):
        import numpy as np
        lo = np.float32(1.0e-9)
        vals = [lo, np.nextafter(lo, np.float32(1)), np.nextafter(np.nextafter(lo, np.float32(1)), np.float32(1))]
        self.assertEqual(len({float(v) for v in vals}), 3)
        art = make_artifact(with_b=False)
        rec = art["clean"]["GCF_A.1"]
        for i, ec in enumerate(("1.1.1.1", "1.1.1.2", "3.1.1.1")):
            rec["target_scores"]["T1"][ec] = float(vals[i])
        rec["target_full"][0, VOCAB.index("1.1.1.1")] = vals[0]
        rec["target_full"][0, VOCAB.index("1.1.1.2")] = vals[1]
        rec["target_full"][0, VOCAB.index("3.1.1.1")] = vals[2]
        rep = self.do_export(art=art, variants=("clean",), manifest=self.manifest_file(rows=ROWS[:-1], name="m6.tsv"))
        got = {c["function_id"]: c["raw_score"] for c in self.claims(rep["claims"]["clean"]["GCF_A.1"]["path"])
               if c["protein_id"] == "P1.1"}
        self.assertEqual([got["1.1.1.1"], got["1.1.1.2"], got["3.1.1.1"]], [float(v) for v in vals])  # bit-exact
        self.assertEqual(len({got["1.1.1.1"], got["1.1.1.2"], got["3.1.1.1"]}), 3)                   # still distinct

    def write_pair(self, raw, post, name):
        claims = os.path.join(self.tmp, f"{name}.claims.jsonl")
        with open(claims, "w") as fh:
            for ec, s in raw.items():
                fh.write(json.dumps({"protein_id": "p", "function_type": "EC", "function_id": ec, "source": "clean",
                                     "raw_score": s}) + "\n")
        out = os.path.join(self.tmp, f"{name}.tsv")
        with open(out, "w") as fh:
            fh.write("protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\n")
            for ec, s in post.items():
                fh.write(f"p\tEC\t{ec}\t\t{s}\n")
        return claims, out

    def test_six_decimal_posteriors_create_false_ties_and_are_detected(self):
        raw = {"1.1.1.1": 6.2e-5, "1.1.1.2": 6.1e-5, "1.1.1.3": 1e-9}
        # the engine's sigmoid(4r-2) maps these to ~0.1192029 apart by ~1e-6/1e-5; %.6f collapses the top two
        claims, out = self.write_pair(raw, {"1.1.1.1": "0.119203", "1.1.1.2": "0.119203", "1.1.1.3": "0.119203"}, "six")
        res = P.rank_roundtrip(claims, out)
        self.assertFalse(res["raw_rank_preserved"])
        self.assertEqual((res["ties_created"], res["top1_posterior_tied"]), (2, 1))

    def test_full_precision_posteriors_pass_and_inversions_fail(self):
        raw = {"1.1.1.1": 6.2e-5, "1.1.1.2": 6.1e-5, "1.1.1.3": 1e-9}
        claims, out = self.write_pair(raw, {"1.1.1.1": "0.11920379", "1.1.1.2": "0.11920355", "1.1.1.3": "0.1192029"}, "full")
        res = P.rank_roundtrip(claims, out)
        self.assertTrue(res["raw_rank_preserved"], res)
        self.assertEqual((res["proteins"], res["adjacent_pairs"]), (1, 2))
        claims, out = self.write_pair(raw, {"1.1.1.1": "0.3", "1.1.1.2": "0.4", "1.1.1.3": "0.1"}, "inv")
        res = P.rank_roundtrip(claims, out)
        self.assertEqual((res["inversions"], res["top1_changed"]), (1, 1))
        self.assertFalse(res["raw_rank_preserved"])

    def test_raw_ties_carry_no_order_and_pooled_sources_are_refused(self):
        claims, out = self.write_pair({"1.1.1.1": 0.5, "1.1.1.2": 0.5}, {"1.1.1.1": "0.7", "1.1.1.2": "0.7"}, "tie")
        self.assertTrue(P.rank_roundtrip(claims, out)["raw_rank_preserved"])
        with open(claims, "a") as fh:
            fh.write(json.dumps({"protein_id": "p", "function_type": "EC", "function_id": "1.1.1.9", "source": "dpz",
                                 "raw_score": 0.1}) + "\n")
        with self.assertRaises(P.ExportError):
            P.rank_roundtrip(claims, out)

    def test_cli_check_roundtrip_exit_codes(self):
        raw = {"1.1.1.1": 6.2e-5, "1.1.1.2": 6.1e-5}
        c1, o1 = self.write_pair(raw, {"1.1.1.1": "0.119203", "1.1.1.2": "0.119203"}, "cli1")
        c2, o2 = self.write_pair(raw, {"1.1.1.1": "0.11920379", "1.1.1.2": "0.11920355"}, "cli2")
        self.assertEqual(P.main(["check-roundtrip", "--claims", c1, "--integrated", o1]), 1)
        self.assertEqual(P.main(["check-roundtrip", "--claims", c2, "--integrated", o2]), 0)


@unittest.skipUnless(os.path.isfile(REAL_ARTIFACT) and os.path.isfile(REAL_MANIFEST), "local Price artifact not present")
class RealArtifact(Base):
    def test_real_full_candidates_and_coverage(self):
        rep = P.export(REAL_ARTIFACT, REAL_MANIFEST, os.path.join(self.tmp, "realfull"), variants=("clean", "topec", "enzbert"))
        cc = rep["candidate_coverage"]
        for p in ("clean", "topec", "enzbert"):
            self.assertTrue(cc[p]["full_candidate_coverage"], p)
        self.assertEqual(cc["clean"]["nonzero_candidates"], 774954 - 0)             # all nonzero, none pruned
        self.assertGreater(cc["clean"]["candidates_a_truncated_artifact_would_lose"], 700000)
        self.assertEqual(cc["topec"]["candidates_a_truncated_artifact_would_lose"], 0)
        self.assertEqual(cc["enzbert"]["candidates_a_truncated_artifact_would_lose"], 0)
        self.assertLess(cc["clean"]["min_score_mass_fraction_of_truncated_rule"], 1.0)
        self.assertEqual(rep["raw_retrieval_manifest"]["files"], 6 * 22)
        self.assertEqual(sum(f["claims"] for f in rep["claims"]["topec"].values()) > 2000, True)

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
