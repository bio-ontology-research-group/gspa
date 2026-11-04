"""Unit tests for benchmark.v3.offgrowth2.score (spec §7)."""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from typing import Any, Dict

from offgrowth2.report import generate_report
from offgrowth2.rules import candidate_index, clear_candidate_index_cache
from offgrowth2.score import (
    PRIMARY_ARMS,
    _normalized_candidate_rows,
    build_genome_context,
    cached_best_by_protein,
    clear_genome_caches,
    fixed_precision_rule,
    genome_block_bootstrap_precision,
    genome_block_bootstrap_recovery,
    paired_mcnemar_comparison,
    permute_candidates_map_within_genome,
    score_genome,
)


def make_toy_genome_payload(genome_id: str = "GCF_00001") -> Dict[str, Any]:
    return {
        "schema": "gspa.offgrowth2.genome/1",
        "genome_id": genome_id,
        "candidate_table_kind": "six",
        "inputs": {},
        "panel": [
            {"panel_id": "PM1:A02", "element": "C", "seed_id": "cpd00027", "exchange_rxn": "EX_cpd00027_e0"},
            {"panel_id": "PM1:A03", "element": "C", "seed_id": "cpd00082", "exchange_rxn": "EX_cpd00082_e0"},
            {"panel_id": "PM1:A04", "element": "C", "seed_id": "cpd00100", "exchange_rxn": "EX_cpd00100_e0"},
        ],
        "transporter_panel": ["PM1:A02", "PM1:A03"],
        "measured_panel": None,
        "minimal_medium": {},
        "default_medium": {},
        "model_genes": ["prot_case1", "prot_case2", "prot_rel", "prot_unrel", "prot_unc", "prot_extra"],
        "draft_reactions": ["rxn00001"],
        "known_labels": {
            "prot_case1": ["1.1.1.10"],
            "prot_case2": ["2.7.1.20"],
            "prot_rel": ["1.3.8.1"],
            "prot_unrel": ["4.1.1.1"],
        },
        "candidates": [
            {"protein_id": "prot_case1", "ec": "1.1.1.10", "score": 0.85},
            {"protein_id": "prot_case2", "ec": "2.7.1.20", "score": 0.80},
            {"protein_id": "prot_rel", "ec": "1.3.8.2", "score": 0.72},
            {"protein_id": "prot_unrel", "ec": "3.1.1.1", "score": 0.65},
            {"protein_id": "prot_unc", "ec": "5.1.1.1", "score": 0.60},
            {"protein_id": "prot_extra", "ec": "6.1.1.1", "score": 0.20},
        ],
        "ec_to_reactions": {
            "1.1.1.10": ["rxn00101"],
            "1.1.1.99": ["rxn00101"],  # reaction-equivalent EC to 1.1.1.10
            "2.7.1.20": ["rxn00102"],
            "1.3.8.1": ["rxn00999"],
            "1.3.8.2": ["rxn00201"],
            "4.1.1.1": ["rxn00998"],
            "3.1.1.1": ["rxn00301"],
            "5.1.1.1": ["rxn00401"],
            "6.1.1.1": ["rxn00501"],
        },
        "universal_stoichiometry": {
            "rxn00101": ["cpd00027", "cpd00555"],
            "rxn00102": ["cpd00555", "cpd00082"],
        },
        "layout": ["prot_case1", "prot_case2", "prot_rel", "prot_unrel", "prot_unc", "prot_extra"],
        "runs": [
            {
                "run_kind": "masked",
                "case": {
                    "case_id": f"{genome_id}:case1",
                    "protein_id": "prot_case1",
                    "truth_ecs": ["1.1.1.10"],
                    "case_source": "price148",
                    "hedged": 0,
                },
                "audit": {},
                "stratum": {"stratum": "off_growth", "wt_biomass": 1.0, "ko_biomass": 1.0, "ko_fraction": 1.0},
                "target_reactions": ["rxn00101"],
                "baseline_fill": {"fill": ["rxn00002"], "replay_fill": ["rxn00002"], "stable_fill": ["rxn00002"], "n_lps": 2, "cpu_seconds": 0.1, "grows_without_fill": False},
                "tests": [
                    {
                        "test_id": "t1_a02",
                        "test_type": "T1",
                        "target_id": "PM1:A02",
                        "panel_id": "PM1:A02",
                        "element": "C",
                        "gapfill_ran": True,
                        "feasible_in_universal": True,
                        "already_satisfied": False,
                        "fill": ["rxn00101", "rxn00102"],
                        "replay_fill": ["rxn00101", "rxn00102"],
                        "stable_fill": ["rxn00101", "rxn00102"],
                        "n_lps": 2,
                        "cpu_seconds": 0.1,
                    },
                    {
                        "test_id": "a0_default",
                        "test_type": "A0",
                        "target_id": "biomass",
                        "panel_id": None,
                        "element": None,
                        "gapfill_ran": True,
                        "feasible_in_universal": True,
                        "already_satisfied": False,
                        "fill": [],
                        "replay_fill": [],
                        "stable_fill": [],
                        "n_lps": 2,
                        "cpu_seconds": 0.1,
                    },
                    {
                        "test_id": "a1old_a04",
                        "test_type": "A1old",
                        "target_id": "PM1:A04",
                        "panel_id": "PM1:A04",
                        "element": "C",
                        "gapfill_ran": True,
                        "feasible_in_universal": True,
                        "already_satisfied": False,
                        "fill": ["rxn00301", "rxn00201", "rxn00401", "rxn00501"],
                        "replay_fill": [],
                        "stable_fill": [],
                        "n_lps": 10,
                        "cpu_seconds": 0.5,
                    },
                ],
                "n_lps": 16,
                "cpu_seconds": 0.8,
                "wall_seconds": 0.8,
                "peak_rss_kb": 1000,
            },
            {
                "run_kind": "masked",
                "case": {
                    "case_id": f"{genome_id}:case2",
                    "protein_id": "prot_case2",
                    "truth_ecs": ["2.7.1.20"],
                    "case_source": "price148",
                    "hedged": 0,
                },
                "audit": {},
                "stratum": {"stratum": "off_growth", "wt_biomass": 1.0, "ko_biomass": 1.0, "ko_fraction": 1.0},
                "target_reactions": ["rxn00102"],
                "baseline_fill": {"fill": ["rxn00002"], "replay_fill": ["rxn00002"], "stable_fill": ["rxn00002"], "n_lps": 2, "cpu_seconds": 0.1, "grows_without_fill": False},
                "tests": [
                    {
                        "test_id": "t2_mod",
                        "test_type": "T2",
                        "target_id": "rxn00102",
                        "panel_id": None,
                        "element": None,
                        "gapfill_ran": True,
                        "feasible_in_universal": True,
                        "already_satisfied": False,
                        "fill": ["rxn00102"],
                        "replay_fill": ["rxn00102"],
                        "stable_fill": ["rxn00102"],
                        "n_lps": 2,
                        "cpu_seconds": 0.1,
                    },
                ],
                "n_lps": 4,
                "cpu_seconds": 0.2,
                "wall_seconds": 0.2,
                "peak_rss_kb": 1000,
            },
            {
                "run_kind": "unmasked",
                "case": None,
                "audit": {},
                "stratum": None,
                "target_reactions": [],
                "baseline_fill": {"fill": ["rxn00002"], "replay_fill": ["rxn00002"], "stable_fill": ["rxn00002"], "n_lps": 2, "cpu_seconds": 0.1, "grows_without_fill": False},
                "tests": [
                    {
                        "test_id": "t1_a02",
                        "test_type": "T1",
                        "target_id": "PM1:A02",
                        "panel_id": "PM1:A02",
                        "element": "C",
                        "gapfill_ran": True,
                        "feasible_in_universal": True,
                        "already_satisfied": False,
                        "fill": ["rxn00101", "rxn00201", "rxn00401"],
                        "replay_fill": ["rxn00101", "rxn00201", "rxn00401"],
                        "stable_fill": ["rxn00101", "rxn00201", "rxn00401"],
                        "n_lps": 2,
                        "cpu_seconds": 0.1,
                    },
                    {
                        "test_id": "t2_dup",
                        "test_type": "T2",
                        "target_id": "rxn00101",
                        "panel_id": None,
                        "element": None,
                        "gapfill_ran": True,
                        "feasible_in_universal": True,
                        "already_satisfied": False,
                        "fill": ["rxn00101", "rxn00501"],
                        "replay_fill": ["rxn00101", "rxn00501"],
                        "stable_fill": ["rxn00101", "rxn00501"],
                        "n_lps": 2,
                        "cpu_seconds": 0.1,
                    },
                    {
                        "test_id": "a1old_big",
                        "test_type": "A1old",
                        "target_id": "PM1:A04",
                        "panel_id": "PM1:A04",
                        "element": "C",
                        "gapfill_ran": True,
                        "feasible_in_universal": True,
                        "already_satisfied": False,
                        "fill": ["rxn00301", "rxn00201", "rxn00401", "rxn00501"],
                        "replay_fill": [],
                        "stable_fill": [],
                        "n_lps": 10,
                        "cpu_seconds": 0.5,
                    },
                ],
                "n_lps": 16,
                "cpu_seconds": 0.8,
                "wall_seconds": 0.8,
                "peak_rss_kb": 1000,
            },
        ],
    }


class TestScore(unittest.TestCase):
    def test_score_genome_recovery_decisions_q_controls_and_stats(self) -> None:
        payload = make_toy_genome_payload("GCF_00001")
        ctx = build_genome_context(payload)
        res = score_genome(payload, ctx=ctx, arms=PRIMARY_ARMS, q_fixed_threshold=0.75)

        # N1 recovers both masked cases; S0 recovers neither
        n1_masked = res.masked_evals_by_arm["N1"]
        s0_masked = res.masked_evals_by_arm["S0"]
        self.assertEqual([e.gene_recovered for e in n1_masked], [True, True])
        self.assertEqual([e.reaction_recovered for e in n1_masked], [True, True])
        self.assertTrue(n1_masked[0].pathway_recovered)
        self.assertEqual([e.gene_recovered for e in s0_masked], [False, False])

        # N2 falls back to P_T and flags measured_panel_fallback when measured_panel is None
        self.assertTrue(res.unmasked_eval_by_arm["N2"].measured_panel_fallback)

        # S1b LP cost is bounded by N1's LP budget
        self.assertLessEqual(
            res.unmasked_eval_by_arm["S1b"].n_lps,
            res.unmasked_eval_by_arm["N1"].n_lps,
        )

        # Distinct decisions counted once per (protein, ec) even though rxn00101 was added by two tests
        n1_decs = res.test_decisions_by_arm["N1"]
        by_prot = {d.protein_id: d for d in n1_decs}
        self.assertEqual(len(n1_decs), 3)
        self.assertEqual(by_prot["prot_case1"].classification, "supported")
        self.assertEqual(by_prot["prot_rel"].classification, "harmful")
        self.assertTrue(by_prot["prot_rel"].related_substrate)
        self.assertEqual(by_prot["prot_unc"].classification, "unconfirmed")

        # Precision excludes unconfirmed: 1 supported / (1 supported + 1 harmful) = 0.50
        fp = fixed_precision_rule(n1_decs, min_precision=0.75, min_decisions=20)
        self.assertEqual(fp["n_supported"], 1)
        self.assertEqual(fp["n_harmful"], 1)
        self.assertEqual(fp["n_related_substrate"], 1)
        self.assertEqual(fp["n_unconfirmed"], 1)
        self.assertAlmostEqual(fp["precision"], 0.50)
        self.assertFalse(fp["passes"])

        # Comparator Q is volume-matched (>= N1 distinct gene calls)
        q_eval = res.unmasked_eval_by_arm["Q"]
        self.assertGreaterEqual(q_eval.n_distinct_gene_calls, res.unmasked_eval_by_arm["N1"].n_distinct_gene_calls)
        self.assertIn("Q_fixed", res.unmasked_eval_by_arm)

        # Orphan step rxn00501 (max candidate score 0.20 < 0.30) is detected
        self.assertEqual([o.rxn_id for o in res.orphan_steps], ["rxn00501"])

        # McNemar and bootstrap run cleanly
        mc = paired_mcnemar_comparison(n1_masked, s0_masked, metric="gene_recovered")
        self.assertEqual(mc["b_better_only"], 2)
        self.assertEqual(mc["c_baseline_only"], 0)
        self.assertAlmostEqual(mc["p_value"], 0.25)

        payload2 = make_toy_genome_payload("GCF_00002")
        res2 = score_genome(payload2)
        boot_rec = genome_block_bootstrap_recovery(
            n1_masked + res2.masked_evals_by_arm["N1"],
            metric="gene_recovered",
            n_resamples=100,
        )
        self.assertAlmostEqual(boot_rec["point"], 1.0)
        boot_prec = genome_block_bootstrap_precision(
            n1_decs + res2.test_decisions_by_arm["N1"],
            all_genome_ids=["GCF_00001", "GCF_00002"],
            n_resamples=100,
        )
        self.assertAlmostEqual(boot_prec["precision"], 0.50)

    def test_genome_seed_caches_reuse_across_calls(self) -> None:
        clear_genome_caches()
        clear_candidate_index_cache()
        payload = make_toy_genome_payload("GCF_00001")
        ctx = build_genome_context(payload)

        perm1 = permute_candidates_map_within_genome(ctx.candidates, 1, extra_proteins=ctx.model_genes, genome_id=ctx.genome_id)
        perm2 = permute_candidates_map_within_genome(ctx.candidates, 1, extra_proteins=list(ctx.model_genes), genome_id=ctx.genome_id)
        self.assertIs(perm1, perm2)

        idx1 = candidate_index(perm1, genome_id=ctx.genome_id, seed=1)
        idx2 = candidate_index(dict(perm1), genome_id=ctx.genome_id, seed=1)
        self.assertIs(idx1, idx2)

        best1 = cached_best_by_protein(perm1, genome_id=ctx.genome_id, seed=1)
        best2 = cached_best_by_protein(dict(perm1), genome_id=ctx.genome_id, seed=1)
        self.assertIs(best1, best2)

        rows1 = _normalized_candidate_rows(perm1, genome_id=ctx.genome_id, seed=1)
        rows2 = _normalized_candidate_rows(dict(perm1), genome_id=ctx.genome_id, seed=1)
        self.assertIs(rows1, rows2)

        clear_genome_caches(ctx.genome_id)

    def test_measured_panel_n2_differs_from_transporter_panel_and_flags_fallback(self) -> None:
        tmp_base = Path("/home/leechuck/tmp")
        tmp_base.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(tmp_base)) as tmpdir:
            root = Path(tmpdir)
            res_dir = root / "results"
            out_dir = root / "out"
            res_dir.mkdir(parents=True, exist_ok=True)

            # Genome 1 (GCA_000006965.1 -> Smeli) has rows in matches.tsv on PM1:A04 (not in P_T = [PM1:A02, PM1:A03])
            # Genome 2 (GCF_000999995.1 -> NoMatchOrg) has no rows in matches.tsv -> falls back to P_T with fallback=1
            for gid in ("GCA_000006965.1", "GCF_000999995.1"):
                payload = make_toy_genome_payload(gid)
                payload.pop("measured_panel", None)
                for run in payload["runs"]:
                    run["tests"].append({
                        "test_id": "t1_a04",
                        "test_type": "T1",
                        "target_id": "bio1",
                        "panel_id": "PM1:A04",
                        "element": "C",
                        "gapfill_ran": True,
                        "feasible_in_universal": True,
                        "already_satisfied": False,
                        "fill": ["rxn00301"],
                        "replay_fill": ["rxn00301"],
                        "stable_fill": ["rxn00301"],
                        "n_lps": 4,
                        "cpu_seconds": 0.2,
                    })
                g_dir = res_dir / gid
                g_dir.mkdir(parents=True, exist_ok=True)
                with gzip.open(g_dir / "genome_results2.json.gz", "wt", encoding="utf-8") as fh:
                    json.dump(payload, fh)

            org_map_tsv = root / "price22_fitness_browser_organisms.tsv"
            org_map_tsv.write_text(
                "assembly\torg_id\tmatched_proteins\tassembly_proteins\tfitness_browser_proteins\tfraction_of_fitness_browser_proteome\n"
                "GCA_000006965.1\tSmeli\t6212\t6173\t6218\t0.999\n"
                "GCF_000999995.1\tNoMatchOrg\t1000\t1000\t1000\t1.0\n",
                encoding="utf-8",
            )
            matches_tsv = root / "experiment_panel_matches.tsv"
            matches_tsv.write_text(
                "org_id\texp_name\texp_group\tcondition\texp_desc\tmedia\tpanel_id\n"
                "Smeli\tset1IT001\tcarbon source\tGlycerol\tGlycerol (C)\tRCH2\tPM1:A04\n",
                encoding="utf-8",
            )

            p1 = json.load(gzip.open(res_dir / "GCA_000006965.1" / "genome_results2.json.gz", "rt"))
            ctx1 = build_genome_context(p1, matches_tsv=matches_tsv, org_map_tsv=org_map_tsv)
            self.assertEqual(ctx1.transporter_panel_ids, ["PM1:A02", "PM1:A03"])
            self.assertEqual(ctx1.measured_panel_ids, ["PM1:A04"])
            self.assertNotEqual(ctx1.measured_panel_ids, ctx1.transporter_panel_ids)

            res1 = score_genome(p1, ctx=ctx1)
            self.assertFalse(res1.unmasked_eval_by_arm["N2"].measured_panel_fallback)
            self.assertTrue(all(not e.measured_panel_fallback for e in res1.masked_evals_by_arm["N2"]))

            p2 = json.load(gzip.open(res_dir / "GCF_000999995.1" / "genome_results2.json.gz", "rt"))
            ctx2 = build_genome_context(p2, matches_tsv=matches_tsv, org_map_tsv=org_map_tsv)
            self.assertIsNone(ctx2.measured_panel_ids)
            res2 = score_genome(p2, ctx=ctx2)
            self.assertTrue(res2.unmasked_eval_by_arm["N2"].measured_panel_fallback)
            self.assertTrue(all(e.measured_panel_fallback for e in res2.masked_evals_by_arm["N2"]))

            generate_report(
                results_dir=res_dir,
                out_dir=out_dir,
                tier="A",
                candidate_table="six",
                matches_tsv=matches_tsv,
                org_map_tsv=org_map_tsv,
                n_resamples=20,
            )
            with (out_dir / "cost.tsv").open("r", encoding="utf-8", newline="") as fh:
                n2_rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["arm"] == "N2"]
            self.assertTrue(any(r["genome_id"] == "GCA_000006965.1" for r in n2_rows))
            self.assertTrue(any(r["genome_id"] == "GCF_000999995.1" for r in n2_rows))
            for r in n2_rows:
                if r["genome_id"] == "GCA_000006965.1":
                    self.assertEqual(int(r["measured_panel_fallback"]), 0)
                else:
                    self.assertEqual(int(r["measured_panel_fallback"]), 1)


if __name__ == "__main__":
    unittest.main()
