#!/usr/bin/env python3
"""Unit tests for `offgrowth2.runner` (end-to-end CLI, schema, masking audit skip, resume, determinism)."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
from typing import Dict, List

import cobra

V3_DIR = Path(__file__).resolve().parents[2]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.mask import mask_cobra_model  # noqa: E402
from offgrowth.panel import write_panel_tsv  # noqa: E402
from offgrowth2.runner import (  # noqa: E402
    SCHEMA_VERSION,
    main as runner_main,
    read_gzip_json,
)
from offgrowth2.tests.test_engine import _build_toy_fixture  # noqa: E402


def _tmp_base_dir() -> str | None:
    cand = Path("/home/leechuck/tmp")
    return str(cand) if cand.is_dir() else None


def _write_toy_genome_files(root: Path, genome_id: str = "toy_genome") -> Dict[str, Path]:
    (
        metabolites,
        reactions,
        _,
        panel,
        default_med,
        draft_model,
        filled_model,
        candidates,
    ) = _build_toy_fixture()

    model_dir = root / "models"
    ref_dir = root / "ref"
    model_dir.mkdir(parents=True, exist_ok=True)
    ref_dir.mkdir(parents=True, exist_ok=True)

    draft_xml = model_dir / f"{genome_id}-draft.xml"
    filled_xml = model_dir / f"{genome_id}-filled.xml"
    medium_csv = model_dir / f"{genome_id}-medium.csv"
    trans_tbl = model_dir / f"{genome_id}-Transporter.tbl"

    cobra.io.write_sbml_model(draft_model, str(draft_xml))
    cobra.io.write_sbml_model(filled_model, str(filled_xml))

    with medium_csv.open("w", encoding="utf-8") as fh:
        fh.write("compounds,name,maxFlux\n")
        for ex_id, flux in sorted(default_med.items()):
            cpd_id = ex_id[len("EX_") : -len("_e0")]
            fh.write(f"{cpd_id},{cpd_id},{flux}\n")

    with trans_tbl.open("w", encoding="utf-8") as fh:
        fh.write("qseqid\tstitle\texid\trea\n")
        fh.write("G_TR_XYL\tG_TR_XYL\tEX_cpd00154_e0\trxn90002\n")

    mets_tsv = ref_dir / "seed_metabolites_edited.tsv"
    with mets_tsv.open("w", encoding="utf-8") as fh:
        fh.write("id\tabbreviation\tname\tformula\tcharge\tkeggID\tbiggID\tis_obsolete\taliases\tbiocycID\n")
        for m in metabolites.values():
            fh.write(
                f"{m.seed_id}\t{m.abbreviation}\t{m.name}\t{m.formula}\t"
                f"{m.charge}\t\t\t0\t\t\n"
            )

    rxns_tsv = ref_dir / "seed_reactions_corrected.tsv"
    with rxns_tsv.open("w", encoding="utf-8") as fh:
        fh.write(
            "id\tname\treversibility\tis_transport\tgapseq.status\t"
            "is_obsolete\tstoichiometry\tec_numbers\n"
        )
        for r in reactions.values():
            stoich_str = ";".join(
                f"{t.coeff}:{t.cpd_id}:{'0' if t.compartment == 'c0' else '1'}:0:\"{t.cpd_id}\""
                for t in r.stoichiometry
            )
            fh.write(
                f"{r.rxn_id}\t{r.name}\t{r.reversibility}\t"
                f"{'1' if r.is_transport else '0'}\t{r.gapseq_status}\t"
                f"{'1' if r.is_obsolete else '0'}\t{stoich_str}\t"
                f"{','.join(r.ec_numbers)}\n"
            )

    panel_tsv = ref_dir / "panel.tsv"
    write_panel_tsv(panel, panel_tsv)

    cands_tsv = ref_dir / "candidates.tsv"
    with cands_tsv.open("w", encoding="utf-8") as fh:
        fh.write("protein_id\tec\tscore\n")
        for c in candidates:
            fh.write(f"{c.protein_id}\t{c.ec}\t{c.score}\n")

    labels_tsv = ref_dir / "fitness_browser_labels.tsv"
    with labels_tsv.open("w", encoding="utf-8") as fh:
        fh.write(
            "genome_id\tprotein_id\torg_id\tlocus_id\tfunction_id\t"
            "evidence\ttier\toverlap\thedged\tnew_annotation\tcomment\n"
        )
        # Overlap == price148 row: included in known_labels, NOT a label case
        fh.write(
            f"{genome_id}\tG_GLC\tToyOrg\tL1\t1.1.1.1\t"
            "fitness_browser_reannotation\tcurated_phenotype\tprice148\t0\tGlc\tseen\n"
        )
        # Overlap == new row: becomes case FB:P_XYL:5.3.1.5 with hedged=1
        fh.write(
            f"{genome_id}\tP_XYL\tToyOrg\tL2\t5.3.1.5\t"
            "fitness_browser_reannotation\tcurated_phenotype\tnew\t1\tXyl\tpossible isomerase\n"
        )

    return {
        "model_dir": model_dir,
        "mets_tsv": mets_tsv,
        "rxns_tsv": rxns_tsv,
        "panel_tsv": panel_tsv,
        "cands_tsv": cands_tsv,
        "labels_tsv": labels_tsv,
    }


def _strip_timings(payload: Dict[str, object]) -> Dict[str, object]:
    out = copy.deepcopy(payload)
    out.pop("inputs", None)
    runs = out.get("runs", [])
    if isinstance(runs, list):
        for r in runs:
            if isinstance(r, dict):
                r.pop("cpu_seconds", None)
                r.pop("wall_seconds", None)
                r.pop("peak_rss_kb", None)
                bf = r.get("baseline_fill")
                if isinstance(bf, dict):
                    bf.pop("cpu_seconds", None)
                tests = r.get("tests", [])
                if isinstance(tests, list):
                    for t in tests:
                        if isinstance(t, dict):
                            t.pop("cpu_seconds", None)
    return out


class TestOffGrowth2Runner(unittest.TestCase):
    """End-to-end unit tests for `offgrowth2.runner`."""

    def test_runner_end_to_end_schema_and_determinism(self) -> None:
        with tempfile.TemporaryDirectory(dir=_tmp_base_dir()) as tmp:
            root = Path(tmp)
            paths = _write_toy_genome_files(root, "toy_genome")
            out1 = root / "out1"
            out2 = root / "out2"

            base_args: List[str] = [
                "genome",
                "--genome-id",
                "toy_genome",
                "--genome-model-dir",
                str(paths["model_dir"]),
                "--candidates-tsv",
                str(paths["cands_tsv"]),
                "--candidate-table-kind",
                "six",
                "--seed-reactions",
                str(paths["rxns_tsv"]),
                "--seed-metabolites",
                str(paths["mets_tsv"]),
                "--panel-tsv",
                str(paths["panel_tsv"]),
                "--label-cases-tsv",
                str(paths["labels_tsv"]),
                "--label-cases-prefix",
                "FB:",
                "--workers",
                "1",
            ]

            rc1 = runner_main(base_args + ["--out-dir", str(out1)])
            self.assertEqual(rc1, 0)
            rc2 = runner_main(base_args + ["--out-dir", str(out2)])
            self.assertEqual(rc2, 0)

            work1 = out1 / "toy_genome"
            for expected_file in (
                "genome_results2.json.gz",
                "receipt.json",
                "timing_per_test.tsv",
                "fill_sizes.tsv",
                "runner.log",
            ):
                self.assertTrue((work1 / expected_file).exists(), f"Missing {expected_file}")

            p1 = read_gzip_json(work1 / "genome_results2.json.gz")
            p2 = read_gzip_json(out2 / "toy_genome" / "genome_results2.json.gz")

            # Verify top-level schema keys
            self.assertEqual(p1["schema"], SCHEMA_VERSION)
            self.assertEqual(p1["genome_id"], "toy_genome")
            self.assertEqual(p1["candidate_table_kind"], "six")
            for req_key in (
                "inputs",
                "panel",
                "minimal_medium",
                "default_medium",
                "model_genes",
                "draft_reactions",
                "known_labels",
                "runs",
            ):
                self.assertIn(req_key, p1)

            # known_labels includes both L0 (overlap=price148) and L1 (overlap=new)
            self.assertEqual(
                p1["known_labels"],
                {"G_GLC": ["1.1.1.1"], "P_XYL": ["5.3.1.5"]},
            )

            runs = p1["runs"]
            assert isinstance(runs, list)
            self.assertEqual(len(runs), 2)
            self.assertEqual(runs[0]["run_kind"], "unmasked")
            self.assertIsNone(runs[0]["case"])
            self.assertEqual(runs[1]["run_kind"], "masked")
            self.assertEqual(
                runs[1]["case"],
                {
                    "case_id": "FB:P_XYL:5.3.1.5",
                    "protein_id": "P_XYL",
                    "truth_ecs": ["5.3.1.5"],
                    "case_source": "fitness_browser",
                    "hedged": 1,
                },
            )
            self.assertTrue(runs[1]["audit"]["passed"])
            self.assertEqual(runs[1]["stratum"]["stratum"], "off_growth")
            self.assertEqual(runs[1]["target_reactions"], ["rxn00002"])

            # Determinism across runs ignoring timings
            self.assertEqual(_strip_timings(p1), _strip_timings(p2))

    def test_masking_audit_failure_skips_case_and_records_it(self) -> None:
        with tempfile.TemporaryDirectory(dir=_tmp_base_dir()) as tmp:
            root = Path(tmp)
            paths = _write_toy_genome_files(root, "toy_genome")
            # Add a second case P_LEAK whose masking we force to leak
            with paths["labels_tsv"].open("a", encoding="utf-8") as fh:
                fh.write(
                    "toy_genome\tP_LEAK\tToyOrg\tL3\t1.1.1.27\t"
                    "fitness_browser_reannotation\tcurated_phenotype\tnew\t0\tLac\ttest\n"
                )

            def _leaky_mask_model(draft_model, protein_id, target_reactions):
                if protein_id == "P_LEAK":
                    # Return model with an unremoved empty-GPR target reaction so audit_case_masking fails
                    leaked = copy.deepcopy(draft_model)
                    leaked.reactions.get_by_id("rxn00020_c0").gene_reaction_rule = ""
                    return leaked, (), ()
                return mask_cobra_model(draft_model, protein_id, target_reactions)

            out_dir = root / "out_audit"
            with mock.patch("offgrowth2.engine.mask_cobra_model", side_effect=_leaky_mask_model):
                rc = runner_main(
                    [
                        "genome",
                        "--genome-id",
                        "toy_genome",
                        "--genome-model-dir",
                        str(paths["model_dir"]),
                        "--candidates-tsv",
                        str(paths["cands_tsv"]),
                        "--candidate-table-kind",
                        "six",
                        "--seed-reactions",
                        str(paths["rxns_tsv"]),
                        "--seed-metabolites",
                        str(paths["mets_tsv"]),
                        "--panel-tsv",
                        str(paths["panel_tsv"]),
                        "--label-cases-tsv",
                        str(paths["labels_tsv"]),
                        "--skip-unmasked",
                        "--workers",
                        "1",
                        "--out-dir",
                        str(out_dir),
                    ]
                )
            self.assertEqual(rc, 0)
            payload = read_gzip_json(out_dir / "toy_genome" / "genome_results2.json.gz")
            runs = payload["runs"]
            assert isinstance(runs, list)
            # Only P_XYL ran; P_LEAK was skipped before any LP
            self.assertEqual(len(runs), 1)
            self.assertEqual(runs[0]["case"]["case_id"], "FB:P_XYL:5.3.1.5")

            skipped = payload["skipped_cases"]
            assert isinstance(skipped, list)
            self.assertEqual(len(skipped), 1)
            self.assertEqual(skipped[0]["case"]["case_id"], "FB:P_LEAK:1.1.1.27")
            self.assertEqual(skipped[0]["reason"], "masking_audit_failed")
            self.assertFalse(skipped[0]["audit"]["passed"])
            self.assertTrue(len(skipped[0]["audit"]["violations"]) > 0)

    def test_resume_skips_finished_cases(self) -> None:
        with tempfile.TemporaryDirectory(dir=_tmp_base_dir()) as tmp:
            root = Path(tmp)
            paths = _write_toy_genome_files(root, "toy_genome")
            out_dir = root / "out_resume"

            args = [
                "genome",
                "--genome-id",
                "toy_genome",
                "--genome-model-dir",
                str(paths["model_dir"]),
                "--candidates-tsv",
                str(paths["cands_tsv"]),
                "--candidate-table-kind",
                "clean",
                "--seed-reactions",
                str(paths["rxns_tsv"]),
                "--seed-metabolites",
                str(paths["mets_tsv"]),
                "--panel-tsv",
                str(paths["panel_tsv"]),
                "--label-cases-tsv",
                str(paths["labels_tsv"]),
                "--workers",
                "1",
                "--out-dir",
                str(out_dir),
            ]
            # First run completes unmasked + FB:P_XYL:5.3.1.5
            self.assertEqual(runner_main(args), 0)

            # Second run with --resume must not call execute_masked_case_superset or run_superset_for_model
            with (
                mock.patch(
                    "offgrowth2.runner.execute_masked_case_superset",
                    side_effect=AssertionError("Finished case should not be re-executed on --resume"),
                ) as mock_case,
                mock.patch(
                    "offgrowth2.runner.run_superset_for_model",
                    side_effect=AssertionError("Finished unmasked run should not be re-executed on --resume"),
                ) as mock_unmasked,
            ):
                self.assertEqual(runner_main(args + ["--resume"]), 0)
                self.assertEqual(mock_case.call_count, 0)
                self.assertEqual(mock_unmasked.call_count, 0)

            payload = read_gzip_json(out_dir / "toy_genome" / "genome_results2.json.gz")
            self.assertEqual(len(payload["runs"]), 2)


if __name__ == "__main__":
    unittest.main()
