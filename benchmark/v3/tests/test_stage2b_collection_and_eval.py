#!/usr/bin/env python3
"""Unit tests for Stage 2B Price collection runner and evaluator."""

from __future__ import annotations

import contextlib
import csv
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

import ec_identity as E  # noqa: E402
from stage2b import evaluate_price_stage2b as ES2B  # noqa: E402
from stage2b import run_price_stage2b_collection as S2B  # noqa: E402


def sample_ec_entries() -> dict[str, dict[str, object]]:
    return {
        "1.1.1.1": dict(status="transferred", targets=["1.1.1.2"]),
        "1.1.1.2": dict(status="active", targets=[]),
        "1.1.1.3": dict(status="transferred", targets=["1.1.1.2", "1.1.1.4"]),
        "1.1.1.4": dict(status="active", targets=[]),
        "1.1.1.5": dict(status="deleted", targets=[]),
        "1.1.1.n1": dict(status="active", targets=[]),
    }


def write_stage_ledger(path: Path, genome_id: str, arm: str, gram: str = "neg", medium: str = "minimal") -> None:
    if arm == "stage0":
        met_stage = dict(name="metabolic_refinement", status="skipped", reason="not_enabled")
    else:
        met_stage = dict(
            name="metabolic_refinement",
            status="executed",
            gram=gram,
            medium=medium,
            n_loci=2,
            n_repaired_reactions=3,
        )
    payload = dict(
        schema_version=2,
        genome_id=genome_id,
        stages=[
            dict(name="evidence_integration", status="executed"),
            met_stage,
        ],
    )
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_metabolism_artifacts(output_tsv: Path) -> dict[str, str]:
    met_dir = Path(str(output_tsv) + ".metabolism")
    met_dir.mkdir(parents=True, exist_ok=True)
    (met_dir / "fba_input.json").write_text(
        json.dumps({"gram": "neg", "medium": "minimal", "loci": {}}, indent=2) + "\n",
        encoding="utf-8",
    )
    (met_dir / "fba_output.json").write_text(
        json.dumps(
            {
                "max_univ_bio": 33.99236,
                "target_bio": 3.39924,
                "repaired_reactions": ["rxn00001"],
                "loci": {},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (met_dir / "metabolic_model.xml").write_text(
        '<?xml version="1.0" encoding="UTF-8"?><sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1"/>\n',
        encoding="utf-8",
    )
    (met_dir / "reactions.tsv").write_text(
        "reaction_id\tname\tequation\tec_numbers\tstatus\tlb\tub\tflux_biomass\trxn_prob\tbest_locus\tbest_prob\n",
        encoding="utf-8",
    )
    files = [
        met_dir / "fba_input.json",
        met_dir / "fba_output.json",
        met_dir / "metabolic_model.xml",
        met_dir / "reactions.tsv",
    ]
    return S2B.V.file_hashes(files)


class Stage2BCollectionAndEvaluationTests(unittest.TestCase):
    def setUp(self) -> None:
        base_tmp = Path(os.environ.get("TMPDIR", str(Path.home() / "tmp")))
        base_tmp.mkdir(parents=True, exist_ok=True)
        self._tmp = tempfile.TemporaryDirectory(dir=base_tmp, prefix="test-stage2b-eval-")
        self.root = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_verify_price22_gram_negative_pseudomonadota(self) -> None:
        assemblies = {
            asm: {"taxid": info["taxid"], "proteins": 100}
            for asm, info in S2B.PRICE22_TAXONOMY.items()
        }
        res = S2B.verify_price22_gram_negative_pseudomonadota({"assemblies": assemblies})
        self.assertTrue(res["all_gram_negative"])
        self.assertTrue(res["all_pseudomonadota"])
        self.assertEqual(res["phylum"], "Pseudomonadota")
        self.assertEqual(res["phylum_taxid"], 1224)
        self.assertEqual(res["gram"], "neg")
        self.assertEqual(res["n_genomes"], 22)

        bad = dict(assemblies)
        bad["GCA_000006965.1"] = {"taxid": 999999, "proteins": 100}
        with self.assertRaisesRegex(ValueError, "taxid mismatch"):
            S2B.verify_price22_gram_negative_pseudomonadota({"assemblies": bad})

        subset = {"GCA_000006965.1": assemblies["GCA_000006965.1"]}
        with self.assertRaisesRegex(ValueError, "do not match the 22 Price genomes"):
            S2B.verify_price22_gram_negative_pseudomonadota({"assemblies": subset})

    def test_verify_arm_command_pair_and_comparable_specs(self) -> None:
        fba_dir = self.root / "fba_data"
        fba_dir.mkdir()
        sidecar = self.root / "stoichiometric_locus_refiner.py"
        sidecar.write_text("# sidecar\n", encoding="utf-8")
        py_bin = self.root / "python"
        py_bin.write_text("# py\n", encoding="utf-8")

        base_s0 = [
            "java",
            "-Djava.io.tmpdir=/tmp/s0",
            "-jar",
            "/path/gspa.jar",
            "integrate",
            "--batch-manifest",
            str(self.root / "production/stage0/manifest.tsv"),
        ]
        base_s2b = [
            "java",
            "-Djava.io.tmpdir=/tmp/s2b",
            "-jar",
            "/path/gspa.jar",
            "integrate",
            "--batch-manifest",
            str(self.root / "production/stage2b/manifest.tsv"),
        ]
        cmd_s0 = base_s0 + ["--enable-priors", "none"]
        cmd_s2b = base_s2b + S2B.stage2b_extra_flags(
            fba_data_dir=fba_dir,
            fba_gram="neg",
            fba_medium="minimal",
            fba_python_bin=py_bin,
            fba_sidecar=sidecar,
        )

        S2B.verify_arm_command_pair(
            cmd_s0,
            cmd_s2b,
            fba_data_dir=str(fba_dir),
            fba_gram="neg",
            fba_medium="minimal",
            fba_python_bin=str(py_bin),
            fba_sidecar=str(sidecar),
        )

        # Reject any undeclared flag difference
        with self.assertRaisesRegex(ValueError, "stage2b command must end with documented Stage 2B flags"):
            S2B.verify_arm_command_pair(
                cmd_s0,
                cmd_s2b + ["--unexpected-flag"],
                fba_data_dir=str(fba_dir),
                fba_gram="neg",
                fba_medium="minimal",
                fba_python_bin=str(py_bin),
                fba_sidecar=str(sidecar),
            )
        with self.assertRaisesRegex(ValueError, "Disallowed command difference"):
            bad_s0 = list(cmd_s0)
            bad_s0[3] = "/other/gspa.jar"
            S2B.verify_arm_command_pair(
                bad_s0,
                cmd_s2b,
                fba_data_dir=str(fba_dir),
                fba_gram="neg",
                fba_medium="minimal",
                fba_python_bin=str(py_bin),
                fba_sidecar=str(sidecar),
            )

    def test_collection_run_and_finalize_verifies_stage2b_metabolism_artifacts(self) -> None:
        col = self.root / "collection"
        col.mkdir()
        fba_dir = self.root / "fba_data"
        fba_dir.mkdir()
        sidecar = self.root / "stoichiometric_locus_refiner.py"
        sidecar.write_text("# sidecar\n", encoding="utf-8")
        py_bin = Path(sys.executable).resolve()

        cell = dict(index=0, part="/part", claims=1, assembly="GCA_000006965.1", predictor="topec", matrix_rows=1)
        plan = dict(
            schema="gspa.price-stage2b-production/1",
            cells=[cell],
            arms=["stage0", "stage2b"],
            layout="/layout",
            jar="/jar",
            reference="/ref",
            fba_data_dir=str(fba_dir),
            fba_gram="neg",
            fba_medium="minimal",
            fba_python_bin=str(py_bin),
            fba_sidecar=str(sidecar),
            fba_sidecar_resolved=str(sidecar),
            taxonomy_verification={"all_gram_negative": True, "all_pseudomonadota": True},
            expected_claims=1,
        )
        S2B.save(col / "plan.json", plan)
        S2B.save(
            col / "registration.json",
            dict(inputs={}, plan_sha256=S2B.V.sha256_file(str(col / "plan.json"))),
        )

        def mock_capacity(part, layout, jar, reference, fba_data_dir, **kwargs):
            out = kwargs["out"]
            out.mkdir(parents=True, exist_ok=False)
            runs = {}
            for arm in kwargs["arms"]:
                d = out / arm
                d.mkdir()
                tsv = d / "GCA_000006965.1.integrated.tsv"
                tsv.write_text(
                    "protein_id\ttype\tfunction_id\taspect\tposterior_prob\np1\tEC\t1.1.1.1\t\t0.8\n",
                    encoding="utf-8",
                )
                write_stage_ledger(
                    Path(str(tsv) + ".stages.json"),
                    "GCA_000006965.1",
                    arm,
                    gram=kwargs["fba_gram"],
                    medium=kwargs["fba_medium"],
                )
                met_hashes = {}
                if arm == "stage0":
                    cmd = ["java", f"--out={d}/out.tsv", "--enable-priors", "none"]
                else:
                    cmd = ["java", f"--out={d}/out.tsv"] + S2B.stage2b_extra_flags(
                        fba_data_dir=fba_data_dir,
                        fba_gram=kwargs["fba_gram"],
                        fba_medium=kwargs["fba_medium"],
                        fba_python_bin=kwargs["fba_python_bin"],
                        fba_sidecar=kwargs["fba_sidecar"],
                    )
                    met_hashes = write_metabolism_artifacts(tsv)
                spec = dict(
                    outputs={"GCA_000006965.1": str(tsv)},
                    required_stages=["evidence_integration"],
                    inputs={},
                    command=cmd,
                )
                receipt = d / S2B.V.RECEIPT
                S2B.save(receipt, dict(spec=spec, returncode=0, artifacts=S2B.V.output_hashes(spec)))
                runs[arm] = dict(
                    prediction_rows=1,
                    input_candidates=1,
                    candidate_universe_preserved=True,
                    receipt_sha256=S2B.V.sha256_file(str(receipt)),
                    metabolism_artifacts=met_hashes,
                )
            res = dict(
                schema="gspa.price-stage2b-capacity/1",
                status="verified",
                runs=runs,
                assembly="GCA_000006965.1",
                predictor="topec",
                matrix_rows=1,
            )
            S2B.save(out / "verification.json", res)
            return res

        with patch.object(S2B, "run_cell_capacity", side_effect=mock_capacity):
            S2B.run(col, 0, "job1", "24g", 100, self.root / "scratch")

        attempt = col / "attempts/0/job1"
        self.assertTrue((attempt / "completed.json").is_file())
        S2B.finalize(col)
        ver = json.loads((col / "collection-verification.json").read_text(encoding="utf-8"))
        self.assertEqual(ver["status"], "verified")
        self.assertEqual(ver["cells"], 1)
        self.assertEqual(ver["arms"], 2)

        # Tampering with a Stage 2B metabolism artifact must fail finalization
        (col / "collection-verification.json").unlink()
        met_xml = attempt / "production/stage2b/GCA_000006965.1.integrated.tsv.metabolism/metabolic_model.xml"
        met_xml.write_text(met_xml.read_text(encoding="utf-8") + "<!-- tampered -->\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Stage 2B metabolism artifacts changed"):
            S2B.finalize(col)

    def _build_full_price148_panel_fixture(self, *, noop_stage2b: bool = False) -> tuple[Path, Path, Path]:
        """Build a 132-cell registered collection with 22 completed topec cells on the 148 Price references."""
        manifest = Path(__file__).resolve().parents[2] / "manifests/suite6_price149_manifest.tsv"
        with manifest.open() as f:
            genomes = sorted({r["host_assembly"] for r in csv.DictReader(f, delimiter="\t")})
        refs, excluded = ES2B.R.references(manifest, set(genomes))
        self.assertEqual((len(refs), len(excluded), sum(map(len, refs.values()))), (148, 1, 151))

        tag = "noop" if noop_stage2b else "active"
        sub = self.root / tag
        sub.mkdir()
        collection = sub / "collection"
        collection.mkdir()
        fba_dir = sub / "fba_data"
        fba_dir.mkdir()
        sidecar = sub / "stoichiometric_locus_refiner.py"
        sidecar.write_text("# sidecar\n", encoding="utf-8")
        py_bin = Path(sys.executable).resolve()

        audit = sub / "audit.json"
        ent = sample_ec_entries()
        ES2B.save(
            audit,
            dict(
                schema="gspa.ec-identity-audit/1",
                inputs={"fixture": "fixture"},
                entries=ent,
                resolutions={x: E.resolve(ent, x) for x in ent},
            ),
        )

        layout = sub / "layout"
        layout.mkdir()
        inventories = {}
        for g in genomes:
            p = layout / (g + ".sequence_hashes.tsv")
            ids = [pid for gid, pid in refs if gid == g] + ["unlabelled"]
            p.write_text(
                "protein_id\tlength\tsequence_sha256\n" + "".join(f"{pid}\t1\tfixture\n" for pid in ids),
                encoding="utf-8",
            )
            inventories[g] = dict(
                taxid=S2B.PRICE22_TAXONOMY[g]["taxid"],
                proteins=len(ids),
                sequence_hashes_sha256=ES2B.sha(p),
            )
        ES2B.save(layout / "report.json", dict(assemblies=inventories))

        cells_list = []
        for pred in ("topec", "mapred", "enzbert", "dpz", "graphec", "clean"):
            for g in genomes:
                count = sum(len(v) + 1 for (gid, p), v in refs.items() if gid == g) + 1
                cells_list.append(
                    dict(
                        index=len(cells_list),
                        predictor=pred,
                        assembly=g,
                        matrix_rows=100,
                        claims=count,
                        part=str(sub / "parts" / pred / g),
                    )
                )

        plan = dict(
            schema="gspa.price-stage2b-production/1",
            arms=["stage0", "stage2b"],
            cells=cells_list,
            layout=str(layout),
            fba_data_dir=str(fba_dir),
            fba_gram="neg",
            fba_medium="minimal",
            fba_python_bin=str(py_bin),
            fba_sidecar=str(sidecar),
            fba_sidecar_resolved=str(sidecar),
            taxonomy_verification=S2B.verify_price22_gram_negative_pseudomonadota({"assemblies": inventories}),
        )
        ES2B.save(collection / "plan.json", plan)
        ES2B.save(
            collection / "registration.json",
            dict(
                inputs={
                    str(layout / "report.json"): ES2B.sha(layout / "report.json"),
                    str(sidecar): ES2B.sha(sidecar),
                    str(py_bin): ES2B.sha(py_bin),
                },
                plan_sha256=ES2B.sha(collection / "plan.json"),
            ),
        )

        for cell in cells_list[:22]:
            g = cell["assembly"]
            part = Path(cell["part"])
            part.mkdir(parents=True)
            claims = [
                dict(protein_id=p, function_type="EC", function_id=term, raw_score=0.2, source="topec")
                for (gid, p), terms in refs.items()
                if gid == g
                for term in sorted(terms)
            ]
            claims += [
                dict(protein_id=p, function_type="EC", function_id="9.9.9.9", raw_score=0.3, source="topec")
                for gid, p in refs
                if gid == g
            ]
            claims += [
                dict(
                    protein_id="unlabelled",
                    function_type="EC",
                    function_id="9.9.9.9",
                    raw_score=0.3,
                    source="topec",
                )
            ]
            raw = part / "claims.jsonl"
            raw.write_text("".join(json.dumps(x) + "\n" for x in claims), encoding="utf-8")

            attempt_dir = collection / "attempts" / str(cell["index"]) / "fixture"
            attempt_dir.mkdir(parents=True)
            ES2B.save(
                attempt_dir / "started.json",
                dict(cell=cell, registration_sha256=ES2B.sha(collection / "registration.json")),
            )
            runs = {}
            for arm in ("stage0", "stage2b"):
                d = attempt_dir / "production" / arm
                d.mkdir(parents=True)
                output = d / "output.tsv"
                use_stage2b_boost = (arm == "stage2b") and (not noop_stage2b)
                lines = []
                for x in claims:
                    if x["function_id"] == "9.9.9.9":
                        prob = "0.7" if not use_stage2b_boost else "0.3"
                    else:
                        prob = "0.4" if not use_stage2b_boost else "0.85"
                    lines.append(f"{x['protein_id']}\tEC\t{x['function_id']}\t{prob}\n")
                output.write_text(
                    "protein_id\ttype\tfunction_id\tposterior_prob\n" + "".join(lines),
                    encoding="utf-8",
                )
                ledger = d / "output.tsv.stages.json"
                write_stage_ledger(ledger, g, arm, gram="neg", medium="minimal")
                met_hashes = {}
                base_cmd = ["java", f"-Djava.io.tmpdir={d}/tmp", "integrate", f"--out={output}"]
                spec_inputs = {str(raw): ES2B.sha(raw)}
                if arm == "stage0":
                    command = base_cmd + ["--enable-priors", "none"]
                else:
                    command = base_cmd + S2B.stage2b_extra_flags(
                        fba_data_dir=fba_dir,
                        fba_gram="neg",
                        fba_medium="minimal",
                        fba_python_bin=py_bin,
                        fba_sidecar=sidecar,
                    )
                    met_hashes = write_metabolism_artifacts(output)
                    spec_inputs[str(sidecar)] = ES2B.sha(sidecar)
                    spec_inputs[str(py_bin)] = ES2B.sha(py_bin)

                ES2B.save(
                    d / "verified_run.json",
                    dict(
                        returncode=0,
                        spec=dict(command=command, outputs={g: str(output)}, inputs=spec_inputs),
                        artifacts={str(output): ES2B.sha(output), str(ledger): ES2B.sha(ledger)},
                    ),
                )
                runs[arm] = dict(
                    candidate_universe_preserved=True,
                    input_candidates=len(claims),
                    prediction_rows=len(claims),
                    receipt_sha256=ES2B.sha(d / "verified_run.json"),
                    metabolism_artifacts=met_hashes,
                )
            verification = attempt_dir / "production/verification.json"
            ES2B.save(
                verification,
                dict(status="verified", assembly=g, predictor="topec", matrix_rows=100, runs=runs),
            )
            ES2B.save(
                attempt_dir / "completed.json",
                dict(status="verified", cell=cell, verification_sha256=ES2B.sha(verification)),
            )

        return collection, manifest, audit

    def test_evaluate_price_stage2b_scores_148_references_and_reports_fmax_and_corrections(self) -> None:
        collection, manifest, audit = self._build_full_price148_panel_fixture(noop_stage2b=False)
        with contextlib.redirect_stdout(io.StringIO()):
            report = ES2B.evaluate(
                collection,
                manifest,
                audit,
                ["topec"],
                self.root / "eval-active",
                n_boot=100,
                seed=13,
            )
        self.assertEqual(report["status"], "verified")
        self.assertEqual(report["cells"], 22)
        self.assertEqual(report["reference_proteins"], 148)
        self.assertEqual(report["original_reference_pairs"], 151)
        self.assertEqual(len(report["exclusions"]), 1)
        self.assertEqual(report["exclusions"][0]["price_entry"], "NP_384885")

        for mode in ("exact", "current_unambiguous"):
            panel = report["summaries"][mode]["topec"]
            self.assertEqual(panel["arms"]["stage0"]["top1_match"], 0)
            self.assertEqual(panel["arms"]["stage2b"]["top1_match"], 148)
            self.assertGreater(
                panel["arms"]["stage2b"]["micro_fmax"]["fmax"],
                panel["arms"]["stage0"]["micro_fmax"]["fmax"],
            )
            self.assertGreater(
                panel["arms"]["stage2b"]["cafa_fmax"]["fmax"],
                panel["arms"]["stage0"]["cafa_fmax"]["fmax"],
            )
            ch = panel["paired"]["stage0_to_stage2b"]["corrections_and_harms"]
            self.assertEqual(ch["lexical_top1"]["corrected"], 148)
            self.assertEqual(ch["lexical_top1"]["harmed"], 0)
            self.assertEqual(ch["lexical_top1"]["net_gain"], 148)
            self.assertGreater(ch["pair_at_0_5"]["corrected"], 0)
            self.assertEqual(ch["pair_at_0_5"]["harmed"], 0)

    def test_noop_arm_stage0_vs_itself_gives_exactly_zero_corrections_and_harms(self) -> None:
        collection, manifest, audit = self._build_full_price148_panel_fixture(noop_stage2b=True)
        with contextlib.redirect_stdout(io.StringIO()):
            report = ES2B.evaluate(
                collection,
                manifest,
                audit,
                ["topec"],
                self.root / "eval-noop",
                n_boot=100,
                seed=13,
            )

        for mode in ("exact", "current_unambiguous"):
            panel = report["summaries"][mode]["topec"]
            for pair_key in ("stage0_to_stage0", "stage0_to_stage2b"):
                paired = panel["paired"][pair_key]
                self.assertEqual(paired["top1_credit_increase"], 0)
                self.assertEqual(paired["top1_credit_decrease"], 0)
                self.assertEqual(paired["top1_gained_reference_match"], 0)
                self.assertEqual(paired["top1_lost_reference_match"], 0)
                self.assertEqual(paired["changed_top1"], 0)

                for boot_key in ("top1_delta_bootstrap", "micro_fmax_delta_bootstrap", "cafa_fmax_delta_bootstrap"):
                    b = paired[boot_key]
                    self.assertEqual((b["estimate"], b["ci_low"], b["ci_high"]), (0.0, 0.0, 0.0))

                ch = paired["corrections_and_harms"]
                for level in ("pair_at_0_5", "lexical_top1"):
                    self.assertEqual(ch[level]["corrected"], 0)
                    self.assertEqual(ch[level]["harmed"], 0)
                    self.assertEqual(ch[level]["net_gain"], 0)
                    self.assertEqual(ch[level]["mcnemar_p"], 1.0)
                    for bname in (
                        "corrected_genome_bootstrap",
                        "harmed_genome_bootstrap",
                        "net_gain_genome_bootstrap",
                    ):
                        b = ch[level][bname]
                        self.assertEqual((b["estimate"], b["ci_low"], b["ci_high"]), (0.0, 0.0, 0.0))

                tie = ch["tie_aware_top1"]
                self.assertEqual(tie["top1_credit_increase"], 0)
                self.assertEqual(tie["top1_credit_decrease"], 0)
                self.assertEqual(tie["net_credit_gain"], 0)


if __name__ == "__main__":
    unittest.main()
