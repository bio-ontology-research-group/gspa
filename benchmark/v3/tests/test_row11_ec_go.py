import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

TMP_ROOT = Path.home() / "tmp"
TMP_ROOT.mkdir(parents=True, exist_ok=True)
os.environ["TMPDIR"] = str(TMP_ROOT)
tempfile.tempdir = str(TMP_ROOT)

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))

import v3_runner as runner  # noqa: E402
from row11_ec_go import compare_full_vs_filtered as C  # noqa: E402
from row11_ec_go import derive_operons as DO  # noqa: E402
from row11_ec_go import merge_claims as MC  # noqa: E402
from row11_ec_go import run_arms as RA  # noqa: E402


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_part(parts_root: Path, predictor: str, assembly: str, lines: list[str], bad_sha: bool = False) -> str:
    pdir = parts_root / predictor / assembly
    pdir.mkdir(parents=True, exist_ok=True)
    payload = "".join(ln if ln.endswith("\n") else ln + "\n" for ln in lines).encode("utf-8")
    (pdir / "claims.jsonl").write_bytes(payload)
    digest = "0" * 64 if bad_sha else _sha256_bytes(payload)
    part = {
        "assembly": assembly,
        "predictor": predictor,
        "claims": {
            "path": "claims.jsonl",
            "sha256": digest,
            "claims": len(lines),
            "bytes": len(payload),
        },
    }
    (pdir / "part.json").write_text(json.dumps(part), encoding="utf-8")
    return _sha256_bytes(payload)


class TestMergeClaims(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=str(TMP_ROOT))
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_merge_keeps_min_raw_no_raw_score_and_extra_in_order(self):
        parts_root = self.root / "parts"
        asm = "GCF_000015245.1"
        p1_lines = [
            json.dumps({"protein_id": "P1", "function_id": "1.1.1.1", "function_type": "EC", "raw_score": 0.01, "source": "clean"}),
            json.dumps({"protein_id": "P2", "function_id": "1.1.1.2", "function_type": "EC", "raw_score": 0.05, "source": "clean"}),
            json.dumps({"protein_id": "P3", "function_id": "1.1.1.3", "function_type": "EC", "source": "clean"}),
        ]
        p2_lines = [
            json.dumps({"protein_id": "P1", "function_id": "2.7.1.1", "function_type": "EC", "raw_score": 0.80, "source": "topec"}),
            json.dumps({"protein_id": "P4", "function_id": "2.7.1.2", "function_type": "EC", "raw_score": 0.049, "source": "topec"}),
            json.dumps({"protein_id": "P5", "function_id": "2.7.1.3", "function_type": "EC", "raw_score": None, "source": "topec"}),
        ]
        sha1 = _write_part(parts_root, "clean", asm, p1_lines)
        sha2 = _write_part(parts_root, "topec", asm, p2_lines)

        extra_path = self.root / "extra.claims.jsonl"
        extra_lines = [
            json.dumps({"protein_id": "P9", "function_id": "3.1.1.1", "function_type": "EC", "raw_score": 0.001, "source": "swissprot_diamond_baseline"}),
            json.dumps({"protein_id": "P9", "function_id": "GO:0008150", "function_type": "GO", "raw_score": 0.42, "source": "swissprot_diamond_baseline"}),
        ]
        extra_bytes = ("".join(ln + "\n" for ln in extra_lines)).encode("utf-8")
        extra_path.write_bytes(extra_bytes)
        extra_sha = _sha256_bytes(extra_bytes)

        out_dir = self.root / "out"
        rec = MC.merge_claims(
            parts_root=parts_root,
            predictors=["clean", "topec"],
            assembly=asm,
            min_raw=0.05,
            out_dir=out_dir,
            extra_claims=extra_path,
            extra_sha256=extra_sha,
            extra_name="swissprot_diamond_baseline",
        )

        self.assertEqual(rec["assembly"], asm)
        self.assertEqual(rec["min_raw"], 0.05)
        self.assertEqual(
            rec["predictors"]["clean"],
            {"claims_in": 3, "claims_kept": 2, "no_raw_score": 1, "sha256_in": sha1, "verified": True},
        )
        self.assertEqual(
            rec["predictors"]["topec"],
            {"claims_in": 3, "claims_kept": 2, "no_raw_score": 1, "sha256_in": sha2, "verified": True},
        )
        self.assertEqual(
            rec["extra"],
            {"name": "swissprot_diamond_baseline", "path": str(extra_path), "sha256": extra_sha, "claims": 2},
        )
        self.assertEqual(rec["claims"]["claims"], 6)
        merged_lines = [json.loads(ln) for ln in (out_dir / "claims.jsonl").read_text(encoding="utf-8").splitlines()]
        self.assertEqual([r["protein_id"] for r in merged_lines], ["P2", "P3", "P1", "P5", "P9", "P9"])
        self.assertEqual(_sha256_bytes((out_dir / "claims.jsonl").read_bytes()), rec["claims"]["sha256"])

    def test_merge_sha256_mismatch_fails_cleanly(self):
        parts_root = self.root / "parts_bad"
        asm = "GCF_000015245.1"
        _write_part(
            parts_root,
            "clean",
            asm,
            [json.dumps({"protein_id": "P1", "function_id": "1.1.1.1", "function_type": "EC", "raw_score": 0.9, "source": "clean"})],
            bad_sha=True,
        )
        out_dir = self.root / "out_bad"
        with self.assertRaises(ValueError):
            MC.merge_claims(parts_root, ["clean"], asm, 0.05, out_dir)
        self.assertFalse((out_dir / "claims.jsonl").exists())
        self.assertFalse((out_dir / "part.json").exists())

    def test_merge_extra_sha256_mismatch_and_fallback_parser(self):
        parts_root = self.root / "parts_fb"
        asm = "GCF_000015245.1"
        # Include a row where "raw_score" also appears in metadata to trigger json.loads fallback
        lines = [
            json.dumps(
                {
                    "protein_id": "P1",
                    "function_id": "1.1.1.1",
                    "function_type": "EC",
                    "metadata": {"note": '"raw_score": 0.99'},
                    "raw_score": 0.01,
                    "source": "clean",
                }
            ),
            json.dumps(
                {
                    "protein_id": "P2",
                    "function_id": "1.1.1.2",
                    "function_type": "EC",
                    "metadata": {"note": '"raw_score": 0.01'},
                    "raw_score": 0.25,
                    "source": "clean",
                }
            ),
        ]
        _write_part(parts_root, "clean", asm, lines)
        out_dir = self.root / "out_fb"
        rec = MC.merge_claims(parts_root, ["clean"], asm, 0.05, out_dir)
        self.assertEqual(rec["predictors"]["clean"]["claims_in"], 2)
        self.assertEqual(rec["predictors"]["clean"]["claims_kept"], 1)
        self.assertIsNone(rec["extra"])

        extra_path = self.root / "bad_extra.jsonl"
        extra_path.write_text('{"protein_id":"PX","function_id":"1.1.1.1","function_type":"EC","raw_score":0.5,"source":"sp"}\n', encoding="utf-8")
        with self.assertRaises(ValueError):
            MC.merge_claims(
                parts_root,
                ["clean"],
                asm,
                0.05,
                self.root / "out_bad_extra",
                extra_claims=extra_path,
                extra_sha256="f" * 64,
            )


class TestDeriveOperons(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=str(TMP_ROOT))
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_derive_operons_two_operons_and_singleton_excluded(self):
        layouts = self.root / "layouts"
        layouts.mkdir()
        # Two operons of sizes 2 (gap 100 bp) and 3 (gaps 100 bp, 50 bp), plus one singleton on c2 excluded
        layout_6 = (
            "protein_id\tcontig\tstart\tend\tstrand\n"
            "P1\tc1\t100\t400\t+\n"
            "P2\tc1\t501\t800\t+\n"      # gap = 501 - 400 - 1 = 100 bp -> op 1 (size 2)
            "P3\tc1\t1500\t1800\t+\n"    # gap = 699 bp -> break, start op 2
            "P4\tc1\t1901\t2200\t+\n"    # gap = 1901 - 1800 - 1 = 100 bp -> join op 2
            "P5\tc1\t2251\t2500\t+\n"    # gap = 2251 - 2200 - 1 = 50 bp -> join op 2 (size 3)
            "P6\tc2\t100\t400\t-\n"      # singleton -> excluded
        )
        # 5-locus layout fixture (operon size 2 with gap 100 bp, operon size 2 with gap 50 bp, 1 singleton excluded)
        layout_5 = (
            "protein_id\tcontig\tstart\tend\tstrand\n"
            "A1\tc1\t100\t400\t+\n"
            "A2\tc1\t501\t800\t+\n"      # gap 100 bp -> op 1 (size 2)
            "B1\tc1\t1300\t1600\t+\n"    # gap 499 bp -> break
            "B2\tc1\t1651\t1900\t+\n"    # gap 50 bp -> op 2 (size 2)
            "S1\tc1\t2500\t2800\t-\n"    # singleton -> excluded
        )
        (layouts / "G1.layout.tsv").write_text(layout_6, encoding="utf-8")
        (layouts / "G2.layout.tsv").write_text(layout_5, encoding="utf-8")

        out_dir = self.root / "operons"
        report = DO.derive_operons(layouts, out_dir)

        self.assertEqual((out_dir / "G1.operons.tsv").read_text(encoding="utf-8"), "P1\tP2\nP3\tP4\tP5\n")
        self.assertEqual(report["G1"]["n_loci"], 6)
        self.assertEqual(report["G1"]["n_operons"], 2)
        self.assertEqual(report["G1"]["n_proteins_in_operons"], 5)

        self.assertEqual((out_dir / "G2.operons.tsv").read_text(encoding="utf-8"), "A1\tA2\nB1\tB2\n")
        self.assertEqual(report["G2"]["n_loci"], 5)
        self.assertEqual(report["G2"]["n_operons"], 2)
        self.assertEqual(report["G2"]["n_proteins_in_operons"], 4)
        self.assertTrue((out_dir / "operons_report.json").is_file())


class TestRunArmsDry(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=str(TMP_ROOT))
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_genome_rows_and_manifest_carry_operons_and_ec_aware_args(self):
        asm = "GCF_000015245.1"
        parts = self.root / "parts"
        _write_part(
            parts,
            "ec6",
            asm,
            [json.dumps({"protein_id": "P1", "function_id": "1.1.1.1", "function_type": "EC", "raw_score": 0.8, "source": "clean"})],
        )
        layouts = self.root / "layouts"
        layouts.mkdir()
        (layouts / "report.json").write_text(
            json.dumps({"assemblies": {asm: {"taxid": 326297}}}),
            encoding="utf-8",
        )
        (layouts / f"{asm}.layout.tsv").write_text(
            "protein_id\tcontig\tstart\tend\tstrand\nP1\tc1\t1\t100\t+\n",
            encoding="utf-8",
        )
        operons = self.root / "operons"
        operons.mkdir()
        operon_file = operons / f"{asm}.operons.tsv"
        operon_file.write_text("P1\tP2\n", encoding="utf-8")

        meta_path = self.root / "cohort_metadata.json"
        meta_path.write_text(
            json.dumps({"genomes": {asm: {"clade": "Shewanella"}}}),
            encoding="utf-8",
        )

        split = RA.genus_blocked_split(meta_path)
        genomes = RA.genomes_for("ec6", parts, layouts, operons, split)
        self.assertEqual(len(genomes), 1)
        self.assertEqual(genomes[0]["operons"], str(operon_file))

        jar = self.root / "gspa-3.0.0.jar"
        jar.write_bytes(b"jar-placeholder")
        ref = self.root / "reference"
        ref.mkdir()
        models = self.root / "models"
        models.mkdir()
        man = RA.manifest(
            jar=jar,
            reference=ref,
            models=models,
            anomalies=self.root / "anomalies.jsonl",
            threads=8,
            heap="48g",
            anomaly_floor=0.20,
            anomaly_max_candidates=10,
        )
        shared = man["gspa"]["shared_args"]
        self.assertIn("--anomaly-ec-aware", shared)
        idx_floor = shared.index("--anomaly-floor")
        self.assertEqual(shared[idx_floor + 1], "0.20")
        idx_max = shared.index("--anomaly-max-candidates")
        self.assertEqual(shared[idx_max + 1], "10")
        self.assertEqual(shared[idx_max + 2 : idx_max + 4], ["--anomaly-max-locus", "150"])
        self.assertIn(str(Path(RA.__file__).resolve()), man["gspa"]["input_files"])

        man_explicit = RA.manifest(
            jar=jar,
            reference=ref,
            models=models,
            anomalies=self.root / "anomalies.jsonl",
            threads=8,
            heap="48g",
            anomaly_floor=0.20,
            anomaly_max_candidates=10,
            anomaly_max_locus=1000000,
        )
        shared_explicit = man_explicit["gspa"]["shared_args"]
        idx_max_explicit = shared_explicit.index("--anomaly-max-candidates")
        self.assertEqual(shared_explicit[idx_max_explicit + 2 : idx_max_explicit + 4], ["--anomaly-max-locus", "1000000"])

        arm = RA.arms_for(man)[0]
        arm_out = self.root / "stage1"
        rec = runner.run_arm(man, arm, str(arm_out), genomes, dry_run=True)
        self.assertEqual(rec["status"], "dry_run")
        bm_lines = (arm_out / "batch_manifest.tsv").read_text(encoding="utf-8").splitlines()
        headers = bm_lines[0].split("\t")
        values = bm_lines[1].split("\t")
        row_map = dict(zip(headers, values))
        self.assertEqual(row_map["operons_path"], os.path.abspath(str(operon_file)))


class TestCompareFullVsFiltered(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=str(TMP_ROOT))
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_compare_script_on_two_tiny_tsvs(self):
        claims_full = self.root / "claims_full.jsonl"
        claims_full.write_text(
            "\n".join(
                [
                    json.dumps({"protein_id": "P1", "function_id": "1.1.1.1", "function_type": "EC", "raw_score": 0.80, "source": "clean"}),
                    json.dumps({"protein_id": "P2", "function_id": "2.7.1.1", "function_type": "EC", "raw_score": 0.20, "source": "clean"}),
                    json.dumps({"protein_id": "P3", "function_id": "3.1.1.1", "function_type": "EC", "raw_score": 0.06, "source": "clean"}),
                    json.dumps({"protein_id": "P4", "function_id": "4.1.1.1", "function_type": "EC", "raw_score": 0.01, "source": "clean"}),
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        full_tsv = self.root / "full.tsv"
        full_tsv.write_text(
            "protein_id\ttype\tfunction_id\tposterior_prob\n"
            "P1\tEC\tEC:1.1.1.1\t0.750000\n"
            "P2\tEC\tEC:2.7.1.1\t0.300000\n"
            "P3\tEC\tEC:3.1.1.1\t0.150000\n"
            "P4\tEC\tEC:4.1.1.1\t0.120000\n",
            encoding="utf-8",
        )
        filtered_tsv = self.root / "filtered.tsv"
        filtered_tsv.write_text(
            "protein_id\ttype\tfunction_id\tposterior_prob\n"
            "P1\tEC\tEC:1.1.1.1\t0.750000\n"
            "P2\tEC\tEC:2.7.1.1\t0.300010\n",
            encoding="utf-8",
        )
        out_json = self.root / "report.json"
        report = C.compare_full_vs_filtered(full_tsv, filtered_tsv, 0.05, claims_full, out_json)
        self.assertEqual(report["kept_pairs"], 3)
        self.assertEqual(report["compared_pairs"], 2)
        self.assertAlmostEqual(report["max_abs_diff"], 1e-5, places=9)
        self.assertEqual(report["count_diff_gt_1e_6"], 1)
        self.assertEqual(report["only_in_full"], 1)
        self.assertEqual(report["only_in_filtered"], 0)
        self.assertEqual(report["rows_only_in_one"], 1)
        self.assertEqual(report["called_full"], 1)
        self.assertEqual(report["called_filtered"], 1)
        self.assertEqual(report["called_status_changes"], 0)


if __name__ == "__main__":
    unittest.main()
