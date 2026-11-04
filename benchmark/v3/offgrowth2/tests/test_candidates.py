"""Unit tests for benchmark.v3.offgrowth2.candidates (spec §7)."""
from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from offgrowth2.candidates import (
    best_by_ec,
    best_by_protein,
    clean_sidecar_to_candidates,
    load_candidates,
)

TMP_BASE = Path("/home/leechuck/tmp")


class TestCandidates(unittest.TestCase):
    def test_clean_sidecar_conversion_and_helpers(self) -> None:
        TMP_BASE.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(TMP_BASE)) as tmpdir:
            sidecar_dir = Path(tmpdir) / "sidecar"
            sidecar_dir.mkdir(parents=True, exist_ok=True)
            clean_tsv = sidecar_dir / "GCF_00001.clean.tsv"
            clean_tsv.write_text(
                "protein_id\tterm\tscore\tannotation_type\tmetadata_json\n"
                "protB\tEC:1.1.1.1\t0.750000\tEC\t{}\n"
                "protA\tEC:2.7.1.1\t0.420000\tEC\t{}\n"
                "protB\tEC:1.1.1.1\t0.810000\tEC\t{}\n"
                "protZero\tEC:1.1.1.1\t0.000000\tEC\t{}\n"
                "protIncomp\tEC:1.1.1.-\t0.900000\tEC\t{}\n",
                encoding="utf-8",
            )
            out_tsv = Path(tmpdir) / "candidates.tsv"
            receipt = clean_sidecar_to_candidates(sidecar_dir, out_tsv)

            self.assertTrue(out_tsv.exists())
            self.assertTrue(Path(receipt["receipt_path"]).exists())
            self.assertEqual(receipt["written_rows"], 2)
            self.assertEqual(receipt["skipped_zero_score_rows"], 1)
            self.assertEqual(receipt["skipped_incomplete_ec_rows"], 1)
            self.assertEqual(receipt["duplicate_pairs_merged"], 1)

            cands = load_candidates(out_tsv)
            self.assertEqual(len(cands), 2)
            self.assertAlmostEqual(cands[("protA", "2.7.1.1")], 0.42)
            self.assertAlmostEqual(cands[("protB", "1.1.1.1")], 0.81)

            by_ec = best_by_ec(cands)
            self.assertEqual(by_ec["1.1.1.1"], [("protB", 0.81)])
            by_prot = best_by_protein(cands)
            self.assertEqual(by_prot["protA"], [("2.7.1.1", 0.42)])

    def test_score_bounds_enforced(self) -> None:
        TMP_BASE.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(TMP_BASE)) as tmpdir:
            bad_high = Path(tmpdir) / "bad_high.clean.tsv"
            bad_high.write_text(
                "protein_id\tterm\tscore\tannotation_type\n"
                "protA\tEC:1.1.1.1\t1.05\tEC\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                clean_sidecar_to_candidates(bad_high, Path(tmpdir) / "out1.tsv")

            bad_neg = Path(tmpdir) / "bad_neg.clean.tsv"
            bad_neg.write_text(
                "protein_id\tterm\tscore\tannotation_type\n"
                "protA\tEC:1.1.1.1\t-0.10\tEC\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                clean_sidecar_to_candidates(bad_neg, Path(tmpdir) / "out2.tsv")


if __name__ == "__main__":
    unittest.main()
