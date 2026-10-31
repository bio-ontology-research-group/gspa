#!/usr/bin/env python3
"""Tests for Swiss-Prot experimental EC labels pipeline (no network)."""

import csv
import json
from pathlib import Path
import tempfile
from typing import Any, Dict, List
import unittest

from offgrowth.runner import load_swissprot_experimental_ecs
from offgrowth.swissprot_labels import (
    EXPERIMENTAL_ECO_CODES,
    SwissprotLabelRecord,
    extract_entry_ecs,
    is_complete_four_level_ec,
    load_proteome_fasta,
    process_genome,
    sha256_file,
)


class TestSwissprotLabels(unittest.TestCase):
    """Unit tests for evidence filtering, EC validation, sequence mapping, and loader compatibility."""

    def test_evidence_filtering_drops_non_experimental(self) -> None:
        """Non-experimental evidence codes (ECO:0000255, ECO:0000250, ECO:0000305) must be dropped."""
        non_exp_entry: Dict[str, Any] = {
            "primaryAccession": "P12345",
            "proteinDescription": {
                "recommendedName": {
                    "ecNumbers": [
                        {
                            "value": "1.1.1.1",
                            "evidences": [
                                {"evidenceCode": "ECO:0000255", "source": "PROSITE-ProRule"},
                                {"evidenceCode": "ECO:0000250", "source": "UniProtKB"},
                                {"evidenceCode": "ECO:0000305", "source": "PubMed", "id": "12345"},
                            ],
                        }
                    ]
                }
            },
        }
        ecs = extract_entry_ecs(non_exp_entry)
        self.assertEqual(ecs, [], "Computational and curator-inferred ECs must be dropped")

        exp_entry: Dict[str, Any] = {
            "primaryAccession": "P67890",
            "proteinDescription": {
                "recommendedName": {
                    "ecNumbers": [
                        {
                            "value": "1.1.1.1",
                            "evidences": [
                                {"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "99999"},
                                {"evidenceCode": "ECO:0000255", "source": "PROSITE-ProRule"},
                            ],
                        },
                        {
                            "value": "2.7.1.1",
                            "evidences": [
                                {"evidenceCode": "ECO:0000314", "source": "PubMed", "id": "88888"},  # direct assay
                            ],
                        },
                    ]
                }
            },
        }
        exp_ecs = extract_entry_ecs(exp_entry)
        self.assertEqual(len(exp_ecs), 2)
        self.assertEqual(exp_ecs[0], ("1.1.1.1", "ECO:0000269"))
        self.assertEqual(exp_ecs[1], ("2.7.1.1", "ECO:0000314"))

    def test_incomplete_ec_dropped(self) -> None:
        """Incomplete or wildcard ECs must be dropped; only 4-level numeric ECs retained."""
        self.assertTrue(is_complete_four_level_ec("1.1.1.1"))
        self.assertTrue(is_complete_four_level_ec("5.4.99.2"))
        self.assertTrue(is_complete_four_level_ec("7.1.1.8"))

        self.assertFalse(is_complete_four_level_ec("4.1.2.-"))
        self.assertFalse(is_complete_four_level_ec("1.13.12.-"))
        self.assertFalse(is_complete_four_level_ec("2.7.4.-"))
        self.assertFalse(is_complete_four_level_ec("1.-.-.-"))
        self.assertFalse(is_complete_four_level_ec("2.7.1"))
        self.assertFalse(is_complete_four_level_ec("EC 1.1.1.1"))
        self.assertFalse(is_complete_four_level_ec("1.1.1.n"))

        incomplete_entry: Dict[str, Any] = {
            "primaryAccession": "P99999",
            "proteinDescription": {
                "recommendedName": {
                    "ecNumbers": [
                        {
                            "value": "4.1.2.-",
                            "evidences": [{"evidenceCode": "ECO:0000269", "source": "Ref"}],
                        },
                        {
                            "value": "1.13.12.-",
                            "evidences": [{"evidenceCode": "ECO:0000269", "source": "Ref"}],
                        },
                    ]
                }
            },
        }
        self.assertEqual(extract_entry_ecs(incomplete_entry), [])

    def test_exact_sequence_mapping(self) -> None:
        """Entries must map by exact 100% sequence identity; mismatches reported as unmatched."""
        with tempfile.TemporaryDirectory() as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)
            raw_dir = tmp_dir / "raw"
            raw_dir.mkdir()
            fasta_path = tmp_dir / "protein.faa"

            # Proteome contains P1 and P2 (identical sequence), and P3 (different sequence)
            fasta_content = (
                ">PROT_1 Alpha subunit\n"
                "MKVILVGDGAVGSSYAFALVNQGIA\n"
                ">PROT_2 Duplicated alpha subunit\n"
                "MKVILVGDGAVGSSYAFALVNQGIA\n"
                ">PROT_3 Beta subunit\n"
                "MAAALVGDGAVGSSYAFALVNQGIA\n"
            )
            fasta_path.write_text(fasta_content, encoding="utf-8")

            # Synthetic Swiss-Prot download:
            # - Entry 1 (SWISS_1): exact match to PROT_1 & PROT_2, EC 1.1.1.1 (ECO:0000269)
            # - Entry 2 (SWISS_2): 1 amino acid mismatch ("...GIZ"), EC 2.2.2.2 (ECO:0000269) -> unmatched
            # - Entry 3 (SWISS_3): match to PROT_3, but computational EC (ECO:0000255) -> dropped
            mock_entries = {
                "results": [
                    {
                        "primaryAccession": "SWISS_1",
                        "sequence": {"value": "MKVILVGDGAVGSSYAFALVNQGIA"},
                        "proteinDescription": {
                            "recommendedName": {
                                "ecNumbers": [
                                    {
                                        "value": "1.1.1.1",
                                        "evidences": [{"evidenceCode": "ECO:0000269"}],
                                    }
                                ]
                            }
                        },
                    },
                    {
                        "primaryAccession": "SWISS_2",
                        "sequence": {"value": "MKVILVGDGAVGSSYAFALVNQGIZ"},  # single mismatch at C-terminus
                        "proteinDescription": {
                            "recommendedName": {
                                "ecNumbers": [
                                    {
                                        "value": "2.2.2.2",
                                        "evidences": [{"evidenceCode": "ECO:0000269"}],
                                    }
                                ]
                            }
                        },
                    },
                    {
                        "primaryAccession": "SWISS_3",
                        "sequence": {"value": "MAAALVGDGAVGSSYAFALVNQGIA"},
                        "proteinDescription": {
                            "recommendedName": {
                                "ecNumbers": [
                                    {
                                        "value": "3.3.3.3",
                                        "evidences": [{"evidenceCode": "ECO:0000255"}],
                                    }
                                ]
                            }
                        },
                    },
                ]
            }
            taxid = 9999
            (raw_dir / f"{taxid}.json").write_text(json.dumps(mock_entries), encoding="utf-8")
            (raw_dir / f"{taxid}.sidecar.json").write_text(
                json.dumps({"taxid": taxid, "query_url": "mock://test", "sha256": "abc"}),
                encoding="utf-8",
            )

            labels, stats = process_genome(
                genome_id="MOCK_GENOME_1",
                organism_taxid=taxid,
                species_fallback_taxid=None,
                raw_dir=raw_dir,
                fasta_path=fasta_path,
                offline=True,
            )

            # SWISS_1 should map to both PROT_1 and PROT_2
            self.assertEqual(len(labels), 2)
            self.assertEqual({lbl.protein_id for lbl in labels}, {"PROT_1", "PROT_2"})
            self.assertTrue(all(lbl.function_id == "1.1.1.1" for lbl in labels))
            self.assertTrue(all(lbl.uniprot_accession == "SWISS_1" for lbl in labels))

            # SWISS_2 should be recorded as unmatched
            self.assertEqual(stats["matched_accessions_count"], 1)
            self.assertEqual(stats["unmatched_accessions_count"], 1)
            self.assertEqual(len(stats["unmatched_entries"]), 1)
            self.assertEqual(stats["unmatched_entries"][0]["uniprot_accession"], "SWISS_2")

    def test_output_loads_with_load_swissprot_experimental_ecs(self) -> None:
        """The output format must be loaded seamlessly by runner.load_swissprot_experimental_ecs."""
        with tempfile.TemporaryDirectory() as tmp_dir_str:
            tsv_path = Path(tmp_dir_str) / "labels.tsv"
            rows = [
                ["genome_id", "protein_id", "function_id", "uniprot_accession", "evidence", "tier"],
                ["GENOME_A", "PROT_1", "1.1.1.1", "P00001", "ECO:0000269", "experimental"],
                ["GENOME_A", "PROT_1", "1.1.1.2", "P00001", "ECO:0000269", "experimental"],  # bifunctional
                ["GENOME_A", "PROT_2", "2.7.1.1", "P00002", "ECO:0000269", "experimental"],
                ["GENOME_B", "PROT_9", "3.1.1.1", "P00009", "ECO:0000269", "experimental"],
                ["GENOME_A", "PROT_X", "4.1.1.1", "P0000X", "ECO:0000255", "computational"],  # non-exp tier
            ]
            with tsv_path.open("w", encoding="utf-8", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t")
                writer.writerows(rows)

            loaded_a = load_swissprot_experimental_ecs(tsv_path, "GENOME_A")
            self.assertIn("PROT_1", loaded_a)
            self.assertEqual(sorted(loaded_a["PROT_1"]), ["1.1.1.1", "1.1.1.2"])
            self.assertEqual(loaded_a["PROT_2"], ["2.7.1.1"])
            self.assertNotIn("PROT_X", loaded_a, "Non-experimental tier must be excluded by loader")

            loaded_b = load_swissprot_experimental_ecs(tsv_path, "GENOME_B")
            self.assertEqual(loaded_b, {"PROT_9": ["3.1.1.1"]})

            loaded_c = load_swissprot_experimental_ecs(tsv_path, "GENOME_C")
            self.assertEqual(loaded_c, {})

    def test_production_price22_dataset_integrity(self) -> None:
        """Validate production data/swissprot_experimental_ecs_price22.tsv against receipt."""
        data_dir = Path(__file__).resolve().parent.parent / "data"
        tsv_path = data_dir / "swissprot_experimental_ecs_price22.tsv"
        receipt_path = data_dir / "swissprot_experimental_ecs_price22.receipt.json"

        self.assertTrue(tsv_path.exists(), f"Production TSV missing at {tsv_path}")
        self.assertTrue(receipt_path.exists(), f"Production receipt missing at {receipt_path}")

        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        actual_sha = sha256_file(tsv_path)
        expected_sha = receipt["output_sha256"]["swissprot_experimental_ecs_price22_tsv"]
        self.assertEqual(actual_sha, expected_sha, "TSV sha256 must match receipt record")

        # Verify all rows have complete 4-level ECs and experimental tier
        with tsv_path.open("r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)

        self.assertEqual(len(rows), receipt["counts"]["total_labels"])
        self.assertGreater(len(rows), 0, "Labels table must not be empty")

        for r in rows:
            self.assertEqual(r["tier"], "experimental")
            self.assertIn(r["evidence"], EXPERIMENTAL_ECO_CODES)
            self.assertTrue(is_complete_four_level_ec(r["function_id"]))

        # Verify zero overlap with Price-148 truth
        self.assertEqual(receipt["price148_truth_overlap"]["total_overlapping_proteins"], 0)


if __name__ == "__main__":
    unittest.main()
