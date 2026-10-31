#!/usr/bin/env python3
"""Unit tests for the Shaheen Empty Quarter (17,055) and benchmark (145) gapsmith pipeline.

Run with:
  TMPDIR=/home/leechuck/tmp /home/leechuck/tmp/gspa-v3-testenv/bin/python -m unittest discover \
    -s /mnt/data2/leechuck/gspa/benchmark/emptyquarter_mags/shaheen/tests -p 'test_*.py'
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import os
import stat
import sys
import tarfile
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

SHAHEEN_DIR = Path(__file__).resolve().parent.parent
if str(SHAHEEN_DIR) not in sys.path:
    sys.path.insert(0, str(SHAHEEN_DIR))

from build_manifest import (  # noqa: E402
    assign_eq_batches,
    assign_fixed_batches,
    build_benchmark_records,
    build_eq_records,
    write_manifests,
)
from check_models import check_models  # noqa: E402
from driver import format_commands, main as driver_main  # noqa: E402
from run_gapsmith_batch import (  # noqa: E402
    clean_medium_csv,
    parse_gnu_time_v,
    parse_sbml_counts,
    run_gapsmith_batch,
)
from stage_batch import (  # noqa: E402
    ensure_merged_data_dir,
    stage_batch,
    verify_pyrodigal_version,
)


def _get_base_tmp() -> str:
    cand = os.environ.get("TMPDIR", "/home/leechuck/tmp")
    os.makedirs(cand, exist_ok=True)
    return cand


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


SAMPLE_SBML_MODEL = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2"
      level="3" version="1">
  <model id="test_model">
    <fbc:listOfGeneProducts>
      <fbc:geneProduct fbc:id="G_prot1" fbc:label="prot1"/>
      <fbc:geneProduct fbc:id="G_prot2" fbc:label="prot2"/>
    </fbc:listOfGeneProducts>
    <listOfSpecies>
      <species id="M_cpd00001_c0" compartment="c0"/>
      <species id="M_cpd00002_c0" compartment="c0"/>
      <species id="M_cpd00003_c0" compartment="c0"/>
    </listOfSpecies>
    <listOfReactions>
      <reaction id="R_rxn00001_c0" reversible="true"/>
      <reaction id="R_rxn00002_c0" reversible="false"/>
      <reaction id="R_rxn00003_c0" reversible="false"/>
      <reaction id="R_EX_cpd00001_e0" reversible="true"/>
    </listOfReactions>
  </model>
</sbml>
"""

SAMPLE_SBML_DRAFT = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2"
      level="3" version="1">
  <model id="test_draft">
    <fbc:listOfGeneProducts>
      <fbc:geneProduct fbc:id="G_prot1" fbc:label="prot1"/>
      <fbc:geneProduct fbc:id="G_prot2" fbc:label="prot2"/>
    </fbc:listOfGeneProducts>
    <listOfSpecies>
      <species id="M_cpd00001_c0" compartment="c0"/>
      <species id="M_cpd00002_c0" compartment="c0"/>
    </listOfSpecies>
    <listOfReactions>
      <reaction id="R_rxn00001_c0" reversible="true"/>
      <reaction id="R_rxn00002_c0" reversible="false"/>
    </listOfReactions>
  </model>
</sbml>
"""


class TestShaheenEqModels(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_dir = Path(tempfile.mkdtemp(prefix="test_shaheen_eq_", dir=_get_base_tmp()))

    def tearDown(self) -> None:
        import shutil

        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def test_build_eq_manifest_full_17055(self) -> None:
        """Verify manifest building across all 17,055 rows of eq_mags_per_bin.tsv."""
        eq_tsv = SHAHEEN_DIR.parent / "results" / "eq_mags_per_bin.tsv"
        self.assertTrue(eq_tsv.is_file(), "Missing eq_mags_per_bin.tsv")

        # Synthesize sites.tsv covering every soil_bin sample in eq_mags_per_bin.tsv
        soil_samples = set()
        with open(eq_tsv, "r", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["dataset"] == "soil_bin":
                    soil_samples.add(row["sample"])

        sites_tsv = self.tmp_dir / "sites.tsv"
        with open(sites_tsv, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample_id", "trip", "site", "feature", "has_wgs", "has_amplicon", "tier1_bytes"])
            for s in sorted(soil_samples):
                w.writerow([s, "trip1", "site10", "D", "1", "0", "1000"])

        records = build_eq_records(
            eq_per_bin_tsv=eq_tsv,
            sites_tsv=sites_tsv,
            sites_tar_dir=Path("/project/k10018/empty-quarter/data/sites-tar"),
            rakpr_dir=Path("/project/k10018/empty-quarter/data/rakpr"),
        )
        self.assertEqual(len(records), 17055)

        n_bakta = sum(1 for r in records if r["protein_source"] == "bakta")
        n_pyro = sum(1 for r in records if r["protein_source"] == "pyrodigal")
        self.assertEqual(n_bakta, 16217)
        self.assertEqual(n_pyro, 838)

        # Ensure pyrodigal records never have a bakta_prefix set
        for r in records:
            if r["protein_source"] == "pyrodigal":
                self.assertEqual(r["bakta_prefix"], "")

        batches = assign_eq_batches(records, batch_size=500)
        self.assertEqual(len(batches), 37)
        self.assertEqual(sum(len(b) for b in batches), 17055)

        # Verify each soil sample tar belongs to exactly 1 batch (no tar split across batches)
        sample_to_batches = {}
        for b in batches:
            bid = b[0]["batch_id"]
            for r in b:
                if r["dataset"] == "soil_bin":
                    sample_to_batches.setdefault(r["sample"], set()).add(bid)
        for sample_id, bids in sample_to_batches.items():
            self.assertEqual(
                len(bids),
                1,
                "Soil sample {0} was split across batches {1}".format(sample_id, bids),
            )

        work_root = self.tmp_dir / "work_eq"
        all_path, idx_path = write_manifests(batches, work_root)
        self.assertTrue(all_path.is_file())
        self.assertTrue(idx_path.is_file())

    def test_build_benchmark_manifest_145(self) -> None:
        """Verify benchmark manifest builder for 22 Price + 18 cohort + 105 peers (145 total)."""
        peers_tsv = SHAHEEN_DIR.parent.parent / "v3" / "results" / "peer-selection-v1" / "peers.tsv"
        self.assertTrue(peers_tsv.is_file())

        with open(peers_tsv, "r", encoding="utf-8") as pfh:
            peer_ids = sorted(
                {row["peer_accession"].strip() for row in csv.DictReader(pfh, delimiter="\t")}
            )
        self.assertEqual(len(peer_ids), 105)

        price_ids = ["GCA_{0:09d}.1".format(i) for i in range(1, 23)]
        cohort_ids = ["GCA_9{0:08d}.1".format(i) for i in range(1, 19)]

        price_dir = self.tmp_dir / "price"
        (price_dir / "layout").mkdir(parents=True, exist_ok=True)
        targets_json = price_dir / "targets.json"
        targets_json.write_text(
            json.dumps({"aliases": {gid: {} for gid in price_ids}}),
            encoding="utf-8",
        )
        report_json = price_dir / "layout" / "report.json"
        report_json.write_text(
            json.dumps(
                {
                    "assemblies": {
                        gid: {
                            "proteins": 3000 + idx,
                            "input_hashes": {"protein.faa": "deadbeef" * 8},
                        }
                        for idx, gid in enumerate(price_ids)
                    }
                }
            ),
            encoding="utf-8",
        )

        wp4_dir = self.tmp_dir / "wp4-search-v1"
        for set_name, gids in (("cohort", cohort_ids), ("peers", peer_ids)):
            arch_gids = gids[: len(gids) // 2]
            bact_gids = gids[len(gids) // 2 :]
            for sk, sub_gids, dom in (
                ("archaea", arch_gids, "Archaea"),
                ("bacteria", bact_gids, "Bacteria"),
            ):
                sdir = wp4_dir / set_name / sk
                (sdir / "proteomes").mkdir(parents=True, exist_ok=True)
                m = {
                    "superkingdom": dom,
                    "genomes": [
                        {
                            "genome_id": g,
                            "proteins": 2000,
                            "proteome": str(sdir / "proteomes" / (g + ".faa")),
                            "sha256": "cafebabe" * 8,
                            "taxonomy": dom,
                        }
                        for g in sub_gids
                    ],
                }
                (sdir / "manifest.json").write_text(json.dumps(m), encoding="utf-8")

        records = build_benchmark_records(
            price_targets_json=targets_json,
            price_genomes_dir=self.tmp_dir / "price_genomes",
            wp4_search_dir=wp4_dir,
            peers_tsv=peers_tsv,
            price_layout_report=report_json,
        )
        self.assertEqual(len(records), 145)
        self.assertEqual(sum(1 for r in records if r["dataset"] == "price_matrix"), 22)
        self.assertEqual(sum(1 for r in records if r["dataset"] == "recurrence_cohort"), 18)
        self.assertEqual(sum(1 for r in records if r["dataset"] == "recurrence_peer"), 105)

        batches = assign_fixed_batches(records, batch_size=500)
        self.assertEqual(len(batches), 1)
        self.assertEqual(len(batches[0]), 145)

    def test_ensure_merged_data_dir_and_pyrodigal_version(self) -> None:
        """Verify symlinked merged reference directory and strict pyrodigal 3.7.1 version check."""
        bact_ref = self.tmp_dir / "bact_ref"
        arch_ref = self.tmp_dir / "arch_ref"
        (bact_ref / "biomass").mkdir(parents=True)
        (bact_ref / "media").mkdir(parents=True)
        (bact_ref / "seq" / "Bacteria" / "rxn").mkdir(parents=True)
        (bact_ref / "seq" / "Archaea").mkdir(parents=True)
        (bact_ref / "seq" / "tcdb.fasta").write_text(">tcdb\nMAAA\n", encoding="utf-8")
        (bact_ref / "seed_reactions_corrected.tsv").write_text("id\tname\n", encoding="utf-8")

        (arch_ref / "seq" / "Archaea" / "rxn").mkdir(parents=True)
        (arch_ref / "seq" / "Archaea" / "rxn" / "rxn00001.fasta").write_text(
            ">a1\nMKKK\n", encoding="utf-8"
        )

        merged_dir = self.tmp_dir / "gapsmith_data_merged"
        ensure_merged_data_dir(bact_ref, arch_ref, merged_dir)

        self.assertTrue((merged_dir / "biomass").is_symlink())
        self.assertTrue((merged_dir / "seq" / "Bacteria").is_symlink())
        self.assertTrue((merged_dir / "seq" / "Archaea").is_symlink())
        self.assertEqual(
            (merged_dir / "seq" / "Archaea").resolve(),
            (arch_ref / "seq" / "Archaea").resolve(),
        )
        self.assertTrue((merged_dir / "seq" / "Archaea" / "rxn" / "rxn00001.fasta").is_file())

        # Test pyrodigal version check
        good_pyro = self.tmp_dir / "pyrodigal_good"
        _write_executable(
            good_pyro,
            "#!/usr/bin/env bash\necho 'pyrodigal v3.7.1'\n",
        )
        self.assertEqual(verify_pyrodigal_version(str(good_pyro), "3.7.1"), "3.7.1")

        bad_pyro = self.tmp_dir / "pyrodigal_bad"
        _write_executable(
            bad_pyro,
            "#!/usr/bin/env bash\necho 'pyrodigal v3.6.3'\n",
        )
        with self.assertRaises(RuntimeError):
            verify_pyrodigal_version(str(bad_pyro), "3.7.1")

    def test_stage_batch_rules_and_end_to_end_pipeline(self) -> None:
        """Test Stage 1 + Stage 2 + Checker end-to-end with mock pyrodigal and gapsmith."""
        # Create mock pyrodigal 3.7.1 that emits 2 proteins from any input FASTA
        mock_pyrodigal = self.tmp_dir / "mock_pyrodigal"
        _write_executable(
            mock_pyrodigal,
            """#!/usr/bin/env bash
set -euo pipefail
if [ "${1:-}" = "--version" ]; then
    echo "pyrodigal v3.7.1"
    exit 0
fi
FAA=""
GFF=""
FFN=""
while [ $# -gt 0 ]; do
    case "$1" in
        -a) FAA="$2"; shift 2 ;;
        -o) GFF="$2"; shift 2 ;;
        -d) FFN="$2"; shift 2 ;;
        -i|-f|-p|-j) shift 2 ;;
        *) echo "pyrodigal: error: unrecognized arguments: $1" >&2; exit 2 ;;   # like pyrodigal 3.7.1 (no -q)
    esac
done
printf ">p1\\nMAAAA\\n>p2\\nMTTTT\\n" > "$FAA"
printf "##gff-version 3\\n" > "$GFF"
if [ -n "$FFN" ]; then printf ">p1\\nATG\\n" > "$FFN"; fi
""",
        )

        # Create a site tar with 4 genomes:
        # 1. S1__Bin_bakta_ok (bakta, Bacteria, expected=3, actual=3) -> staged_ok -> computed_gapsmith
        # 2. S1__Bin_pyro_arch (pyrodigal, Archaea, expected=2, actual=2; also has a 9-protein decoy
        #    bakta .faa in the tar to prove pyrodigal never uses bakta!) -> staged_ok -> computed_gapsmith_refill
        # 3. S1__Bin_mismatch (bakta, Bacteria, expected=5, actual=3) -> protein_count_mismatch (not built)
        # 4. S1__Bin_missing (bakta, Bacteria, expected=4, missing bakta .faa even though bin .fa.gz exists)
        #    -> missing_input (no silent fallback to pyrodigal!)
        tar_build_dir = self.tmp_dir / "tar_content"
        sample_rel = Path("trip1/site10/S1")
        bins_dir = tar_build_dir / sample_rel / "wgs/binning/single_sample/output_bins"
        bakta_root = tar_build_dir / sample_rel / "wgs/annotation/single_sample/bakta"
        bins_dir.mkdir(parents=True)

        for bname in ("Bin_bakta_ok", "Bin_pyro_arch", "Bin_mismatch", "Bin_missing"):
            with gzip.open(bins_dir / (bname + ".fa.gz"), "wt", encoding="utf-8") as gz:
                gz.write(">contig1\nATGCATGCATGC\n")

        # Bin_bakta_ok: 3 proteins
        b1 = bakta_root / "Bin_bakta_ok"
        b1.mkdir(parents=True)
        (b1 / "Bin_bakta_ok.faa").write_text(">b1\nMAA\n>b2\nMCC\n>b3\nMGG\n", encoding="utf-8")
        with gzip.open(b1 / "Bin_bakta_ok.gff3.gz", "wt", encoding="utf-8") as gz:
            gz.write("##gff-version 3\n")

        # Bin_pyro_arch: decoy bakta with 9 proteins (must NOT be used!)
        b2 = bakta_root / "Bin_pyro_arch"
        b2.mkdir(parents=True)
        (b2 / "Bin_pyro_arch.faa").write_text(
            "".join(">d{0}\nMAA\n".format(i) for i in range(9)), encoding="utf-8"
        )

        # Bin_mismatch: 3 proteins in .faa, but expected_n_proteins = 5
        b3 = bakta_root / "Bin_mismatch"
        b3.mkdir(parents=True)
        (b3 / "Bin_mismatch.faa").write_text(">m1\nMAA\n>m2\nMCC\n>m3\nMGG\n", encoding="utf-8")

        site_tar_path = self.tmp_dir / "S1.tar"
        with tarfile.open(site_tar_path, "w") as tf:
            tf.add(tar_build_dir / "trip1", arcname="trip1")

        # Build reference table and batch manifest
        ref_tsv = self.tmp_dir / "eq_mags_per_bin_test.tsv"
        ref_rows = [
            {
                "mag_id": "S1__Bin_bakta_ok",
                "dataset": "soil_bin",
                "sample": "S1",
                "bin": "Bin_bakta_ok",
                "domain": "Bacteria",
                "protein_source": "bakta",
                "n_proteins": "3",
                "gapseq_source": "computed_gapsmith",
            },
            {
                "mag_id": "S1__Bin_pyro_arch",
                "dataset": "soil_bin",
                "sample": "S1",
                "bin": "Bin_pyro_arch",
                "domain": "Archaea",
                "protein_source": "pyrodigal",
                "n_proteins": "2",
                "gapseq_source": "computed_gapsmith_refill",
            },
            {
                "mag_id": "S1__Bin_mismatch",
                "dataset": "soil_bin",
                "sample": "S1",
                "bin": "Bin_mismatch",
                "domain": "Bacteria",
                "protein_source": "bakta",
                "n_proteins": "5",
                "gapseq_source": "computed_gapsmith",
            },
            {
                "mag_id": "S1__Bin_missing",
                "dataset": "soil_bin",
                "sample": "S1",
                "bin": "Bin_missing",
                "domain": "Bacteria",
                "protein_source": "bakta",
                "n_proteins": "4",
                "gapseq_source": "computed_gapsmith",
            },
        ]
        with open(ref_tsv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(ref_rows[0].keys()), delimiter="\t")
            w.writeheader()
            for r in ref_rows:
                w.writerow(r)

        sites_tsv = self.tmp_dir / "sites_test.tsv"
        sites_tsv.write_text(
            "sample_id\ttrip\tsite\nS1\ttrip1\tsite10\n", encoding="utf-8"
        )
        sites_tar_root = self.tmp_dir / "sites-tar"
        (sites_tar_root / "trip1/site10").mkdir(parents=True)
        os.symlink(site_tar_path, sites_tar_root / "trip1/site10/S1.tar")

        work_root = self.tmp_dir / "work"
        eq_recs = build_eq_records(
            eq_per_bin_tsv=ref_tsv,
            sites_tsv=sites_tsv,
            sites_tar_dir=sites_tar_root,
            rakpr_dir=self.tmp_dir / "rakpr",
        )
        batches = assign_eq_batches(eq_recs, batch_size=500)
        write_manifests(batches, work_root)

        batch_dir = work_root / "batches" / "batch_000"
        statuses = stage_batch(
            batch_dir=batch_dir,
            pyrodigal_bin=str(mock_pyrodigal),
            expected_pyrodigal_version="3.7.1",
            workers=2,
            scratch_tmp=self.tmp_dir,
        )
        s_map = {r["mag_id"]: r for r in statuses}
        self.assertEqual(s_map["S1__Bin_bakta_ok"]["status"], "staged_ok")
        self.assertEqual(s_map["S1__Bin_bakta_ok"]["ready_for_gapsmith"], "1")
        self.assertEqual(s_map["S1__Bin_pyro_arch"]["status"], "staged_ok")
        self.assertEqual(s_map["S1__Bin_pyro_arch"]["actual_n_proteins"], "2")
        self.assertEqual(s_map["S1__Bin_pyro_arch"]["ready_for_gapsmith"], "1")
        self.assertEqual(s_map["S1__Bin_mismatch"]["status"], "protein_count_mismatch")
        self.assertEqual(s_map["S1__Bin_mismatch"]["ready_for_gapsmith"], "0")
        self.assertEqual(s_map["S1__Bin_missing"]["status"], "missing_input")
        self.assertEqual(s_map["S1__Bin_missing"]["ready_for_gapsmith"], "0")

        # Verify inputs.tar ONLY contains the 2 verified genomes
        with tarfile.open(batch_dir / "inputs.tar", "r:*") as tf:
            top_dirs = {m.name.split("/", 1)[0] for m in tf.getmembers()}
        self.assertEqual(top_dirs, {"S1__Bin_bakta_ok", "S1__Bin_pyro_arch"})

        # Create mock gapsmith-linux (0.1.0) that:
        # - For S1__Bin_bakta_ok: succeeds on `doall` and writes SBML model + draft
        # - For S1__Bin_pyro_arch: fails on `doall` after writing draft.gmod.cbor + medium.csv
        #   with an embedded comma, then succeeds on `fill --media <cleaned.csv>`
        call_counter = self.tmp_dir / "gapsmith_calls.log"
        mock_gapsmith = self.tmp_dir / "mock_gapsmith"
        _write_executable(
            mock_gapsmith,
            """#!/usr/bin/env python3
import os, sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args or "-V" in args:
    print("gapsmith 0.1.0")
    sys.exit(0)

counter_file = Path({counter!r})
with open(counter_file, "a", encoding="utf-8") as cf:
    cf.write(" ".join(args) + "\\n")

if "doall" in args:
    out_dir = Path(args[args.index("-f") + 1])
    tax = args[args.index("-t") + 1]
    faa = Path(args[-1])
    mag = faa.stem
    print("running doall tax=" + tax + " mag=" + mag)
    (out_dir / (mag + "-draft.xml")).write_text({draft_xml!r}, encoding="utf-8")
    (out_dir / (mag + "-draft.gmod.cbor")).write_bytes(b"cbor_draft")
    (out_dir / (mag + "-Reactions.tbl")).write_text("rxn\\tname\\nrxn00001\\tr1\\nrxn00002\\tr2\\n", encoding="utf-8")
    if mag == "S1__Bin_pyro_arch":
        (out_dir / (mag + "-medium.csv")).write_text(
            "compounds,name,maxFlux\\ncpd001,2,3-dihydroxybenzoate,100\\n", encoding="utf-8"
        )
        sys.exit(1)
    else:
        (out_dir / (mag + "-medium.csv")).write_text("compounds,name,maxFlux\\ncpd001,glucose,100\\n", encoding="utf-8")
        (out_dir / (mag + "-filled.xml")).write_text({model_xml!r}, encoding="utf-8")
        sys.exit(0)
elif "fill" in args:
    out_dir = Path(args[args.index("-o") + 1])
    med = Path(args[args.index("--media") + 1])
    med_text = med.read_text(encoding="utf-8")
    for line in med_text.strip().splitlines()[1:]:
        if len(line.split(",")) != 3:
            sys.exit(2)
    cbor = Path(args[-1])
    mag = cbor.name.replace("-draft.gmod.cbor", "")
    (out_dir / (mag + "-filled.xml")).write_text({model_xml!r}, encoding="utf-8")
    sys.exit(0)
""".format(
                counter=str(call_counter),
                draft_xml=SAMPLE_SBML_DRAFT,
                model_xml=SAMPLE_SBML_MODEL,
            ),
        )

        data_dir = self.tmp_dir / "data_merged"
        data_dir.mkdir(parents=True)

        summary1 = run_gapsmith_batch(
            batch_dir=batch_dir,
            gapsmith_bin=str(mock_gapsmith),
            data_dir=data_dir,
            diamond_bin="diamond",
            expected_gapsmith_version="0.1.0",
            workers=2,
            threads_per_genome=4,
            scratch_tmp=self.tmp_dir,
        )
        sum_map = {r["mag_id"]: r for r in summary1}
        self.assertEqual(sum_map["S1__Bin_bakta_ok"]["gapseq_source_shaheen"], "computed_gapsmith")
        self.assertEqual(sum_map["S1__Bin_bakta_ok"]["n_reactions_model"], "4")
        self.assertEqual(sum_map["S1__Bin_bakta_ok"]["n_genes_model"], "2")

        self.assertEqual(
            sum_map["S1__Bin_pyro_arch"]["gapseq_source_shaheen"], "computed_gapsmith_refill"
        )
        self.assertEqual(sum_map["S1__Bin_pyro_arch"]["n_reactions_model"], "4")
        self.assertEqual(sum_map["S1__Bin_mismatch"]["gapseq_source_shaheen"], "protein_count_mismatch")
        self.assertEqual(sum_map["S1__Bin_missing"]["gapseq_source_shaheen"], "missing_input")

        calls_before_resume = call_counter.read_text(encoding="utf-8").splitlines()
        # 2 doall calls + 1 fill refill call = 3 invocations
        self.assertEqual(len(calls_before_resume), 3)
        self.assertTrue(any("-t Archaea" in c for c in calls_before_resume))
        self.assertTrue(any("-t Bacteria" in c for c in calls_before_resume))

        # Test resumability: running run_gapsmith_batch a second time skips both DONE genomes
        summary2 = run_gapsmith_batch(
            batch_dir=batch_dir,
            gapsmith_bin=str(mock_gapsmith),
            data_dir=data_dir,
            diamond_bin="diamond",
            expected_gapsmith_version="0.1.0",
            workers=2,
            threads_per_genome=4,
            scratch_tmp=self.tmp_dir,
        )
        self.assertEqual(len(summary2), 4)
        calls_after_resume = call_counter.read_text(encoding="utf-8").splitlines()
        self.assertEqual(
            len(calls_after_resume),
            3,
            "Resumed batch should not re-invoke gapsmith doall/fill for DONE genomes",
        )

        # Test check_models.py
        report = check_models(
            ref_tsv=ref_tsv,
            work_root=work_root,
            inspect_tars=True,
        )
        self.assertEqual(report["n_reference_genomes"], 4)
        self.assertTrue(report["all_reference_present_in_manifests"])
        self.assertTrue(report["no_mismatched_genome_built"])
        self.assertEqual(report["n_built_with_protein_mismatch"], 0)
        self.assertEqual(report["status_counts"]["computed_gapsmith"], 1)
        self.assertEqual(report["status_counts"]["computed_gapsmith_refill"], 1)
        self.assertEqual(report["status_counts"]["protein_count_mismatch"], 1)
        self.assertEqual(report["status_counts"]["missing_input"], 1)
        self.assertEqual(report["model_metrics"]["overall"]["reactions_filled"]["count"], 2)
        self.assertEqual(report["model_metrics"]["overall"]["reactions_filled"]["mean"], 4.0)
        self.assertEqual(report["model_metrics"]["overall"]["genes_filled"]["mean"], 2.0)

    def test_clean_medium_csv_and_sbml_and_time_parsers(self) -> None:
        """Unit test clean_medium_csv, parse_sbml_counts, and parse_gnu_time_v."""
        raw_med = self.tmp_dir / "med.csv"
        clean_med = self.tmp_dir / "med.clean.csv"
        raw_med.write_text(
            "compounds,name,maxFlux\ncpd00027,D-Glucose,100\ncpd00123,2,3-dihydroxybenzoate,50\n",
            encoding="utf-8",
        )
        n_fixed = clean_medium_csv(raw_med, clean_med)
        self.assertEqual(n_fixed, 1)
        lines = clean_med.read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(lines[1], "cpd00027,D-Glucose,100")
        self.assertEqual(lines[2], "cpd00123,23-dihydroxybenzoate,50")

        xml_file = self.tmp_dir / "m.xml"
        xml_file.write_text(SAMPLE_SBML_MODEL, encoding="utf-8")
        self.assertEqual(parse_sbml_counts(xml_file), (4, 2, 3))

        time_file = self.tmp_dir / "time_v.txt"
        time_file.write_text(
            "\tElapsed (wall clock) time (h:mm:ss or m:ss): 2:34.50\n"
            "\tMaximum resident set size (kbytes): 1456000\n",
            encoding="utf-8",
        )
        wall, rss = parse_gnu_time_v(time_file)
        self.assertAlmostEqual(wall, 154.50, places=2)
        self.assertEqual(rss, 1456000)

    def test_driver_dry_run_prints_sbatch_commands(self) -> None:
        """Verify driver.py prints sbatch commands with required Shaheen flags without submitting."""
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = driver_main(["--mode", "eq", "--work-root", str(self.tmp_dir / "drv")])
        self.assertEqual(rc, 0)
        out = buf.getvalue()
        self.assertIn("--account=k10018", out)
        self.assertIn("--partition=workq", out)
        self.assertIn("--exclude=nid00022", out)
        self.assertIn("--time=24:00:00", out)
        self.assertIn("--cpus-per-task=192", out)
        self.assertIn("--exclusive", out)
        self.assertIn("--array=0-36", out)
        self.assertIn("Dry-run only", out)


if __name__ == "__main__":
    unittest.main()
