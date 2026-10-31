#!/usr/bin/env python3
"""Run controlled cohort & assembly perturbations on swissprot_label_rich_v1 (Items 6 and 8)."""
from __future__ import annotations

import copy
import csv
import json
import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))

import homology_baseline as H  # noqa: E402
import run_training_cohort_ablations as base  # noqa: E402
from v3_metrics import load_obo  # noqa: E402

CACHE = base.CACHE
RES = ROOT / "benchmark" / "v3" / "results" / "training-cohort"
PERT_DIR = CACHE / "perturbation_runs"
OUT_JSON = RES / "perturbation_controls.json"


def build_shuffled_layouts(out_dir: Path, seed: int = 42) -> Path:
    lay_dir = out_dir / "layouts_shuffled"
    lay_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    for src in sorted((CACHE / "layout").glob("*.layout.tsv")):
        lines = [ln.rstrip("\n") for ln in src.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
        header = lines[0]
        rows = [ln.split("\t") for ln in lines[1:]]
        pids = [r[0] for r in rows]
        rng.shuffle(pids)
        out_lines = [header]
        for pid, r in zip(pids, rows):
            out_lines.append("\t".join([pid] + r[1:]))
        (lay_dir / src.name).write_text("\n".join(out_lines) + "\n")
    return lay_dir


def build_fragmented_layouts(out_dir: Path, seed: int = 42) -> tuple[Path, Path]:
    """Split each genome into 5-gene draft contigs and drop 50% of contigs (unplaced proteins declared in metadata)."""
    lay_dir = out_dir / "layouts_fragmented"
    lay_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    meta = json.loads((CACHE / "cohort_inputs" / "cohort_metadata.json").read_text())
    frag_meta = copy.deepcopy(meta)

    for gid in sorted(frag_meta["genomes"].keys()):
        src = CACHE / "layout" / f"{gid}.layout.tsv"
        lines = [ln.rstrip("\n") for ln in src.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
        header = lines[0]
        rows = [ln.split("\t") for ln in lines[1:]]
        chunks = [rows[i : i + 5] for i in range(0, len(rows), 5)]
        kept_rows = []
        unplaced = set(frag_meta["genomes"][gid].get("unplacedProteins", []))
        contigs_meta = {}
        for c_idx, chunk in enumerate(chunks):
            cid = f"frag_{c_idx:04d}"
            if rng.random() < 0.5:
                for r in chunk:
                    unplaced.add(r[0])
            else:
                min_s = min(int(r[2]) for r in chunk)
                max_e = max(int(r[3]) for r in chunk)
                for r in chunk:
                    s_rel = int(r[2]) - min_s + 1
                    e_rel = int(r[3]) - min_s + 1
                    kept_rows.append([r[0], cid, str(s_rel), str(e_rel), r[4]])
                contigs_meta[cid] = {"length": max_e - min_s + 100, "closed": False}
        out_lines = [header] + ["\t".join(r) for r in kept_rows]
        (lay_dir / src.name).write_text("\n".join(out_lines) + "\n")
        frag_meta["genomes"][gid]["assemblyStatus"] = "DRAFT"
        frag_meta["genomes"][gid]["completenessFraction"] = 0.5
        frag_meta["genomes"][gid]["contigs"] = contigs_meta
        frag_meta["genomes"][gid]["unplacedProteins"] = sorted(unplaced)

    meta_path = out_dir / "cohort_metadata_fragmented.json"
    meta_path.write_text(json.dumps(frag_meta, indent=2) + "\n")
    return lay_dir, meta_path


def build_duplicated_strain_cohort(out_dir: Path, n_dups: int = 5) -> tuple[Path, Path, str, list[str]]:
    """Duplicate one training genome n_dups times in the same clade and strain cluster."""
    meta = json.loads((CACHE / "cohort_inputs" / "cohort_metadata.json").read_text())
    dup_meta = copy.deepcopy(meta)
    src_gid = sorted(meta["genomes"].keys())[0]
    dup_gids = []
    for i in range(1, n_dups + 1):
        dg = f"{src_gid}_dup{i}"
        dup_gids.append(dg)
        dup_meta["genomes"][dg] = copy.deepcopy(meta["genomes"][src_gid])

    meta_path = out_dir / "cohort_metadata_duplicated.json"
    meta_path.write_text(json.dumps(dup_meta, indent=2) + "\n")

    og_src = CACHE / "cohort_inputs" / "cohort_orthogroups.tsv"
    og_dst = out_dir / "cohort_orthogroups_duplicated.tsv"
    lines = og_src.read_text().splitlines()
    dup_lines = list(lines)
    src_rows = [ln for ln in lines[1:] if ln.startswith(src_gid + "\t")]
    for dg in dup_gids:
        for ln in src_rows:
            dup_lines.append(dg + ln[len(src_gid) :])
    og_dst.write_text("\n".join(dup_lines) + "\n")
    return meta_path, og_dst, src_gid, dup_gids


def run_perturbation_batch(
    name: str,
    genomes: list[str],
    layout_dir: Path,
    cohort_meta_path: Path,
    cohort_og_path: Path,
    learned_model_path: Path | None,
    claims_override: dict[str, Path] | None = None,
    layout_override: dict[str, Path] | None = None,
    rlc_suggester: bool = True,
    sources: tuple[Path, ...] = (),
) -> Path:
    run_dir = PERT_DIR / name
    out_dir = run_dir / "out"
    run_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    prov = json.loads((CACHE / "provenance.json").read_text())
    acc_to_tax = {
        r["assembly"]["used_accession"]: t
        for t, r in prov["organisms"].items()
        if r.get("assembly", {}).get("status") == "resolved"
    }
    sel = json.loads((ROOT / "benchmark/v3/results/label-rich-selection/selection_manifest.json").read_text())
    tax_to_sk = {o["taxid"]: o["superkingdom"].lower() for o in sel["organisms"]}

    manifest_path = run_dir / "batch_manifest.tsv"
    rows = [
        "genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\tgenome_layout_path\torthogroups_path\tcluster_consensus_path\tfba_medium\tfba_gram"
    ]
    for gid in genomes:
        base_gid = gid.split("_dup")[0]
        taxid = acc_to_tax[base_gid]
        sk = tax_to_sk[taxid]
        preset = "archaea" if "archaea" in sk else "bacteria"
        c_path = (claims_override or {}).get(gid, CACHE / "baseline" / "claims" / "homology" / f"{base_gid}.claims.jsonl")
        l_path = (layout_override or {}).get(gid, layout_dir / f"{base_gid}.layout.tsv")
        out_tsv = out_dir / f"{gid}.integrated.tsv"
        rows.append(
            f"{gid}\t{c_path}\t{out_tsv}\tNCBITaxon_{taxid}\t{preset}\t-\t-\t-\t-\t-\t-\t{l_path}\t-\t-\t-\t-"
        )
    manifest_path.write_text("\n".join(rows) + "\n")

    cmd = [
        "java", "-Xmx12g", "-jar", str(base.production_jar()), "integrate",
        "--batch-manifest", str(manifest_path),
        "--batch-threads", "4",
        "--go-owl", str(ROOT / "benchmark/whole_genome_cafa/data/reference/go.owl"),
        "--ec2go", str(ROOT / "benchmark/whole_genome_cafa/data/reference/ec2go.txt"),
        "--pathways", str(ROOT / "benchmark/whole_genome_cafa/data/reference/kegg_pathways.tsv"),
        "--modules", str(ROOT / "benchmark/whole_genome_cafa/data/reference/kegg_modules.tsv"),
        "--reasoner-cache", str(CACHE / "reasoner-cache"),
        "--reaction-graph", "/home/leechuck/Public/software/gapsmith/data/seed_reactions.tsv",
        "--diffusion-mets", "/home/leechuck/Public/software/gapsmith/data/diffusion_mets.tsv",
        "--reaction-ec-aliases", "/home/leechuck/Public/software/gapsmith/data/seed_Enzyme_Class_Reactions_Aliases_unique.tsv",
        "--cohort",
        "--cohort-metadata", str(cohort_meta_path),
        "--cohort-orthogroups", str(cohort_og_path),
        "--cohort-holdout", "CLADE",
    ]
    if rlc_suggester:
        cmd.append("--rlc-suggester")
    if learned_model_path is not None:
        cmd.extend(["--learned-model", str(learned_model_path)])

    print(f"[running perturbation] {name} ...", flush=True)
    base.run_batch(cmd, manifest_path, out_dir, sources=[Path(__file__), Path(base.__file__), *sources])
    return out_dir


def compare_test_cohort_features(dir_a: Path, dir_b: Path, test_genomes: list[str]) -> dict:
    report = dict(n_test_candidates_compared=0, max_abs_probability_diff=None,
                  unsupported_both=0, support_transitions=0, missing_before=0, missing_after=0)
    def candidates(folder, gid):
        rows = json.loads((folder / f'{gid}.integrated.tsv.learned.json').read_text())['candidates']
        result = {}
        for row in rows:
            key, value = row['candidate_key'], row['probability']
            if key in result or not isinstance(key, str) or not key:
                raise ValueError('invalid or duplicate learned candidate identity')
            if value is not None and (type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1):
                raise ValueError('invalid learned candidate probability')
            result[key] = value
        return result
    for gid in test_genomes:
        a, b = candidates(dir_a, gid), candidates(dir_b, gid)
        report['missing_before'] += len(b.keys() - a.keys())
        report['missing_after'] += len(a.keys() - b.keys())
        for key in a.keys() & b.keys():
            if a[key] is None and b[key] is None:
                report['unsupported_both'] += 1
            elif a[key] is None or b[key] is None:
                report['support_transitions'] += 1
            else:
                report['n_test_candidates_compared'] += 1
                report['max_abs_probability_diff'] = max(report['max_abs_probability_diff'] or 0, abs(a[key] - b[key]))
    return report


def main():
    # Keep the entry point, but require the verified prepare/run/calibrate/evaluate
    # interface. The historical no-argument main reused unverified baseline outputs
    # and chose Fmax thresholds on test outcomes.
    from run_cohort_controls import main as verified_main
    verified_main()


if __name__ == '__main__':
    main()
