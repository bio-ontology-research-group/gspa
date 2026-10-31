#!/usr/bin/env python3
"""Run controlled cohort & assembly perturbations on swissprot_label_rich_v1 (Items 6 and 8)."""
from __future__ import annotations

import copy
import csv
import json
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
    max_diff = 0.0
    n_compared = 0
    for gid in test_genomes:
        la = json.loads((dir_a / f"{gid}.integrated.tsv.learned.json").read_text())
        lb = json.loads((dir_b / f"{gid}.integrated.tsv.learned.json").read_text())
        map_a = {c["candidate_key"]: c["probability"] for c in la.get("candidates", [])}
        map_b = {c["candidate_key"]: c["probability"] for c in lb.get("candidates", [])}
        for k, pa in map_a.items():
            if k in map_b:
                n_compared += 1
                max_diff = max(max_diff, abs(pa - map_b[k]))
    return {"n_test_candidates_compared": n_compared, "max_abs_probability_diff": max_diff}


def main():
    PERT_DIR.mkdir(parents=True, exist_ok=True)
    pid_to_seq = {}
    with open(CACHE / "prepared/proteins.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pid_to_seq[(r["genome_id"], r["genome_protein_id"])] = r["protein_id"]

    split_of = {}
    with open(RES / "splits.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            split_of[r["genome_id"]] = r["split"]

    go = load_obo(str(CACHE / "go-basic.obo"))
    truth = H.read_truth(RES / "labels.tsv", CACHE / "prepared/proteins.tsv")
    base_genomes = sorted(split_of.keys())
    test_genomes = sorted(g for g, s in split_of.items() if s == "test")
    logistic_model = RES / "learned_models" / "model_logistic.json"

    shuf_lay = build_shuffled_layouts(PERT_DIR, seed=42)
    frag_lay, frag_meta = build_fragmented_layouts(PERT_DIR, seed=42)
    dup_meta, dup_og, src_gid, dup_gids = build_duplicated_strain_cohort(PERT_DIR, n_dups=5)

    unperturbed_dir = CACHE / "ablation_runs" / "homology_learned_logistic"
    shuf_dir = run_perturbation_batch(
        "shuffled_neighborhoods_logistic",
        base_genomes,
        shuf_lay,
        CACHE / "cohort_inputs" / "cohort_metadata.json",
        CACHE / "cohort_inputs" / "cohort_orthogroups.tsv",
        logistic_model,
    )
    frag_dir = run_perturbation_batch(
        "fragmented_assembly_50pct_logistic",
        base_genomes,
        frag_lay,
        frag_meta,
        CACHE / "cohort_inputs" / "cohort_orthogroups.tsv",
        logistic_model,
    )
    dup_dir = run_perturbation_batch(
        "duplicated_strain_5x_logistic",
        base_genomes + dup_gids,
        CACHE / "layout",
        dup_meta,
        dup_og,
        logistic_model,
    )

    def eval_dir(d: Path) -> dict:
        preds = base.load_integrated_predictions(d, pid_to_seq)
        return H.evaluate(preds, truth, split_of, go=go, n_boot=200)["test"]

    res = {
        "schema": "gspa.training_cohort_perturbations/1",
        "baseline_unperturbed": eval_dir(unperturbed_dir),
        "shuffled_neighborhoods": eval_dir(shuf_dir),
        "fragmented_assembly_50pct": eval_dir(frag_dir),
        "duplicated_strain_5x": eval_dir(dup_dir),
        "duplicated_strain_invariance_check": compare_test_cohort_features(unperturbed_dir, dup_dir, test_genomes),
    }
    OUT_JSON.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    print("Wrote", OUT_JSON)
    for k in ("baseline_unperturbed", "shuffled_neighborhoods", "fragmented_assembly_50pct", "duplicated_strain_5x"):
        m = res[k]
        print(f"  {k:28s}: test GO Fmax={m['go_fmax']['estimate']:.4f}, EC top1={m['ec_top1_accuracy']['estimate']:.4f}")
    print("  duplicated_strain_invariance_check:", res["duplicated_strain_invariance_check"])


if __name__ == "__main__":
    main()
