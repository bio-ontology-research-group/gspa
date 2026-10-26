#!/usr/bin/env python3
"""Turn the acquired cohort into audited split/protein/family inputs for ``prepare_training_data.py``.

Reads the provenance manifest (``acquire_training_cohort.py``), the selection, and the label matches
(``match_swissprot_labels.py``). Writes, into ``--out-dir``:

``splits.tsv``            genome_id split group source species genus family - only organisms whose assembly
                          resolved AND that the Price audit did not exclude. ``group`` is the family unit the
                          selection split on; ``genus``/``family`` come from NCBI Taxonomy, so the audit checks
                          disjointness on both the selection's and NCBI's family.
``proteins.tsv``          every acquired protein. The protein id is a cross-genome identity (sequence hash);
                          labelled proteins are ``target``, all others ``context``.
``protein_clusters.tsv``  MMseqs2 sequence-similarity clusters (default 30 % identity, 80 % coverage, the same
                          setting as the Price families) over the unique acquired sequences.
``example_clusters.tsv``  the same clusters restricted to labelled proteins (what ``build-table`` masks on).
``split_audit.json``      the ``audit-splits`` report (exit status of the audit is propagated).

A cluster is a sequence-similarity family, not an ortholog group, and says nothing about functional
independence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prepare_price_layout as PL  # noqa: E402

SOURCE = "swissprot_label_rich_v1"
SELECTION = HERE / "results" / "label-rich-selection" / "selection_manifest.json"


class CohortError(ValueError):
    pass


def seq_id(sequence: str) -> str:
    return "seq_" + hashlib.sha256(sequence.encode()).hexdigest()[:20]


def eligible_genomes(provenance: dict, selection: dict) -> list:
    """Organisms that may enter the cohort: resolved assembly, not removed by the Price audit."""
    by_tax = {o["taxid"]: o for o in selection["organisms"]}
    out = []
    for taxid, rec in sorted(provenance["organisms"].items()):
        asm = rec.get("assembly", {})
        audit = rec.get("price_audit")
        if audit is None:
            raise CohortError(f"{taxid}: the Price audit has not been run; refusing to build a cohort without it")
        if asm.get("status") != "resolved" or audit["excluded"]:
            continue
        sel, ncbi = by_tax[taxid], audit["ncbi_classification"]
        out.append({"taxid": taxid, "genome_id": asm["used_accession"], "directory": asm["directory"],
                    "split": sel["proposed_split"], "group": sel["family"] or f"genus:{sel['genus']}",
                    "species": sel["species"], "genus": ncbi["genus"], "family": ncbi["family"]})
    return out


def write_splits(genomes, path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("genome_id\tsplit\tgroup\tsource\tspecies\tgenus\tfamily\n")
        for g in genomes:
            fh.write("\t".join((g["genome_id"], g["split"], g["group"], SOURCE, g["species"], g["genus"], g["family"])) + "\n")


def read_labelled(labels_path: Path) -> dict:
    labelled = defaultdict(set)
    with open(labels_path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            labelled[r["genome_id"]].add(r["protein_id"])
    return labelled


def write_proteins(genomes, labelled, proteins_path: Path, fasta_path: Path) -> dict:
    """Write proteins.tsv and the unique-sequence FASTA; returns ``{(genome, protein): seq_id}``."""
    ids, unique = {}, {}
    with open(proteins_path, "w") as fh:
        fh.write("genome_id\tprotein_id\tkind\tgenome_protein_id\n")
        for g in genomes:
            seqs = PL.fasta_sequences(Path(g["directory"]) / "protein.faa")
            for pid, seq in sorted(seqs.items()):
                sid = seq_id(seq)
                ids[(g["genome_id"], pid)] = sid
                unique[sid] = seq
                kind = "target" if pid in labelled.get(g["genome_id"], ()) else "context"
                fh.write(f"{g['genome_id']}\t{sid}\t{kind}\t{pid}\n")
    with open(fasta_path, "w") as fh:
        for sid, seq in sorted(unique.items()):
            fh.write(f">{sid}\n{seq}\n")
    return ids


def run_mmseqs(fasta: Path, work: Path, min_id=0.3, cov=0.8, threads=4) -> dict:
    """``easy-cluster`` over ``fasta``; returns the command, version and ``{member: representative}``."""
    work.mkdir(parents=True, exist_ok=True)
    prefix, tmp = work / "clusters", work / "tmp"
    cmd = ["mmseqs", "easy-cluster", str(fasta), str(prefix), str(tmp), "--min-seq-id", str(min_id), "-c", str(cov),
           "--cov-mode", "0", "--cluster-mode", "2", "--cluster-reassign", "1", "-s", "7.5", "--max-seqs", "1000",
           "--threads", str(min(threads, 4))]
    log = subprocess.run(cmd, capture_output=True, text=True)
    (work / "mmseqs.log").write_text(log.stdout + log.stderr)
    if log.returncode:
        raise CohortError(f"mmseqs failed ({log.returncode}); see {work / 'mmseqs.log'}")
    version = subprocess.run(["mmseqs", "version"], capture_output=True, text=True).stdout.strip()
    members = {}
    with open(f"{prefix}_cluster.tsv") as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")
            if member in members:
                raise CohortError(f"cluster membership repeats {member}")
            members[member] = rep
    return {"command": " ".join(cmd), "version": version, "members": members}


def write_clusters(members: dict, unique_ids: set, ids: dict, labelled, genomes, out_dir: Path) -> dict:
    missing = sorted(unique_ids - set(members))
    if missing:
        raise CohortError(f"clustering lost {len(missing)} sequence(s), e.g. {missing[:3]}")
    with open(out_dir / "protein_clusters.tsv", "w") as fh:
        fh.write("protein_id\tcluster\n")
        for sid in sorted(unique_ids):
            fh.write(f"{sid}\t{members[sid]}\n")
    n = 0
    with open(out_dir / "example_clusters.tsv", "w") as fh:
        fh.write("protein_id\tgenome_id\tcluster\n")
        for g in genomes:
            for pid in sorted(labelled.get(g["genome_id"], ())):
                fh.write(f"{pid}\t{g['genome_id']}\t{members[ids[(g['genome_id'], pid)]]}\n")
                n += 1
    return {"unique_sequences": len(unique_ids), "clusters": len(set(members[s] for s in unique_ids)), "labelled_examples": n}


def family_overlap_preview(example_clusters: Path, genomes) -> dict:
    """What ``build-table`` masking would remove, by its precedence (test > calibration > train).

    Counts labelled PROTEINS (sequence identities) per split, and those whose family also holds a labelled protein
    of a higher-precedence split. Test proteins are never dropped; they are reported as ``test_with_family_in_train``.
    """
    split_of = {g["genome_id"]: g["split"] for g in genomes}
    members = defaultdict(lambda: defaultdict(set))          # cluster -> split -> proteins
    with open(example_clusters, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            members[r["cluster"]][split_of[r["genome_id"]]].add(r["protein_id"])
    per = {s: {"labelled": 0, "masked": 0} for s in ("train", "calibration", "test")}
    higher = {"train": ("calibration", "test"), "calibration": ("test",), "test": ()}
    test_in_train = set()
    for by_split in members.values():
        for split, prots in by_split.items():
            per[split]["labelled"] += len(prots)
            if any(by_split.get(h) for h in higher[split]):
                per[split]["masked"] += len(prots)
        if by_split.get("test") and by_split.get("train"):
            test_in_train |= by_split["test"]
    per["test"]["with_family_in_train_kept"] = len(test_in_train)
    for split in ("train", "calibration"):
        per[split]["kept"] = per[split]["labelled"] - per[split]["masked"]
    per["test"]["kept"] = per["test"]["labelled"]
    return per


def prepare(provenance_path: Path, labels_path: Path, out_dir: Path, selection_path=SELECTION, threads=4,
            disjoint="genome_id,group,genus,family", cluster=True) -> dict:
    prov = json.loads(provenance_path.read_text())
    selection = json.loads(Path(selection_path).read_text())
    genomes = eligible_genomes(prov, selection)
    out_dir.mkdir(parents=True, exist_ok=True)
    labelled = read_labelled(labels_path)
    write_splits(genomes, out_dir / "splits.tsv")
    ids = write_proteins(genomes, labelled, out_dir / "proteins.tsv", out_dir / "unique_proteins.faa")
    unique_ids = set(ids.values())
    summary = {"genomes": len(genomes), "per_split": dict(Counter(g["split"] for g in genomes)),
               "excluded_unresolved": sorted(t for t, r in prov["organisms"].items() if r.get("assembly", {}).get("status") != "resolved"),
               "excluded_price": prov.get("excluded_by_price_audit", []),
               "proteins": len(ids), "targets": sum(len(v) for v in labelled.values())}
    if cluster:
        clu = run_mmseqs(out_dir / "unique_proteins.faa", out_dir / "mmseqs", threads=threads)
        summary["clustering"] = {"command": clu["command"], "version": clu["version"],
                                 **write_clusters(clu["members"], unique_ids, ids, labelled, genomes, out_dir)}
    if cluster:
        summary["family_masking_preview"] = family_overlap_preview(out_dir / "example_clusters.tsv", genomes)
    audit_cmd = [sys.executable, "-I", str(HERE / "prepare_training_data.py"), "audit-splits", "--splits", str(out_dir / "splits.tsv"),
                 "--proteins", str(out_dir / "proteins.tsv"), "--disjoint", disjoint, "--out", str(out_dir / "split_audit.json")]
    if cluster:
        audit_cmd += ["--protein-clusters", str(out_dir / "protein_clusters.tsv"), "--cluster-policy", "mask"]
    proc = subprocess.run(audit_cmd, capture_output=True, text=True)
    summary["audit_exit_status"] = proc.returncode
    summary["audit_stderr"] = proc.stderr.strip()[-800:]
    (out_dir / "cohort_summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--provenance", type=Path, required=True)
    ap.add_argument("--labels", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--no-cluster", action="store_true")
    args = ap.parse_args(argv)
    s = prepare(args.provenance, args.labels, args.out_dir, threads=args.threads, cluster=not args.no_cluster)
    print(json.dumps(s, indent=1))
    return s["audit_exit_status"]


if __name__ == "__main__":
    sys.exit(main())
