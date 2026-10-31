#!/usr/bin/env python3
"""Write the WP4 fenced gapsmith search manifests (design §6: every genome is searched with the fenced reference).

One ``gspa.cohort-pathway-search/1`` manifest per (genome set, superkingdom), consumed unchanged by
``observe_gapsmith_cohort.py`` (one genome per invocation, ``--index``) and, after the searches, by
``prepare_wp4_batches.py gaps|context --search-root SET=ROOT`` and ``audit_cohort_pathways.audit``.

``--set NAME=KIND:PATH`` (repeatable) names one genome set:

``cohort:DIR``        ``DIR/<genome_id>.faa`` -- the 10406 cohort proteomes. Each file must hash to the sha256 of
                      its row in the ``--cohort-manifest`` search manifests (cohort-pathway-search-v1 and
                      -archaea-v2), which also give ``taxid`` and ``taxonomy`` (Bacteria|Archaea). The directory
                      must hold exactly the cohort of those manifests. Headers are the original bytes (the
                      consumers require the 10406 proteome bytes), so they are not rewritten or checked.
``peers:DIR``         ``DIR/<accession>/{protein.faa,taxonomy.json,provenance.json}`` from ``acquire_peer_genomes.py``.
                      ``protein.faa`` must hash to ``provenance.json`` ``file_sha256s`` (and the acquisition must be
                      ``verified`` when ``DIR/acquisition.json`` exists). ``taxid`` = species taxid. Superkingdom:
                      ``taxonomy.json`` ``superkingdom``/``domain`` when present, else the ``domain``/``superkingdom``
                      rank of the organism in the hash-checked NCBI Taxonomy cache named by ``provenance.json``;
                      it must agree with the superkingdom of the cohort target in ``selection_rows``.
``perturbation:DIR``  a ``prepare_wp4_perturbations.py`` output (``DIR/gapsmith_inputs.tsv``) or a mode directory
                      of per-genome outputs (``DIR/*/gapsmith_inputs.tsv``). Each copy's proteome must hash to its
                      ``gapsmith_inputs.tsv`` row and to the output's ``receipt.json``; taxid and superkingdom are
                      inherited from ``source_genome`` (cohort manifests, or a peers set).

Peer and copy FASTA headers must be bare ids (``>id`` with no whitespace). A genome in two sets, a missing or
unsupported superkingdom, a proteome differing from its source receipt, or a superkingdom without a verified
fenced reference all fail before anything is written.

Outputs (``--out`` must not exist): ``<set>/<superkingdom>/{manifest.json, staging-hashes.json,
proteomes/<genome_id>.faa}`` (byte-identical copies named by genome id: gapsmith names its tables after the
proteome stem and the consumers expect ``<genome_id>-all-Reactions.tbl``), ``search-plan.tsv`` and
``receipt.json``. Paths written into a manifest (proteomes, ``output_root``, runner scripts) are rewritten with
``--path-map LOCAL=REMOTE``; ``--binary``, ``--diamond`` and the references are given as cluster paths with
their sha256. ``output_root`` is ``<mapped set/superkingdom dir>/observations`` so the staged directory is a
complete search root. No label file is read.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = "gspa.cohort-pathway-search/1"
RECEIPT_SCHEMA = "gspa.wp4-search-manifests/1"
FENCE_SCHEMA = "gspa.wp4-fenced-gapsmith-reference/1"
FENCE_VERIFY_SCHEMA = "gspa.wp4-fenced-gapsmith-reference-verification/1"
SUPERKINGDOMS = ("Bacteria", "Archaea")
KINDS = ("cohort", "peers", "perturbation")
PLAN_FIELDS = ["set", "superkingdom", "genome_id", "manifest", "proteins"]
INPUT_FIELDS = ["genome_id", "proteome", "sha256", "source_genome", "perturbation"]
RUNNERS = ("observe_gapsmith.py", "observe_gapsmith_cohort.py")
HEX64 = re.compile(r"[0-9a-f]{64}")
BARE_HEADER = re.compile(r">\S+")
LABEL_SUFFIXES = ("labels.tsv", "swissprot_matches.tsv", "label_policy.json", "label_coverage.json")
LABEL_PREFIXES = ("uniprot_sprot.dat", "role_truth")
HERE = Path(__file__).resolve().parent


class ManifestError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ManifestError(message)


def digest(path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def save_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def refuse_label_path(path):
    name = Path(path).name
    require(not (name.endswith(LABEL_SUFFIXES) or name.startswith(LABEL_PREFIXES)),
            f"refusing a label file as input: {path}")


def normalize_superkingdom(value, genome_id):
    if isinstance(value, dict):
        value = value.get("name")
    require(isinstance(value, str) and value.strip(), f"superkingdom missing for {genome_id}")
    normalized = value.strip().capitalize()
    require(normalized in SUPERKINGDOMS, f"unsupported superkingdom {value!r} for {genome_id}")
    return normalized


# --------------------------------------------------------------------------------------------- path maps
def parse_path_maps(values):
    result = []
    for value in values or ():
        require("=" in value, "invalid --path-map (expected LOCAL=REMOTE): " + value)
        local, remote = value.split("=", 1)
        require(local and remote, "invalid empty --path-map endpoint: " + value)
        result.append((local.rstrip("/") or "/", remote.rstrip("/") or "/"))
    return result


def map_path(path, path_maps):
    """Longest-prefix rewrite of the absolute path, then of its resolved form; (text, mapped?)."""
    candidates = [os.path.abspath(str(path)), str(Path(path).resolve())]
    prefixes = [(local, remote) for local, remote in path_maps]
    prefixes += [(str(Path(local).resolve()), remote) for local, remote in path_maps]  # symlinked local roots
    for text in candidates:
        for local, remote in sorted(prefixes, key=lambda pair: len(pair[0]), reverse=True):
            if text == local or text.startswith(local.rstrip("/") + "/"):
                return remote + text[len(local.rstrip("/")):], True
    return candidates[0], False


# ------------------------------------------------------------------------------------------------ inputs
class Inputs:
    """Every file read, with its sha256 (provenance_inputs and the receipt)."""

    def __init__(self):
        self.hashes = {}

    def record(self, path):
        path = Path(path).resolve()
        refuse_label_path(path)
        require(path.is_file(), f"missing input: {path}")
        key = str(path)
        if key not in self.hashes:
            self.hashes[key] = digest(path)
        return self.hashes[key]

    def json(self, path):
        self.record(path)
        return read_json(path)

    def verify_unchanged(self):
        for path, expected in self.hashes.items():
            require(digest(path) == expected, f"input changed while preparing manifests: {path}")


def fasta_headers(path):
    with Path(path).open() as handle:
        return [line.rstrip("\n\r") for line in handle if line.startswith(">")]


def check_bare_headers(path, genome_id):
    headers = fasta_headers(path)
    require(headers, f"empty proteome for {genome_id}: {path}")
    bad = [h for h in headers if not BARE_HEADER.fullmatch(h)]
    if bad:
        raise ManifestError(f"{genome_id}: FASTA header carries more than a bare id (product name?): {bad[0][:80]!r}")
    ids = [h[1:] for h in headers]
    require(len(ids) == len(set(ids)), f"{genome_id}: duplicate protein ids in {path}")
    return len(headers)


def cohort_index(manifests, inputs):
    """genome -> {sha256, taxid, taxonomy, proteins} from the 10406 cohort search manifests (must agree)."""
    index = {}
    for path in manifests or ():
        data = inputs.json(path)
        require(data.get("schema") == SCHEMA, f"not a cohort search manifest: {path}")
        for row in data.get("genomes", []):
            gid = row.get("genome_id")
            require(isinstance(gid, str) and gid, f"invalid cohort manifest row in {path}")
            record = {"sha256": row.get("sha256"), "taxid": str(row.get("taxid")),
                      "taxonomy": normalize_superkingdom(row.get("taxonomy"), gid),
                      "proteins": row.get("proteins")}
            require(isinstance(record["sha256"], str) and HEX64.fullmatch(record["sha256"]),
                    f"invalid sha256 for {gid} in {path}")
            if gid in index:
                previous = {k: v for k, v in index[gid].items() if k != "manifests"}
                require(previous == record, f"cohort manifests disagree on {gid}")
                index[gid]["manifests"].append(str(Path(path).resolve()))
            else:
                index[gid] = {**record, "manifests": [str(Path(path).resolve())]}
    return index


class TaxonomyCache:
    def __init__(self, inputs):
        self.inputs = inputs
        self.parsed = {}

    def domain(self, source, taxid):
        path = Path(source.get("cache_file", ""))
        require(path.is_file(), f"NCBI taxonomy cache missing: {path}")
        sha = self.inputs.record(path)
        require(source.get("sha256") in (None, sha), f"NCBI taxonomy cache changed: {path}")
        key = str(path.resolve())
        if key not in self.parsed:
            domains = {}
            for taxon in ET.fromstring(path.read_bytes()).findall("Taxon"):
                for parent in taxon.findall("LineageEx/Taxon"):
                    if parent.findtext("Rank") in ("domain", "superkingdom"):
                        domains[taxon.findtext("TaxId")] = parent.findtext("ScientificName")
            self.parsed[key] = domains
        return self.parsed[key].get(str(taxid))


def collect_cohort(name, directory, index, inputs):
    directory = Path(directory).resolve()
    require(directory.is_dir(), f"cohort proteome directory missing: {directory}")
    require(index, f"set {name}: cohort sets require --cohort-manifest")
    files = {p.name[:-len(".faa")]: p for p in sorted(directory.glob("*.faa"))}
    require(set(files) == set(index), f"set {name}: cohort proteomes differ from the cohort manifests "
            f"(only in dir: {sorted(set(files) - set(index))[:3]}, missing: {sorted(set(index) - set(files))[:3]})")
    rows = []
    for gid, path in files.items():
        sha = inputs.record(path)
        require(sha == index[gid]["sha256"], f"{gid}: cohort proteome differs from its search manifest receipt")
        proteins = len(fasta_headers(path))
        require(index[gid]["proteins"] in (None, proteins), f"{gid}: protein count differs from cohort manifest")
        rows.append({"genome_id": gid, "source": path, "sha256": sha, "proteins": proteins,
                     "taxid": index[gid]["taxid"], "taxonomy": index[gid]["taxonomy"],
                     "origin": {"kind": "cohort", "receipts": index[gid]["manifests"],
                                "headers": "original 10406 bytes (not rewritten)"}})
    return rows


def peer_records(directory, inputs, index):
    """accession -> peer row (shared by peers sets and copies inheriting from a peer)."""
    directory = Path(directory).resolve()
    require(directory.is_dir(), f"peers directory missing: {directory}")
    acquisition = {}
    if (directory / "acquisition.json").is_file():
        acquisition = inputs.json(directory / "acquisition.json").get("genomes", {})
    taxonomy_cache = TaxonomyCache(inputs)
    rows = {}
    for sub in sorted(p for p in directory.iterdir() if p.is_dir()):
        gid = sub.name
        provenance_path = sub / "provenance.json"
        require(provenance_path.is_file(), f"peer {gid}: incomplete acquisition (no provenance.json)")
        if acquisition:
            require(acquisition.get(gid, {}).get("status") == "verified", f"peer {gid}: acquisition not verified")
        provenance = inputs.json(provenance_path)
        taxonomy = inputs.json(sub / "taxonomy.json")
        proteome = sub / "protein.faa"
        sha = inputs.record(proteome)
        require(sha == provenance.get("file_sha256s", {}).get("protein.faa"),
                f"peer {gid}: protein.faa differs from its acquisition receipt")
        require(taxonomy.get("accession") in (None, gid), f"peer {gid}: taxonomy.json names another accession")
        proteins = check_bare_headers(proteome, gid)
        declared = taxonomy.get("superkingdom") or taxonomy.get("domain")
        evidence = []
        if declared:
            evidence.append(("taxonomy.json", normalize_superkingdom(declared, gid)))
        else:
            for source in provenance.get("taxonomy_sources") or ():
                found = taxonomy_cache.domain(source, taxonomy.get("organism_taxid"))
                if found:
                    evidence.append(("ncbi-taxonomy-cache:" + str(source.get("cache_file")),
                                     normalize_superkingdom(found, gid)))
        targets = sorted({r.get("target_accession") for r in provenance.get("selection_rows") or ()
                          if r.get("target_accession") in index})
        for target in targets:
            evidence.append(("cohort-target:" + target, index[target]["taxonomy"]))
        values = {value for _, value in evidence}
        require(values, f"peer {gid}: no superkingdom (taxonomy.json, NCBI cache and cohort target all silent)")
        require(len(values) == 1, f"peer {gid}: conflicting superkingdom evidence {evidence}")
        taxid = taxonomy.get("species_taxid") or (taxonomy.get("species") or {}).get("taxid")
        require(taxid, f"peer {gid}: taxonomy.json has no species taxid")
        rows[gid] = {"genome_id": gid, "source": proteome, "sha256": sha, "proteins": proteins,
                     "taxid": str(taxid), "taxonomy": values.pop(),
                     "origin": {"kind": "peer", "receipt": str(provenance_path),
                                "superkingdom_evidence": [e for e, _ in evidence]}}
    require(rows, f"no peer genomes in {directory}")
    return rows


def perturbation_outputs(directory):
    directory = Path(directory).resolve()
    if (directory / "gapsmith_inputs.tsv").is_file():
        return [directory]
    subdirs = sorted(p for p in directory.iterdir() if p.is_dir())
    require(subdirs, f"no gapsmith_inputs.tsv in {directory} or its subdirectories")
    incomplete = [p.name for p in subdirs
                  if not ((p / "gapsmith_inputs.tsv").is_file() and (p / "receipt.json").is_file())]
    require(not incomplete, f"incomplete perturbation outputs in {directory} (no gapsmith_inputs.tsv/receipt.json): "
            + ", ".join(incomplete[:5]))
    return subdirs


def collect_perturbation(name, directory, index, peers, inputs):
    mode_dir = Path(directory).resolve()
    rows = []
    for output in perturbation_outputs(mode_dir):
        receipt = inputs.json(output / "receipt.json")
        outputs = receipt.get("outputs") or {}
        with (output / "gapsmith_inputs.tsv").open(newline="") as handle:
            inputs.record(output / "gapsmith_inputs.tsv")
            reader = csv.DictReader(handle, delimiter="\t")
            require(reader.fieldnames == INPUT_FIELDS, f"unexpected gapsmith_inputs.tsv header in {output}")
            listed = list(reader)
        require(outputs.get("gapsmith_inputs.tsv") in (None, inputs.hashes[str(output / "gapsmith_inputs.tsv")]),
                f"gapsmith_inputs.tsv differs from receipt.json in {output}")
        for row in listed:
            gid, relative = row["genome_id"], row["proteome"]
            require(gid and not Path(relative).is_absolute() and ".." not in Path(relative).parts,
                    f"invalid gapsmith input row in {output}: {row}")
            proteome = output / relative
            sha = inputs.record(proteome)
            require(sha == row["sha256"], f"{gid}: proteome differs from gapsmith_inputs.tsv")
            require(outputs.get(Path(relative).as_posix()) == sha,
                    f"{gid}: proteome differs from (or is absent in) its perturbation receipt")
            source_id = row["source_genome"]
            source = index.get(source_id) or (peers or {}).get(source_id)
            require(source, f"{gid}: source genome {source_id} is neither a cohort nor a peer genome")
            rows.append({"genome_id": gid, "source": proteome, "sha256": sha,
                         "proteins": check_bare_headers(proteome, gid),
                         "taxid": str(source["taxid"]), "taxonomy": source["taxonomy"],
                         "origin": {"kind": "copy", "source_genome": source_id,
                                    "perturbation": row["perturbation"],
                                    "subdirectory": output.relative_to(mode_dir).as_posix() or ".",
                                    "receipt": str(output / "receipt.json")}})
    return rows


# ---------------------------------------------------------------------------------------------- fence
def fence_block(superkingdom, reference, report_path, verify_path, inputs):
    require(reference, f"--reference-{superkingdom.lower()} is required (genomes of that superkingdom)")
    require(report_path and verify_path,
            f"--fence-report-{superkingdom.lower()} and --fence-verify-{superkingdom.lower()} are required")
    report = inputs.json(report_path)
    require(report.get("schema") == FENCE_SCHEMA and report.get("status") == "built"
            and report.get("label_blind") is True, f"not a built label-blind fenced reference report: {report_path}")
    verify = inputs.json(verify_path)
    report_sha = inputs.hashes[str(Path(report_path).resolve())]
    require(verify.get("schema") == FENCE_VERIFY_SCHEMA and verify.get("passed") is True
            and verify.get("report_sha256") == report_sha,
            f"fenced reference verification did not pass for this report: {verify_path}")
    reference = reference.rstrip("/")
    return {"superkingdom": superkingdom, "reference": reference,
            "report": str(Path(reference).parent / "fenced_reference_report.json"),
            "report_sha256": report_sha,
            "verify_report": str(Path(reference).parent / "verify" / "verify_report.json"),
            "verify_report_sha256": inputs.hashes[str(Path(verify_path).resolve())],
            "local_report": str(Path(report_path).resolve()), "local_verify_report": str(Path(verify_path).resolve()),
            "thresholds": report.get("thresholds"), "counts": report.get("counts", {}).get("removed_by_reason")}


def tool_hash(name, sha, local, inputs):
    require(bool(sha) != bool(local), f"give exactly one of --{name}-sha256 and --{name}-local")
    if local:
        return inputs.record(local)
    require(HEX64.fullmatch(sha), f"--{name}-sha256 must be 64 lowercase hex characters")
    return sha


def parse_sets(values):
    result = {}
    for value in values or ():
        require("=" in value and ":" in value.split("=", 1)[1], f"invalid --set (NAME=KIND:PATH): {value}")
        name, spec = value.split("=", 1)
        kind, path = spec.split(":", 1)
        require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name), f"invalid set name: {name}")
        require(kind in KINDS, f"unknown set kind {kind!r} (expected {', '.join(KINDS)})")
        require(name not in result, f"duplicate set name: {name}")
        result[name] = (kind, Path(path))
    require(result, "at least one --set NAME=KIND:PATH is required")
    return result


# ----------------------------------------------------------------------------------------------- main
def prepare(sets, out, cohort_manifests=(), peers_dir=None, binary=None, binary_sha256=None, binary_local=None,
            diamond=None, diamond_sha256=None, diamond_local=None, references=None, fence_reports=None,
            fence_verifies=None, threads=4, path_maps=(), selection_policy=None):
    out = Path(out).resolve()
    require(not out.exists(), f"refusing to overwrite an existing output directory: {out}")
    require(isinstance(threads, int) and threads >= 1, "--threads must be a positive integer")
    require(binary and diamond, "--binary and --diamond (cluster paths) are required")
    sets = parse_sets(sets) if not isinstance(sets, dict) else sets
    path_maps = list(path_maps)
    inputs = Inputs()
    index = cohort_index(cohort_manifests, inputs)
    peer_sets = [path for kind, path in sets.values() if kind == "peers"]
    peers_source = peers_dir or (peer_sets[0] if len(peer_sets) == 1 else None)
    peers = peer_records(peers_source, inputs, index) if peers_source else {}
    references, fence_reports, fence_verifies = references or {}, fence_reports or {}, fence_verifies or {}

    collected, owner = {}, {}
    for name, (kind, path) in sets.items():
        if kind == "cohort":
            rows = collect_cohort(name, path, index, inputs)
        elif kind == "peers":
            rows = list(peers.values()) if Path(path).resolve() == Path(peers_source).resolve() \
                else list(peer_records(path, inputs, index).values())
        else:
            rows = collect_perturbation(name, path, index, peers, inputs)
        require(rows, f"set {name} has no genomes")
        for row in rows:
            gid = row["genome_id"]
            require(gid not in owner, f"genome {gid} appears in two sets ({owner.get(gid)}, {name})")
            owner[gid] = name
        collected[name] = rows

    runner_hashes = {}
    unmapped = []
    for script in RUNNERS:
        local = HERE / script
        runner_hashes[script] = (map_path(local, path_maps), inputs.record(local))
    binary_hash = tool_hash("binary", binary_sha256, binary_local, inputs)
    diamond_hash = tool_hash("diamond", diamond_sha256, diamond_local, inputs)
    sources = {binary: binary_hash, diamond: diamond_hash}
    for script, ((remote, mapped), sha) in runner_hashes.items():
        sources[remote] = sha
        if not mapped:
            unmapped.append(remote)
    require(len(sources) == 4, "binary, diamond and runner scripts must be four distinct paths")

    fences = {}
    for superkingdom in sorted({row["taxonomy"] for rows in collected.values() for row in rows}):
        fences[superkingdom] = fence_block(superkingdom, references.get(superkingdom),
                                           fence_reports.get(superkingdom), fence_verifies.get(superkingdom),
                                           inputs)

    # Every check passed: write.
    plan, manifests = [], {}
    for name, rows in collected.items():
        for superkingdom in SUPERKINGDOMS:
            members = sorted((r for r in rows if r["taxonomy"] == superkingdom), key=lambda r: r["genome_id"])
            if not members:
                continue
            root = out / name / superkingdom.lower()
            (root / "proteomes").mkdir(parents=True)
            genomes, provenance = [], {}
            for row in members:
                staged = root / "proteomes" / f"{row['genome_id']}.faa"
                shutil.copyfile(row["source"], staged)
                require(digest(staged) == row["sha256"], f"{row['genome_id']}: staged copy differs")
                remote, mapped = map_path(staged, path_maps)
                if not mapped:
                    unmapped.append(remote)
                genomes.append({"genome_id": row["genome_id"], "proteins": row["proteins"], "proteome": remote,
                                "sha256": row["sha256"], "taxid": row["taxid"], "taxonomy": superkingdom})
                provenance[str(Path(row["source"]).resolve())] = row["sha256"]
                if row["origin"].get("receipt"):
                    provenance[row["origin"]["receipt"]] = inputs.hashes[row["origin"]["receipt"]]
                for receipt in row["origin"].get("receipts", ()):
                    provenance[receipt] = inputs.hashes[receipt]
            output_root, mapped = map_path(root / "observations", path_maps)
            if not mapped:
                unmapped.append(output_root)
            manifest = {
                "schema": SCHEMA, "set": name, "superkingdom": superkingdom,
                "binary": binary, "diamond": diamond, "reference": fences[superkingdom]["reference"],
                "threads": threads, "output_root": output_root, "genomes": genomes, "sources": sources,
                "provenance_inputs": dict(sorted(provenance.items())),
                "fence": {k: v for k, v in fences[superkingdom].items() if not k.startswith("local_")},
                "origins": {r["genome_id"]: r["origin"] for r in members},
                "selection_policy": selection_policy or (
                    f"WP4 set {name}: every {superkingdom} genome of the set, searched with the fenced gapsmith "
                    "reference (design §6); no function labels consulted."),
                "generated_by": {"script": str(Path(__file__).resolve()), "sha256": digest(Path(__file__))},
            }
            save_json(root / "manifest.json", manifest)
            staging = {"manifest.json": digest(root / "manifest.json")}
            staging.update({f"proteomes/{g['genome_id']}.faa": g["sha256"] for g in genomes})
            save_json(root / "staging-hashes.json", staging)
            manifests[f"{name}/{superkingdom.lower()}"] = {
                "manifest": f"{name}/{superkingdom.lower()}/manifest.json",
                "manifest_sha256": staging["manifest.json"], "genomes": len(genomes),
                "output_root": output_root}
            plan += [[name, superkingdom, g["genome_id"], f"{name}/{superkingdom.lower()}/manifest.json",
                      str(g["proteins"])] for g in genomes]
    require(not (path_maps and unmapped), "paths not covered by any --path-map: " + ", ".join(sorted(set(unmapped))[:5]))
    with (out / "search-plan.tsv").open("w") as handle:
        handle.write("\t".join(PLAN_FIELDS) + "\n")
        handle.writelines("\t".join(row) + "\n" for row in plan)
    inputs.verify_unchanged()
    outputs = {p.relative_to(out).as_posix(): digest(p) for p in sorted(out.rglob("*"))
               if p.is_file() and "/proteomes/" not in "/" + p.relative_to(out).as_posix()}
    receipt = {
        "schema": RECEIPT_SCHEMA, "status": "prepared", "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": {"path": str(Path(__file__).resolve()), "sha256": digest(Path(__file__))},
        "sets": {name: {"kind": kind, "path": str(Path(path).resolve()), "genomes": len(collected[name])}
                 for name, (kind, path) in sets.items()},
        "manifests": manifests, "fences": fences, "sources": sources, "threads": threads,
        "path_maps": [{"local": local, "remote": remote} for local, remote in path_maps],
        "unmapped_paths": sorted(set(unmapped)),
        "inputs": dict(sorted(inputs.hashes.items())), "outputs": outputs,
        "label_blind": True,
    }
    save_json(out / "receipt.json", receipt)
    return receipt


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", action="append", default=[], metavar="NAME=KIND:PATH", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--cohort-manifest", action="append", default=[], type=Path,
                   help="10406 cohort search manifest (repeatable): cohort sha256/taxid/taxonomy")
    p.add_argument("--peers-dir", type=Path, help="peer acquisitions, for copies of peers (default: the peers set)")
    p.add_argument("--binary", required=True, help="cluster path of gapsmith")
    p.add_argument("--binary-sha256")
    p.add_argument("--binary-local", type=Path, help="local copy to hash instead of --binary-sha256")
    p.add_argument("--diamond", required=True, help="cluster path of DIAMOND")
    p.add_argument("--diamond-sha256")
    p.add_argument("--diamond-local", type=Path)
    for sk in ("bacteria", "archaea"):
        p.add_argument(f"--reference-{sk}", help=f"cluster path of the fenced {sk} gapsmith reference")
        p.add_argument(f"--fence-report-{sk}", type=Path, help="local copy of its fenced_reference_report.json")
        p.add_argument(f"--fence-verify-{sk}", type=Path, help="local copy of its verify/verify_report.json")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--path-map", action="append", default=[], metavar="LOCAL=REMOTE")
    p.add_argument("--selection-policy")
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    try:
        receipt = prepare(
            a.set, a.out, a.cohort_manifest, a.peers_dir, a.binary, a.binary_sha256, a.binary_local,
            a.diamond, a.diamond_sha256, a.diamond_local,
            {"Bacteria": a.reference_bacteria, "Archaea": a.reference_archaea},
            {"Bacteria": a.fence_report_bacteria, "Archaea": a.fence_report_archaea},
            {"Bacteria": a.fence_verify_bacteria, "Archaea": a.fence_verify_archaea},
            a.threads, parse_path_maps(a.path_map), a.selection_policy)
    except ManifestError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({k: v["genomes"] for k, v in receipt["manifests"].items()}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
