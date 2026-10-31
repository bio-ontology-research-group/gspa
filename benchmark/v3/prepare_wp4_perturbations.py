#!/usr/bin/env python3
"""Prepare the WP4 perturbation inputs (design §4): removal, fragmentation, shuffles and duplicates.

``removal``   one copy of a target genome per removed protein (protein, claims, layout row and family row
              removed). Removed: every strict-EC-labelled protein plus, per genome, an equal number of random
              (seed 42) unlabelled proteins with a complete EC claim. This is the only step that reads labels,
              and only to choose the removed proteins (the documented exception of design §6); copy names carry
              no label information, the label-derived category is kept in ``removal_selection.tsv``.
``fragment``  target genomes cut at contig lengths drawn i.i.d. from fixed log-normals (medians 100/25/5 kb,
              sigma 1, min 500 bp, seeds 1-2). One draw stream per (median, seed), identical for every target.
              CDSs crossing a cut are dropped from the proteome, claims, layout and family rows (the
              ``prepare_sequence_loss`` rule); pieces are open contigs of a DRAFT assembly.
``shuffle``   S1: in peer genomes only, family assignments are permuted among the placed proteins of one stratum
              (contig, window gene count within +-halfWidthBp, window state, contig closed flag: the engine's
              matched-background bucket). Family multisets per genome and stratum, layouts and claims are kept;
              cohort (target) rows stay byte-identical. S2: degree-preserving swaps of ``reaction_id`` between
              ``pathway-groups.tsv`` groups (group sizes and per-reaction group counts kept).
``duplicate`` per family block: (i) exact copy of the top congeneric peer in its own species cluster, (ii) the
              same copy as a new species cluster, (iii) the near-duplicate relabelled as its own species cluster,
              (iv) an exact copy of the target relabelled as a congeneric species, (v) the block without the
              near-duplicate.

Arm variant names follow the analysis conventions: ``removal``, ``fragment:<tier>:<seed>`` (tier = median, e.g.
``25k``), ``S1:<seed>``, ``S2:<seed>`` and ``dup:i`` .. ``dup:v``; each output lists them in ``variants.tsv`` with
``sealed_genomes`` (the target plus every copy of it as ``copy=target``). ``removal_selection.tsv`` is the T8
removals manifest (copy_genome, target_genome, removed_protein, class, function_id, function_source): one row per
strict EC label of a labelled protein and per complete EC claim of a random one; it is the only output that carries
label-derived content and must not be read by any feature-building step.

Copies of a target share its clade and strain cluster (removal, fragment) so that the holdout keeps them out of
each other's peers. Removal and fragment copies need a fresh gapsmith search: their ``plausibleReactions`` are
dropped from the metadata and the proteomes are listed in ``gapsmith_inputs.tsv``. Every output directory gets a
``receipt.json`` with input and output hashes, seeds and counts. Outputs carry bare protein ids and no product
names; ``verify`` re-checks a receipt.
"""
from __future__ import annotations

import argparse
import bisect
import copy
import csv
import hashlib
import json
import math
import random
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = "gspa.wp4-perturbation/1"
LAYOUT_FIELDS = ["protein_id", "contig", "start", "end", "strand"]
ORTHO_FIELDS = ["genome_id", "protein_id", "orthogroup_id", "subfamily_id", "weight"]
GROUP_FIELDS = ["pathway_id", "group_id", "reaction_id", "source"]
REMOVAL_SEED = 42
FRAGMENT_MEDIANS = (100000, 25000, 5000)
FRAGMENT_SIGMA = 1.0
FRAGMENT_MIN_BP = 500
FRAGMENT_SEEDS = (1, 2)
S1_SEEDS = (1, 2, 3, 4, 5)
S2_SEEDS = (1, 2, 3)
HALF_WIDTH_BP = 5000          # CohortConfig.halfWidthBp default (engine defaults are not tuned, D5)
SWAP_FACTOR = 10
PRODUCT_KEYS = {"product", "product_name", "description", "protein_name", "name", "gene_product"}
LABEL_SUFFIXES = ("labels.tsv", "swissprot_matches.tsv", "label_policy.json", "label_coverage.json")
LABEL_PREFIXES = ("uniprot_sprot.dat", "role_truth")
LABEL_EXCEPTION = ("prepare_wp4_perturbations.py removal reads labels only to choose the removed proteins "
                   "(design §6 documented exception); no label is written into a genome, claim or family file.")
DUPLICATE_VARIANTS = ("i_top_c_copy_same_cluster", "ii_top_c_copy_new_species", "iii_nd_own_cluster",
                      "iv_target_copy_congeneric", "v_without_nd")
VARIANT_FIELDS = ["variant", "target", "sealed_genomes", "path"]
REMOVAL_FIELDS = ["copy_genome", "target_genome", "removed_protein", "class", "function_id", "function_source"]


def duplicate_variant_name(slug):
    """``dup:<kind>`` arm variant (kind = roman numeral of design §4) for a duplicate directory slug."""
    return "dup:" + slug.split("_", 1)[0]


def fragment_tier(median):
    return f"{median // 1000}k" if median % 1000 == 0 else f"{median}bp"


def sealed(target, copies=()):
    """``sealed_genomes`` value of an arm: the target and every copy of it as ``copy=target``."""
    return ",".join([target] + [f"{c}={target}" for c in copies])


class PerturbationError(ValueError):
    """An input, consistency or output contract was not satisfied."""


def require(condition, message):
    if not condition:
        raise PerturbationError(message)


# ----------------------------------------------------------------------------------------------- helpers
def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def is_label_path(path) -> bool:
    path = Path(path)
    name = path.name.lower()
    return (any(name.endswith(s) for s in LABEL_SUFFIXES) or any(name.startswith(p) for p in LABEL_PREFIXES)
            or path.parent.name.lower() == "labels")


def refuse_label_paths(paths):
    for path in paths:
        if path is not None and is_label_path(path):
            raise PerturbationError(f"refusing a label-file input in a label-blind step: {path}")


def is_complete_ec(function_id: str) -> bool:
    value = function_id[3:] if function_id.startswith("EC:") else function_id
    parts = value.split(".")
    return len(parts) == 4 and all(p and p != "-" for p in parts)


def stable_rng(*parts) -> random.Random:
    return random.Random(hashlib.sha256("\t".join(str(p) for p in parts).encode()).hexdigest())


class Inputs:
    """Reads the joint WP4 inputs and records the sha256 of every file read."""

    def __init__(self, sequence_dir, cohort_layout, target_claims=None, peers_dir=None):
        self.sequence_dir = Path(sequence_dir).resolve()
        self.cohort_layout_dir = Path(cohort_layout).resolve()
        self.target_claims_dir = Path(target_claims).resolve() if target_claims else None
        self.peers_dir = Path(peers_dir).resolve() if peers_dir else None
        self.hashes = {}
        self.ortho_path = self.sequence_dir / "cohort_orthogroups.tsv"
        self.ortho_header, self.ortho_rows = self.read_orthogroups(self.ortho_path)
        metadata = json.loads(self.read_text(self.sequence_dir / "cohort_metadata.json"))
        require(isinstance(metadata.get("genomes"), dict), "joint metadata has no genomes object")
        self.metadata = metadata["genomes"]
        self.cohort_proteins = self.read_namespaced_fasta(self.cohort_layout_dir / "all_proteins.faa")
        self.cohort = set(self.cohort_proteins)
        require(self.cohort <= set(self.metadata), "cohort genomes missing from joint metadata")
        self.rows_by_genome = defaultdict(list)
        for row in self.ortho_rows:
            self.rows_by_genome[row["genome_id"]].append(row)
        cohort_flags = [row["genome_id"] in self.cohort for row in self.ortho_rows]
        require(cohort_flags == sorted(cohort_flags, reverse=True),
                "cohort orthogroup rows must form the prefix of the joint orthogroups")

    def record(self, path):
        path = Path(path).resolve()
        require(path.is_file(), f"missing input file: {path}")
        self.hashes[str(path)] = sha256_file(path)
        return path

    def read_text(self, path):
        return self.record(path).read_text()

    def read_bytes(self, path):
        return self.record(path).read_bytes()

    def read_orthogroups(self, path):
        data = self.read_bytes(path)
        require(data.endswith(b"\n"), f"orthogroups must end with a newline: {path}")
        lines = data.split(b"\n")[:-1]
        require(lines and lines[0].decode().split("\t") == ORTHO_FIELDS, f"unexpected orthogroup header: {path}")
        rows, seen = [], set()
        for number, raw in enumerate(lines[1:], 2):
            fields = raw.decode().split("\t")
            require(len(fields) == 5 and all(fields[:3]), f"invalid orthogroup row at {path}:{number}")
            key = (fields[0], fields[1])
            require(key not in seen, f"protein has two family rows: {key}")
            seen.add(key)
            rows.append({"genome_id": fields[0], "protein_id": fields[1], "fields": fields, "raw": raw + b"\n"})
        return lines[0] + b"\n", rows

    def read_namespaced_fasta(self, path):
        proteins = defaultdict(dict)
        for identifier, sequence in parse_fasta(self.read_text(path), path):
            require("__" in identifier, f"cohort FASTA id is not namespaced <genome>__<protein>: {identifier}")
            genome, protein = identifier.split("__", 1)
            require(protein not in proteins[genome], f"duplicate cohort protein {identifier}")
            proteins[genome][protein] = sequence
        return dict(proteins)

    def cohort_layout(self, genome):
        return read_layout(self.record(self.cohort_layout_dir / f"{genome}.layout.tsv"))

    def target_claims(self, genome):
        require(self.target_claims_dir is not None, "--target-claims is required")
        return read_claims(self.record(self.target_claims_dir / f"{genome}.claims.jsonl"))

    def peer_claims(self, accession):
        return read_claims(self.record(self.sequence_dir / "claims" / f"{accession}.claims.jsonl"))

    def peer_file(self, accession, name):
        require(self.peers_dir is not None, "--peers-dir is required")
        return self.record(self.peers_dir / accession / name)

    def peer_layout(self, accession):
        return read_layout(self.peer_file(accession, "layout.tsv"))

    def verify_unchanged(self):
        changed = [p for p, h in self.hashes.items() if sha256_file(p) != h]
        require(not changed, "input changed during preparation: " + ", ".join(changed))


def parse_fasta(text, path="FASTA"):
    entries, identifier, chunks = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if identifier is not None:
                entries.append((identifier, "".join(chunks)))
            fields = line[1:].split()
            require(fields, f"empty FASTA header in {path}")
            identifier, chunks = fields[0], []
        elif line.strip():
            require(identifier is not None, f"sequence before the first header in {path}")
            chunks.append(line.strip())
    if identifier is not None:
        entries.append((identifier, "".join(chunks)))
    return entries


def read_layout(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == LAYOUT_FIELDS, f"layout header must be {LAYOUT_FIELDS}: {path}")
        rows = []
        for number, row in enumerate(reader, 2):
            require(None not in row and all(row[f] for f in LAYOUT_FIELDS), f"malformed layout row {path}:{number}")
            start, end = int(row["start"]), int(row["end"])
            require(0 < start <= end, f"invalid CDS coordinates at {path}:{number}")
            rows.append({**row, "start": start, "end": end})
    ids = [r["protein_id"] for r in rows]
    require(len(ids) == len(set(ids)), f"duplicate layout protein in {path}")
    return rows


def read_claims(path):
    """``[(protein_id, raw line, parsed object)]``, refusing claims that carry product names."""
    claims = []
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            obj = json.loads(line)
            require(isinstance(obj, dict) and obj.get("protein_id"), f"claim without protein_id at {path}:{number}")
            require(not _product_keys(obj), f"claim carries a product-name field at {path}:{number}")
            claims.append((obj["protein_id"], line if line.endswith("\n") else line + "\n", obj))
    return claims


def _product_keys(value):
    if isinstance(value, dict):
        return {k for k in value if k.lower() in PRODUCT_KEYS} | set().union(*(_product_keys(v) for v in value.values()))
    if isinstance(value, list):
        return set().union(*(_product_keys(v) for v in value)) if value else set()
    return set()


def check_partition(genome, proteins, layout_rows, unplaced, contigs):
    placed = {r["protein_id"] for r in layout_rows}
    unplaced = set(unplaced)
    require(not placed & unplaced, f"{genome}: protein both placed and unplaced")
    require(placed | unplaced == set(proteins), f"{genome}: FASTA ids differ from layout + unplaced")
    for row in layout_rows:
        require(row["contig"] in contigs, f"{genome}: layout contig absent from metadata: {row['contig']}")
        length = contigs[row["contig"]].get("length")
        require(length is None or row["end"] <= length, f"{genome}: CDS beyond contig end: {row['protein_id']}")


# ------------------------------------------------------------------------------------------------ writer
class Output:
    def __init__(self, out):
        self.out = Path(out).resolve()
        require(not self.out.exists(), f"refusing to overwrite an existing output directory: {self.out}")
        self.out.mkdir(parents=True)

    def path(self, *parts):
        path = self.out.joinpath(*parts)
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def fasta(self, relative, proteins):
        with self.path(relative).open("w") as handle:
            for protein in sorted(proteins):
                require(protein and not any(c.isspace() for c in protein), f"invalid protein id {protein!r}")
                handle.write(f">{protein}\n{proteins[protein]}\n")

    def layout(self, relative, rows):
        with self.path(relative).open("w") as handle:
            handle.write("\t".join(LAYOUT_FIELDS) + "\n")
            for row in sorted(rows, key=lambda r: (r["contig"], r["start"], r["end"], r["protein_id"])):
                handle.write("\t".join(str(row[f]) for f in LAYOUT_FIELDS) + "\n")

    def lines(self, relative, lines):
        with self.path(relative).open("w") as handle:
            handle.writelines(lines)

    def genome(self, copy_id, proteins, layout_rows, unplaced, contigs, taxonomy, claim_lines):
        check_partition(copy_id, proteins, layout_rows, unplaced, contigs)
        bad = [p for p in (json.loads(line)["protein_id"] for line in claim_lines) if p not in proteins]
        require(not bad, f"{copy_id}: claim for a protein outside the proteome: {bad[:3]}")
        base = Path("genomes") / copy_id
        self.fasta(base / "protein.faa", proteins)
        self.layout(base / "layout.tsv", layout_rows)
        self.lines(base / "unplaced.txt", [p + "\n" for p in sorted(unplaced)])
        save_json(self.path(base / "contigs.json"), contigs)
        save_json(self.path(base / "taxonomy.json"), taxonomy)
        self.lines(Path("claims") / f"{copy_id}.claims.jsonl", claim_lines)
        return base / "protein.faa"

    def joint(self, prefix, inputs, population, extra_rows, extra_metadata, metadata_overrides=None):
        """Joint orthogroups (input rows of the population, unchanged bytes, then new rows) and metadata."""
        population = set(population)
        lines = [inputs.ortho_header]
        lines += [row["raw"] for row in inputs.ortho_rows if row["genome_id"] in population]
        for fields in extra_rows:
            lines.append(("\t".join(fields) + "\n").encode())
        self.path(Path(prefix) / "cohort_orthogroups.tsv").write_bytes(b"".join(lines))
        genomes = {g: inputs.metadata[g] for g in population}
        genomes.update(metadata_overrides or {})
        genomes.update(extra_metadata)
        save_json(self.path(Path(prefix) / "cohort_metadata.json"),
                  {"genomes": {g: genomes[g] for g in sorted(genomes)}})
        return len(lines) - 1

    def scan_no_product_names(self):
        for path in sorted(self.out.rglob("*")):
            if not path.is_file():
                continue
            if path.suffix in (".faa", ".fasta"):
                for line in path.read_text().splitlines():
                    require(not line.startswith(">") or (len(line) > 1 and not any(c.isspace() for c in line[1:])),
                            f"non-bare FASTA header in output {path}: {line[:80]!r}")
            elif path.suffix == ".jsonl":
                for line in path.read_text().splitlines():
                    require(not _product_keys(json.loads(line)), f"product-name field in output {path}")
            elif path.name == "layout.tsv":
                require(path.read_text().split("\n", 1)[0].split("\t") == LAYOUT_FIELDS,
                        f"layout header changed in output {path}")

    def receipt(self, mode, inputs, parameters, counts, label_inputs=None, extra=None):
        inputs.verify_unchanged()
        self.scan_no_product_names()
        script = Path(__file__).resolve()
        label_inputs = label_inputs or {}
        outputs = {p.relative_to(self.out).as_posix(): sha256_file(p)
                   for p in sorted(self.out.rglob("*")) if p.is_file()}
        receipt = {"schema": SCHEMA, "mode": mode, "created_utc": datetime.now(timezone.utc).isoformat(),
                   "script": {"path": str(script), "sha256": sha256_file(script)},
                   "inputs": {p: h for p, h in sorted(inputs.hashes.items()) if p not in label_inputs},
                   "label_inputs": label_inputs, "label_exception": LABEL_EXCEPTION if label_inputs else None,
                   "parameters": parameters, "counts": counts, "outputs": outputs, **(extra or {})}
        save_json(self.out / "receipt.json", receipt)
        return receipt


def write_variants(output, variants):
    """``variants.tsv``: the analysis arm variant names (Task Q conventions) and the genomes each seals."""
    lines = ["\t".join(VARIANT_FIELDS) + "\n"]
    for row in variants:
        lines.append("\t".join(str(row[f]) for f in VARIANT_FIELDS) + "\n")
    output.lines("variants.tsv", lines)


def gapsmith_inputs(output, rows):
    lines = ["genome_id\tproteome\tsha256\tsource_genome\tperturbation\n"]
    for genome, relative, source, kind in sorted(rows):
        lines.append(f"{genome}\t{relative.as_posix()}\t{sha256_file(output.out / relative)}\t{source}\t{kind}\n")
    output.lines("gapsmith_inputs.tsv", lines)


def copy_metadata(source, assembly=None, drop_plausible=True, **updates):
    record = copy.deepcopy(source)
    if drop_plausible:
        record.pop("plausibleReactions", None)
    if assembly:
        record.update(assembly)
    record.update(updates)
    return record


def target_genome(inputs, genome):
    require(genome in inputs.cohort, f"not a cohort genome: {genome}")
    meta = inputs.metadata[genome]
    proteins = inputs.cohort_proteins[genome]
    layout_rows = inputs.cohort_layout(genome)
    unplaced = list(meta.get("unplacedProteins", []))
    check_partition(genome, proteins, layout_rows, unplaced, meta["contigs"])
    claims = inputs.target_claims(genome)
    bad = sorted({p for p, _, _ in claims if p not in proteins})
    require(not bad, f"{genome}: claims name undeclared proteins: {bad[:3]}")
    return meta, proteins, layout_rows, unplaced, claims


def taxonomy_record(meta, source, kind, **extra):
    return {"clade": meta.get("clade"), "strainCluster": meta.get("strainCluster"),
            "assemblyStatus": meta.get("assemblyStatus"), "source_genome": source, "perturbation": kind, **extra}


# ----------------------------------------------------------------------------------------------- removal
def read_label_rows(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome_id", "protein_id", "function_id", "tier"}
        require(required <= set(reader.fieldnames or ()), f"labels file lacks {sorted(required)}: {path}")
        return [{k: row[k] for k in required} for row in reader]


def ec_value(function_id):
    value = str(function_id).strip()
    return value[3:] if value.startswith("EC:") else value


def select_removals(genome, proteins, claims, label_rows, seed):
    """Strict-EC-labelled proteins + an equal number of random unlabelled complete-EC-claimed proteins.

    Returns the labelled proteins, the random proteins, the pool size and ``{protein: [(EC, source)]}``: the
    strict EC labels of a labelled protein, the complete EC claims of a random one (its base-rate function).
    """
    strict = defaultdict(set)
    for r in label_rows:
        if (r["genome_id"] == genome and r["tier"] == "experimental" and r["function_id"].startswith("EC:")
                and is_complete_ec(r["function_id"])):
            strict[r["protein_id"]].add(ec_value(r["function_id"]))
    labelled = sorted(strict)
    missing = [p for p in labelled if p not in proteins]
    require(not missing, f"{genome}: labelled proteins absent from the proteome: {missing[:3]}")
    any_label = {r["protein_id"] for r in label_rows if r["genome_id"] == genome}
    ec_claimed = {p for p, _, obj in claims if obj.get("function_type") == "EC" and is_complete_ec(str(obj.get("function_id", "")))}
    pool = sorted(ec_claimed - any_label)
    require(len(pool) >= len(labelled), f"{genome}: only {len(pool)} unlabelled EC-claimed proteins for "
                                        f"{len(labelled)} labelled removals")
    chosen = sorted(random.Random(f"{seed}\t{genome}").sample(pool, len(labelled)))
    functions = {p: [(ec, "label") for ec in sorted(strict[p])] for p in labelled}
    for p in chosen:
        functions[p] = [(ec, "claim") for ec in sorted({ec_value(obj.get("function_id")) for q, _, obj in claims
                                                        if q == p and obj.get("function_type") == "EC"
                                                        and is_complete_ec(str(obj.get("function_id", "")))})]
    return labelled, chosen, len(pool), functions


def run_removal(sequence_dir, cohort_layout, target_claims, labels, genomes, out, seed=REMOVAL_SEED):
    refuse_label_paths([sequence_dir, cohort_layout, target_claims, out])
    require(genomes, "at least one --genome is required")
    output = Output(out)
    inputs = Inputs(sequence_dir, cohort_layout, target_claims)
    label_rows = read_label_rows(inputs.record(labels))
    label_inputs = {str(Path(labels).resolve()): inputs.hashes[str(Path(labels).resolve())]}
    extra_rows, extra_meta, search_rows, selection = [], {}, [], ["\t".join(REMOVAL_FIELDS) + "\n"]
    counts, variants = {}, []
    for genome in sorted(set(genomes)):
        meta, proteins, layout_rows, unplaced, claims = target_genome(inputs, genome)
        labelled, randoms, pool, functions = select_removals(genome, proteins, claims, label_rows, seed)
        copies = []
        category = {p: "labelled" for p in labelled}
        category.update({p: "random" for p in randoms})
        family_rows = inputs.rows_by_genome[genome]
        for index, protein in enumerate(sorted(category), 1):
            copy_id = f"{genome}-rm{index:03d}"
            require(copy_id not in inputs.metadata, f"copy id collision: {copy_id}")
            kept = {p: s for p, s in proteins.items() if p != protein}
            relative = output.genome(
                copy_id, kept, [r for r in layout_rows if r["protein_id"] != protein],
                [p for p in unplaced if p != protein], meta["contigs"],
                taxonomy_record(meta, genome, "removal"),
                [line for p, line, _ in claims if p != protein])
            extra_rows += [[copy_id, *r["fields"][1:]] for r in family_rows if r["protein_id"] != protein]
            extra_meta[copy_id] = copy_metadata(meta, unplacedProteins=sorted(p for p in unplaced if p != protein))
            search_rows.append((copy_id, relative, genome, "removal"))
            copies.append(copy_id)
            for ec, source in functions[protein]:
                selection.append("\t".join([copy_id, genome, protein, category[protein], f"EC:{ec}", source]) + "\n")
        variants.append({"variant": "removal", "target": genome, "sealed_genomes": sealed(genome, copies),
                         "path": "."})
        counts[genome] = {"labelled_removed": len(labelled), "random_removed": len(randoms),
                          "unlabelled_ec_claimed_pool": pool, "copies": len(category),
                          "proteins": len(proteins), "claims": len(claims)}
    rows = output.joint(".", inputs, inputs.metadata, extra_rows, extra_meta)
    output.lines("removal_selection.tsv", selection)
    write_variants(output, variants)
    gapsmith_inputs(output, search_rows)
    return output.receipt("removal", inputs, {
        "seed": seed, "genomes": sorted(set(genomes)),
        "selection": ("every protein with a strict experimental complete-EC label; plus, per genome, the same number "
                      "drawn with random.Random('<seed>\\t<genome>').sample from the sorted proteins that have a "
                      "complete EC claim and no labels.tsv row of any tier"),
        "copy_rule": "one copy per removed protein; clade and strainCluster of the target; plausibleReactions dropped "
                     "(fresh gapsmith search on the reduced proteome)"},
        {"genomes": counts, "copies": len(search_rows), "joint_orthogroup_rows": rows},
        label_inputs=label_inputs, extra={"variants": variants, "removals_manifest": "removal_selection.tsv"})


# ---------------------------------------------------------------------------------------------- fragment
def draw_seed(median, sigma, min_bp, seed):
    return f"wp4-fragment|median={median}|sigma={sigma}|min={min_bp}|seed={seed}"


def fragment_draws(median, sigma=FRAGMENT_SIGMA, min_bp=FRAGMENT_MIN_BP, seed=1):
    """The fixed i.i.d. log-normal length stream for one (median, seed); identical for every genome."""
    rng = random.Random(draw_seed(median, sigma, min_bp, seed))
    mu = math.log(median)
    while True:
        yield max(min_bp, int(round(rng.lognormvariate(mu, sigma))))


def cut_contigs(contigs, draws):
    """Cut each contig (sorted by id) with consecutive draws; the last piece of a contig ends at its end."""
    pieces, used = {}, []
    for contig in sorted(contigs):
        length = contigs[contig]["length"]
        require(type(length) is int and length > 0, f"contig without a positive length: {contig}")
        position, spans = 0, []
        while position < length:
            draw = next(draws)
            used.append(draw)
            end = min(length, position + draw)
            spans.append((position + 1, end))
            position = end
        pieces[contig] = spans
    return pieces, used


def fragment_layout(layout_rows, pieces):
    """Keep CDSs wholly inside one piece (re-based coordinates); CDSs crossing a cut are lost."""
    kept, lost = [], []
    starts = {contig: [a for a, _ in spans] for contig, spans in pieces.items()}
    for row in layout_rows:
        spans = pieces[row["contig"]]
        index = bisect.bisect_right(starts[row["contig"]], row["start"]) - 1
        a, b = spans[index]
        if row["end"] <= b:
            kept.append({**row, "contig": f"{row['contig']}__frag_{a}_{b}",
                         "start": row["start"] - a + 1, "end": row["end"] - a + 1})
        else:
            lost.append(row["protein_id"])
    return kept, sorted(lost)


def run_fragment(sequence_dir, cohort_layout, target_claims, genomes, out, medians=FRAGMENT_MEDIANS,
                 seeds=FRAGMENT_SEEDS, sigma=FRAGMENT_SIGMA, min_bp=FRAGMENT_MIN_BP):
    refuse_label_paths([sequence_dir, cohort_layout, target_claims, out])
    require(genomes, "at least one --genome is required")
    require(all(type(m) is int and m >= min_bp for m in medians), "medians must be integers >= min_bp")
    require(sigma > 0 and min_bp > 0, "sigma and min_bp must be positive")
    output = Output(out)
    inputs = Inputs(sequence_dir, cohort_layout, target_claims)
    extra_rows, extra_meta, search_rows, counts, draws_record = [], {}, [], {}, {}
    variants = []
    for genome in sorted(set(genomes)):
        meta, proteins, layout_rows, unplaced, claims = target_genome(inputs, genome)
        family_rows = inputs.rows_by_genome[genome]
        counts[genome] = {}
        for median in medians:
            for seed in seeds:
                tag = fragment_tier(median)
                copy_id = f"{genome}-frag{tag}-s{seed}"
                require(copy_id not in inputs.metadata, f"copy id collision: {copy_id}")
                pieces, used = cut_contigs(meta["contigs"], fragment_draws(median, sigma, min_bp, seed))
                key = f"median={median};seed={seed}"
                record = draws_record.setdefault(key, {"median_bp": median, "sigma": sigma, "min_bp": min_bp,
                                                       "seed": seed, "rng_seed": draw_seed(median, sigma, min_bp, seed),
                                                       "draws": [], "draws_used": {}})
                require(record["draws"][:len(used)] == used[:len(record["draws"])], "draw stream differs by genome")
                if len(used) > len(record["draws"]):
                    record["draws"] = used
                record["draws_used"][genome] = len(used)
                kept_rows, lost = fragment_layout(layout_rows, pieces)
                lost_set = set(lost)
                contigs = {f"{c}__frag_{a}_{b}": {"length": b - a + 1, "closed": False}
                           for c, spans in pieces.items() for a, b in spans}
                kept = {p: s for p, s in proteins.items() if p not in lost_set}
                kept_claims = [line for p, line, _ in claims if p not in lost_set]
                relative = output.genome(copy_id, kept, kept_rows, unplaced, contigs,
                                         taxonomy_record(meta, genome, "fragment", assemblyStatus="DRAFT",
                                                         fragment=key), kept_claims)
                extra_rows += [[copy_id, *r["fields"][1:]] for r in family_rows if r["protein_id"] not in lost_set]
                extra_meta[copy_id] = copy_metadata(meta, {"assemblyStatus": "DRAFT", "contigs": contigs},
                                                    unplacedProteins=sorted(unplaced))
                search_rows.append((copy_id, relative, genome, "fragment"))
                variants.append({"variant": f"fragment:{tag}:{seed}", "target": genome,
                                 "sealed_genomes": sealed(genome, [copy_id]), "path": "."})
                counts[genome][copy_id] = {"median_bp": median, "seed": seed, "pieces": len(contigs),
                                           "draws_used": len(used), "boundary_cds_removed": len(lost),
                                           "proteins_kept": len(kept), "claims_kept": len(kept_claims),
                                           "claims_removed": len(claims) - len(kept_claims),
                                           "removed_proteins": lost}
    for record in draws_record.values():
        record["draws_sha256"] = hashlib.sha256(json.dumps(record["draws"]).encode()).hexdigest()
    save_json(output.path("fragment_draws.json"), draws_record)
    rows = output.joint(".", inputs, inputs.metadata, extra_rows, extra_meta)
    write_variants(output, variants)
    gapsmith_inputs(output, search_rows)
    return output.receipt("fragment", inputs, {
        "medians_bp": list(medians), "sigma": sigma, "min_bp": min_bp, "seeds": list(seeds),
        "genomes": sorted(set(genomes)),
        "draw_rule": ("one random.Random(rng_seed) stream per (median, seed); lengths = max(min_bp, round("
                      "lognormvariate(ln median, sigma))); contigs in sorted order consume consecutive draws, the "
                      "last piece of each contig is truncated at the contig end; the same stream for every genome"),
        "boundary_rule": "CDS crossing a cut is removed from proteome, claims, layout and family rows",
        "copy_rule": "DRAFT, pieces closed=false, unplaced proteins retained; clade and strainCluster of the target; "
                     "plausibleReactions dropped (fresh gapsmith search)"},
        {"genomes": counts, "copies": len(search_rows), "joint_orthogroup_rows": rows},
        extra={"variants": variants})


# ----------------------------------------------------------------------------------------------- shuffle
def window_strata(layout_rows, contigs, half_width=HALF_WIDTH_BP):
    """Engine bucket per placed protein: (contig, window count, window state, closed key).

    Mirrors ``CohortInference.windowsOf``: midpoint = (start + end) / 2 truncated, window = loci on the same
    contig with midpoint within +-halfWidthBp (centre included), state FULL/TRUNCATED/UNKNOWN as
    ``windowState`` and closed key closed/open/unk as ``closedKey``.
    """
    midpoints = defaultdict(list)
    for row in layout_rows:
        midpoints[row["contig"]].append((row["start"] + row["end"]) // 2)
    for values in midpoints.values():
        values.sort()
    strata = {}
    for row in layout_rows:
        contig, m = row["contig"], (row["start"] + row["end"]) // 2
        values = midpoints[contig]
        count = bisect.bisect_right(values, m + half_width) - bisect.bisect_left(values, m - half_width)
        info = contigs.get(contig) or {}
        closed, length = info.get("closed"), info.get("length")
        left, right = m - half_width < 1, length is not None and m + half_width > length
        if closed is True:
            state = "FULL"
        elif closed is False and (left or right):
            state = "TRUNCATED"
        elif left or right:
            state = "UNKNOWN"
        else:
            state = "FULL" if length is not None else "UNKNOWN"
        closed_key = "unk" if closed is None else ("closed" if closed else "open")
        strata[row["protein_id"]] = (contig, count, state, closed_key)
    return strata


def shuffle_families(rows, strata, seed, genome):
    """Permute (orthogroup, subfamily, weight) payloads within each stratum; ``{protein: payload}``."""
    by_stratum = defaultdict(list)
    for row in rows:
        stratum = strata.get(row["protein_id"])
        if stratum is not None:
            by_stratum[stratum].append(row)
    assigned = {}
    for stratum in sorted(by_stratum):
        members = sorted(by_stratum[stratum], key=lambda r: r["protein_id"])
        payloads = [tuple(r["fields"][2:]) for r in members]
        stable_rng("wp4-s1", seed, genome, *stratum).shuffle(payloads)
        for row, payload in zip(members, payloads):
            assigned[row["protein_id"]] = payload
    return assigned, len(by_stratum)


def read_pathway_groups(path):
    with Path(path).open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        require(next(reader, None) == GROUP_FIELDS, f"pathway-groups header must be {GROUP_FIELDS}: {path}")
        edges, seen = [], set()
        for number, row in enumerate(reader, 2):
            require(len(row) == 4 and all(row), f"invalid pathway-groups row at {path}:{number}")
            key = (row[0], row[1], row[2])
            require(key not in seen, f"duplicate pathway-groups row at {path}:{number}")
            seen.add(key)
            edges.append([(row[0], row[1]), row[2], row[3]])
    require(edges, "pathway-groups file is empty")
    return edges


def swap_pathway_groups(edges, seed, swap_factor=SWAP_FACTOR):
    """Degree-preserving reaction swaps between groups; returns new edges and accepted swap count."""
    edges = [list(e) for e in edges]
    present = {(g, r) for g, r, _ in edges}
    rng = random.Random(f"wp4-s2|seed={seed}")
    accepted, n = 0, len(edges)
    for _ in range(swap_factor * n if n > 1 else 0):
        i, j = rng.randrange(n), rng.randrange(n)
        (g1, r1, _), (g2, r2, _) = edges[i], edges[j]
        if g1 == g2 or r1 == r2 or (g1, r2) in present or (g2, r1) in present:
            continue
        present -= {(g1, r1), (g2, r2)}
        present |= {(g1, r2), (g2, r1)}
        edges[i][1], edges[j][1] = r2, r1
        accepted += 1
    return edges, accepted


def run_shuffle(sequence_dir, cohort_layout, peers_dir, pathway_groups, out, s1_seeds=S1_SEEDS,
                s2_seeds=S2_SEEDS, half_width=HALF_WIDTH_BP, peers=None, swap_factor=SWAP_FACTOR):
    refuse_label_paths([sequence_dir, cohort_layout, peers_dir, pathway_groups, out])
    require(type(half_width) is int and half_width > 0, "halfWidthBp must be a positive integer")
    output = Output(out)
    inputs = Inputs(sequence_dir, cohort_layout, peers_dir=peers_dir)
    available = sorted(g for g in inputs.rows_by_genome if g not in inputs.cohort)
    selected = sorted(set(peers)) if peers else available
    require(not set(selected) - set(available), "requested peers lack family rows: "
            + ", ".join(sorted(set(selected) - set(available))))
    strata_by_genome, stratum_counts = {}, {}
    for genome in selected:
        layout_rows = inputs.peer_layout(genome)
        contigs = inputs.metadata[genome]["contigs"]
        proteins = {r["protein_id"] for r in inputs.rows_by_genome[genome]}
        # Strata are computed on the full layout (the engine's window counts include proteins without a family
        # row); only proteins with a family row exchange payloads.
        strata_by_genome[genome] = window_strata(layout_rows, contigs, half_width)
        stratum_counts[genome] = {"strata": len(set(strata_by_genome[genome].values())),
                                  "placed_without_family_row": len({r["protein_id"] for r in layout_rows} - proteins)}
    counts = {"S1": {}, "S2": {}}
    variants = [{"variant": f"S1:{seed}", "target": "*", "sealed_genomes": "",
                 "path": f"S1/seed{seed}/cohort_orthogroups.tsv"} for seed in s1_seeds]
    cohort_prefix = b"".join([inputs.ortho_header] + [r["raw"] for r in inputs.ortho_rows
                                                       if r["genome_id"] in inputs.cohort])
    for seed in s1_seeds:
        assigned, moved, strata_used = {}, Counter(), {}
        for genome in selected:
            payloads, used = shuffle_families(inputs.rows_by_genome[genome], strata_by_genome[genome], seed, genome)
            strata_used[genome] = used
            for row in inputs.rows_by_genome[genome]:
                payload = payloads.get(row["protein_id"])
                if payload is not None and payload != tuple(row["fields"][2:]):
                    assigned[(genome, row["protein_id"])] = payload
                    moved[genome] += 1
        lines = [inputs.ortho_header]
        for row in inputs.ortho_rows:
            payload = assigned.get((row["genome_id"], row["protein_id"]))
            lines.append(row["raw"] if payload is None else
                         ("\t".join(row["fields"][:2] + list(payload)) + "\n").encode())
        data = b"".join(lines)
        require(data.startswith(cohort_prefix), "S1 changed a cohort (target) family row")
        output.path("S1", f"seed{seed}", "cohort_orthogroups.tsv").write_bytes(data)
        counts["S1"][f"seed{seed}"] = {"proteins_moved": dict(sorted(moved.items())),
                                       "strata": strata_used, "rows": len(lines) - 1}
    if s2_seeds:
        require(pathway_groups is not None, "--pathway-groups is required for S2")
        edges = read_pathway_groups(inputs.record(pathway_groups))
        group_degree = Counter(g for g, _, _ in edges)
        reaction_degree = Counter(r for _, r, _ in edges)
        for seed in s2_seeds:
            swapped, accepted = swap_pathway_groups(edges, seed, swap_factor)
            require(Counter(g for g, _, _ in swapped) == group_degree, "S2 changed a group size")
            require(Counter(r for _, r, _ in swapped) == reaction_degree, "S2 changed a reaction group count")
            require(len({(g, r) for g, r, _ in swapped}) == len(swapped), "S2 produced a duplicate edge")
            rows = sorted((g[0], g[1], r, s) for g, r, s in swapped)
            output.lines(Path("S2") / f"seed{seed}" / "pathway-groups.tsv",
                         ["\t".join(GROUP_FIELDS) + "\n"] + ["\t".join(row) + "\n" for row in rows])
            original = {(g, r) for g, r, _ in edges}
            variants.append({"variant": f"S2:{seed}", "target": "*", "sealed_genomes": "",
                             "path": f"S2/seed{seed}/pathway-groups.tsv"})
            counts["S2"][f"seed{seed}"] = {"edges": len(swapped), "accepted_swaps": accepted,
                                           "edges_changed": sum((g, r) not in original for g, r, _ in swapped)}
    write_variants(output, variants)
    return output.receipt("shuffle", inputs, {
        "s1_seeds": list(s1_seeds), "s2_seeds": list(s2_seeds), "half_width_bp": half_width,
        "swap_factor": swap_factor, "peers": selected,
        "s1_rule": "peer genomes only; payload (orthogroup_id, subfamily_id, weight) permuted within the stratum "
                   "(contig, window gene count, window state, closed key); unplaced proteins and cohort rows "
                   "unchanged; random.Random(sha256('wp4-s1', seed, genome, stratum))",
        "s2_rule": "edge swaps (group_a, r1), (group_b, r2) -> (group_a, r2), (group_b, r1) rejected when an edge "
                   "would repeat; swap_factor x edges attempts; random.Random('wp4-s2|seed=<seed>')"},
        {**counts, "peer_strata": stratum_counts}, extra={"variants": variants})


# --------------------------------------------------------------------------------------------- duplicate
def read_peer_rows(path, target):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"target_accession", "peer_accession", "role"}
        require(required <= set(reader.fieldnames or ()), f"peers.tsv lacks {sorted(required)}")
        return [row for row in reader if row["target_accession"] == target]


def peer_copy(inputs, accession):
    """Exact files of an acquired peer genome."""
    proteins = dict(parse_fasta(inputs.peer_file(accession, "protein.faa").read_text(), accession))
    layout_rows = inputs.peer_layout(accession)
    unplaced = [l.strip() for l in inputs.peer_file(accession, "unplaced.txt").read_text().splitlines() if l.strip()]
    contigs = json.loads(inputs.peer_file(accession, "contigs.json").read_text())
    claims = inputs.peer_claims(accession)
    return proteins, layout_rows, unplaced, contigs, [line for _, line, _ in claims]


def run_duplicate(sequence_dir, cohort_layout, target_claims, peers_dir, peers_tsv, targets, out):
    refuse_label_paths([sequence_dir, cohort_layout, target_claims, peers_dir, peers_tsv, out])
    require(targets, "at least one --target is required")
    output = Output(out)
    inputs = Inputs(sequence_dir, cohort_layout, target_claims, peers_dir)
    inputs.record(peers_tsv)
    existing_clusters = {m.get("strainCluster") for m in inputs.metadata.values()}
    counts, variants = {}, []
    for target in sorted(set(targets)):
        require(target in inputs.cohort, f"duplicate target is not a cohort genome: {target}")
        meta_t = inputs.metadata[target]
        rows = read_peer_rows(peers_tsv, target)
        block = [r["peer_accession"] for r in rows]
        require(block, f"{target}: no peers.tsv rows")
        missing = [a for a in block if a not in inputs.metadata]
        require(not missing, f"{target}: peers absent from joint metadata: {missing}")
        tiers = {}
        for row in rows:
            accession = row["peer_accession"]
            if row["role"] == "near_duplicate":
                tiers[accession] = "ND"
            elif inputs.metadata[accession].get("clade") == meta_t.get("clade"):
                tiers[accession] = "C"
            else:
                tiers[accession] = "F"
        close = [a for a in block if tiers[a] == "C"]
        near = [a for a in block if tiers[a] == "ND"]
        require(len(near) <= 1, f"{target}: more than one near-duplicate")
        base = list(inputs.cohort) + block
        summary = {"tiers": tiers, "top_c_peer": close[0] if close else None,
                   "near_duplicate": near[0] if near else None, "variants": {}}

        def new_cluster(label):
            value = f"wp4dup_species_{label}"
            require(value not in existing_clusters, f"strain cluster collision: {value}")
            return value

        def write_peer_copy(variant, source, copy_id, strain):
            proteins, layout_rows, unplaced, contigs, claim_lines = peer_copy(inputs, source)
            meta = copy_metadata(inputs.metadata[source], drop_plausible=False, strainCluster=strain)
            output.genome(Path(target) / variant / copy_id, proteins, layout_rows, unplaced, contigs,
                          taxonomy_record(meta, source, "duplicate", variant=duplicate_variant_name(variant)), claim_lines)
            family = [[copy_id, *r["fields"][1:]] for r in inputs.rows_by_genome[source]]
            return family, {copy_id: meta}

        for variant in DUPLICATE_VARIANTS:
            prefix = Path(target) / variant
            population, extra_rows, extra_meta, overrides, note = list(base), [], {}, {}, None
            top, nd = summary["top_c_peer"], summary["near_duplicate"]
            if variant.startswith(("i_", "ii_")) and top is None:
                note = "skipped: no congeneric (C) peer"
            elif variant.startswith(("iii_", "v_")) and nd is None:
                note = "skipped: no near-duplicate peer"
            elif variant == "i_top_c_copy_same_cluster":
                extra_rows, extra_meta = write_peer_copy(variant, top, f"{top}-dup1",
                                                         inputs.metadata[top]["strainCluster"])
            elif variant == "ii_top_c_copy_new_species":
                extra_rows, extra_meta = write_peer_copy(variant, top, f"{top}-dup2", new_cluster(f"{top}-dup2"))
            elif variant == "iii_nd_own_cluster":
                overrides = {nd: copy_metadata(inputs.metadata[nd], drop_plausible=False,
                                               strainCluster=new_cluster(f"{nd}-own"))}
            elif variant == "iv_target_copy_congeneric":
                meta, proteins, layout_rows, unplaced, claims = target_genome(inputs, target)
                copy_id = f"{target}-dup4"
                record = copy_metadata(meta, drop_plausible=False, strainCluster=new_cluster(copy_id))
                output.genome(prefix / copy_id, proteins, layout_rows, unplaced, meta["contigs"],
                              taxonomy_record(record, target, "duplicate", variant=duplicate_variant_name(variant)),
                              [line for _, line, _ in claims])
                extra_rows = [[copy_id, *r["fields"][1:]] for r in inputs.rows_by_genome[target]]
                extra_meta = {copy_id: record}
            elif variant == "v_without_nd":
                population = [g for g in base if g != nd]
            if note:
                summary["variants"][variant] = {"status": note}
                continue
            rows_written = output.joint(prefix, inputs, population, extra_rows, extra_meta, overrides)
            genomes = sorted(set(population) | set(extra_meta))
            output.lines(prefix / "population.txt", [g + "\n" for g in genomes])
            target_copies = [c for c in extra_meta if c.startswith(f"{target}-")]
            variants.append({"variant": duplicate_variant_name(variant), "target": target,
                             "sealed_genomes": sealed(target, target_copies), "path": prefix.as_posix()})
            summary["variants"][variant] = {"status": "written", "variant": duplicate_variant_name(variant),
                                            "genomes": len(genomes),
                                            "copies": sorted(extra_meta), "relabelled": sorted(overrides),
                                            "orthogroup_rows": rows_written}
        counts[target] = summary
    write_variants(output, variants)
    return output.receipt("duplicate", inputs, {
        "targets": sorted(set(targets)), "variants": list(DUPLICATE_VARIANTS),
        "tier_rule": "ND = near_duplicate role; C = close_peer with the target's clade (genus); F = other close_peer",
        "top_c_rule": "first C row of the target in peers.tsv order (selection ranking)",
        "copy_rule": "exact copies keep protein ids, layout, unplaced, contigs, claims and family rows; only the "
                     "genome id (and the strain cluster where the variant says so) changes",
        "population": "18 cohort genomes + the family block's peers (+ copy, - near-duplicate for v)"},
        {"targets": counts}, extra={"variants": variants})


# ------------------------------------------------------------------------------------------------ verify
def verify(out):
    out = Path(out).resolve()
    receipt = json.loads((out / "receipt.json").read_text())
    require(receipt.get("schema") == SCHEMA, "not a WP4 perturbation receipt")
    failures = [p for p, h in {**receipt["inputs"], **receipt["label_inputs"]}.items()
                if not Path(p).is_file() or sha256_file(p) != h]
    actual = {p.relative_to(out).as_posix(): sha256_file(p) for p in sorted(out.rglob("*")) if p.is_file()}
    actual.pop("receipt.json", None)
    if actual != receipt["outputs"]:
        failures.append("outputs differ from the receipt")
    return {"passed": not failures, "failures": failures, "mode": receipt["mode"]}


# --------------------------------------------------------------------------------------------------- CLI
def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    steps = parser.add_subparsers(dest="step", required=True)

    def common(p, claims=True, peers=False):
        p.add_argument("--sequence-dir", type=Path, required=True,
                       help="WP4 sequence dir (cohort_orthogroups.tsv, cohort_metadata.json, claims/)")
        p.add_argument("--cohort-layout", type=Path, required=True,
                       help="cohort layout dir (<genome>.layout.tsv, all_proteins.faa)")
        if claims:
            p.add_argument("--target-claims", type=Path, required=True,
                           help="10406 diamond_homology claim dir (<genome>.claims.jsonl)")
        if peers:
            p.add_argument("--peers-dir", type=Path, required=True)
        p.add_argument("--out", type=Path, required=True)

    p = steps.add_parser("removal")
    common(p)
    p.add_argument("--labels", type=Path, required=True, help="labels.tsv (read only to choose removals)")
    p.add_argument("--genome", action="append", required=True)
    p.add_argument("--seed", type=int, default=REMOVAL_SEED)
    p = steps.add_parser("fragment")
    common(p)
    p.add_argument("--genome", action="append", required=True)
    p.add_argument("--medians", type=int, nargs="+", default=list(FRAGMENT_MEDIANS))
    p.add_argument("--seeds", type=int, nargs="+", default=list(FRAGMENT_SEEDS))
    p.add_argument("--sigma", type=float, default=FRAGMENT_SIGMA)
    p.add_argument("--min-bp", type=int, default=FRAGMENT_MIN_BP)
    p = steps.add_parser("shuffle")
    common(p, claims=False, peers=True)
    p.add_argument("--pathway-groups", type=Path)
    p.add_argument("--s1-seeds", type=int, nargs="*", default=list(S1_SEEDS))
    p.add_argument("--s2-seeds", type=int, nargs="*", default=list(S2_SEEDS))
    p.add_argument("--half-width-bp", type=int, default=HALF_WIDTH_BP)
    p.add_argument("--swap-factor", type=int, default=SWAP_FACTOR)
    p.add_argument("--peer", action="append")
    p = steps.add_parser("duplicate")
    common(p, peers=True)
    p.add_argument("--peers-tsv", type=Path, required=True)
    p.add_argument("--target", action="append", required=True)
    p = steps.add_parser("verify")
    p.add_argument("--out", type=Path, required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.step == "removal":
        result = run_removal(args.sequence_dir, args.cohort_layout, args.target_claims, args.labels,
                             args.genome, args.out, args.seed)
    elif args.step == "fragment":
        result = run_fragment(args.sequence_dir, args.cohort_layout, args.target_claims, args.genome, args.out,
                              tuple(args.medians), tuple(args.seeds), args.sigma, args.min_bp)
    elif args.step == "shuffle":
        result = run_shuffle(args.sequence_dir, args.cohort_layout, args.peers_dir, args.pathway_groups, args.out,
                             tuple(args.s1_seeds), tuple(args.s2_seeds), args.half_width_bp, args.peer,
                             args.swap_factor)
    elif args.step == "duplicate":
        result = run_duplicate(args.sequence_dir, args.cohort_layout, args.target_claims, args.peers_dir,
                               args.peers_tsv, args.target, args.out)
    else:
        result = verify(args.out)
        print(json.dumps(result, indent=2))
        return 0 if result["passed"] else 1
    print(json.dumps({"mode": result["mode"], "counts": result["counts"]}, indent=1, sort_keys=True)[:4000])
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except PerturbationError as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
