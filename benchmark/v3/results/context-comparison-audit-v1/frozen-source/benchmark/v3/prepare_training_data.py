#!/usr/bin/env python3
"""Training-table preparation and split audit for the learned context model.

Two subcommands, both read-only on their inputs:

``audit-splits``  validates a genome split table against the cohort registry
                  (``registry/training_cohort_registry.json``) and, optionally, a protein
                  list / homology-cluster map. Exit status 2 on any violation. The report names
                  the axes that were CHECKED and the ones that were NOT; it never turns the
                  registry's role labels into an "untouched test" claim.
``build-table``   joins production feature rows (``gspa integrate --feature-export-dir``)
                  with independently curated labels and the audited split table, and writes
                  the ``scripts/train_context_model.py`` input plus a manifest of exact counts.

Context proteins vs labelled examples
-------------------------------------
Whole-genome context (every protein of a genome, used only to compute features) CANNOT be made
cluster-disjoint across genomes: conserved core families occur in every genome. The protocol
therefore separates two populations:

* *context proteins* - unlabelled genome content. May share homologs across splits; sharing is
  counted and reported, never a violation.
* *labelled examples* - rows of the training table (positives, declared negatives) that the model
  is fitted on or evaluated against. With ``--protein-clusters`` the builder removes every
  labelled example of a lower-precedence split whose family also occurs among the labelled
  examples of a higher-precedence one (test > calibration > train: test targets are never
  dropped). Every removed example is written to ``masked_examples.tsv`` with its cluster and the
  split it clashes with, counted per split in the manifest, and a split left without positives
  stops the build. Hard test cases (family seen in train) are retained and counted.

Label policy (nothing is invented)
----------------------------------
A candidate row is labelled 1 only when an allowed-evidence truth assertion names the same
protein and function (GO optionally propagated over ``is_a``/``part_of``; obsolete GO terms are
dropped and counted). Explicit ``label=0`` or ``qualifier=NOT`` assertions require an
``explicit_negative_truth`` policy source and a row ``reference``; their polarity is preserved,
GO negatives stay exact, and contradictory assertions fail. Absence of an assertion is NOT a
negative. An unconfirmed row becomes 0 only when its dataset declares
``negative_policy: closed_world_declared`` for that function type with a written
justification, only for the type's default role, and only for proteins that carry at least one
truth assertion of that type (labelled proteins); every such row is stamped
``label_basis: closed_world_unconfirmed``. Everything else is excluded and counted
(``no_truth_not_negative``).

Roles: an EC assertion labels only CATALYTIC candidates and a GO assertion only UNSPECIFIED ones.
A COMPLEX_SUBUNIT (or any other) role is labelled only from role truth that says so: a labels-file
``role`` column whose role the dataset's policy declares under ``role_truth`` with its source, or,
when the policy opts in via ``complex_membership_from_cc`` (with a source) and a GO OBO is supplied,
an experimental cellular-component assertion to a protein-containing-complex term (GO:0032991
descendant) for a candidate naming that same complex term. An unrelated MF/BP label never labels a
subunit role, and no role negatives are ever inferred. ``locus_group_id`` is only set from an
explicit locus-group file whose group declares complete truth. Every exported candidate of the
function/role/genome must be declared, independently labelled and retained after family masking.
Incomplete groups retain marginal rows but lose locus supervision; unconfirmed closed-world rows
cannot certify a null outcome. Complete groups carry membership and source certificates.

Evidence: ``allowed_evidence`` is ``experimental`` (strict, the default: EXP/IDA/IPI/IMP/IGI/IEP/H*
and ECO:0000269), ``cafa_eligible`` (adds the curated-inferential TAS/IC/ECO:0000303/305, which are
not experimental) or ``all``/an explicit list. See ``evidence_codes.py``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import evidence_codes as EV  # noqa: E402
from v3_metrics import GO_ROOTS, load_obo  # noqa: E402

REGISTRY_DEFAULT = HERE / "registry" / "training_cohort_registry.json"
SPLITS = ("train", "calibration", "test")
REGISTRY_SCHEMA = "gspa.training-cohort-registry/1"
MANIFEST_SCHEMA = "gspa.training-table-manifest/1"
EXPORT_MANIFEST_SCHEMA = "gspa.context-feature-export/1"
EXPORT_ROW_SCHEMA = "gspa.context-feature-row/1"
EXPORT_ROWS_FILE, EXPORT_CONTRACT_FILE, EXPORT_MANIFEST_FILE = (
    "candidate_features.jsonl", "feature_contract.json", "export_manifest.json")
EXPORT_CONTEXT_KEYS = ("sequence_evidence", "predictor_status", "operon", "layout", "gap", "neighborhood",
                       "cohort", "assembly", "assembly_completeness", "fba", "fba_suppressed")
CANDIDATE_ROLES = ("CATALYTIC", "ACCESSORY", "TRANSPORTER", "REGULATOR", "COMPLEX_SUBUNIT", "UNSPECIFIED")
NEGATIVE_POLICIES = ("none", "closed_world_declared")
# One default (positive-labelled) role per function type. Other roles need their own role truth.
DEFAULT_ROLE_MAP = {"EC": ("CATALYTIC",), "GO": ("UNSPECIFIED",)}
COMPLEX_ROLE = "COMPLEX_SUBUNIT"
GO_COMPLEX, GO_CC_ROOT = "GO:0032991", "GO:0005575"
TAXON_KEYS = ("species", "genus", "family")
DISJOINT_KEYS = ("genome_id", "group") + TAXON_KEYS
# Same split disjointness the trainer enforces. Taxonomy keys are opt-in (and then mandatory).
DEFAULT_DISJOINT = ("genome_id", "group")
# Precedence for family masking: a lower-precedence split loses examples that clash with these.
MASK_PRECEDENCE = {"train": ("calibration", "test"), "calibration": ("test",), "test": ()}


class PreparationError(ValueError):
    """Inputs violate the training-data protocol; nothing may be written."""


# ------------------------------------------------------------------------------------ helpers

def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path, required):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in required if c not in (reader.fieldnames or [])]
        if missing:
            raise PreparationError(f"{path}: missing column(s) {missing}")
        return list(reader)


def load_registry(path=REGISTRY_DEFAULT) -> dict:
    reg = json.loads(Path(path).read_text())
    if reg.get("schema") != REGISTRY_SCHEMA:
        raise PreparationError(f"{path}: schema must be {REGISTRY_SCHEMA}")
    out = {}
    for d in reg["datasets"]:
        for key in ("id", "role", "allowed_splits"):
            if key not in d:
                raise PreparationError(f"registry dataset without '{key}': {d}")
        bad = [s for s in d["allowed_splits"] if s not in SPLITS]
        if bad:
            raise PreparationError(f"registry dataset {d['id']}: unknown split(s) {bad}")
        out[d["id"]] = d
    return out


def function_id(kind: str, value: str) -> str:
    v = value.strip().replace("_", ":")
    if kind == "EC" and not v.upper().startswith("EC:"):
        v = "EC:" + v
    return v


def contract_fingerprint(contract: dict) -> str:
    md = hashlib.sha256(contract["id"].encode())
    for f in contract["features"]:
        md.update(("\n" + f["name"] + "\t" + f["group"]).encode())
    return md.hexdigest()


# ------------------------------------------------------------------------------------ split audit

_SPLIT_COLUMNS = ("genome_id", "split", "group", "source")


def read_split_table(path) -> list[dict]:
    rows = read_tsv(path, _SPLIT_COLUMNS)
    for i, r in enumerate(rows, 2):
        for c in _SPLIT_COLUMNS:
            if not (r[c] or "").strip():
                raise PreparationError(f"{path} line {i}: empty '{c}'")
            r[c] = r[c].strip()
        for c in TAXON_KEYS:
            r[c] = (r.get(c) or "").strip()
    return rows


def _normalise_proteins(proteins):
    """``(genome, protein[, kind])`` -> ``(genome, protein, kind)``; kind is ``target`` or ``context``."""
    out = []
    for t in proteins:
        genome, protein, kind = (tuple(t) + ("target",))[:3]
        kind = (kind or "target").strip() or "target"
        if kind not in ("target", "context"):
            raise PreparationError(f"protein {protein!r}: kind {kind!r} must be 'target' or 'context'")
        out.append((genome, protein, kind))
    return out


def audit_splits(rows, registry, *, disjoint=DEFAULT_DISJOINT, proteins=None, protein_clusters=None,
                 min_test_genomes=5, cluster_policy="refuse") -> dict:
    """Return ``{violations, warnings, counts, claims, ...}``; never raises on a finding.

    ``proteins``: ``(genome_id, protein_id[, kind])``. ``protein_id`` must be a cross-genome identity
    (accession / sequence hash), not a per-genome sequence id. ``kind='target'`` proteins are the
    labelled/evaluated examples; ``kind='context'`` proteins are unlabelled genome content, whose
    cross-split sharing is counted (``context_sharing``) but is not a violation. Target homology
    clusters spanning splits are a violation under ``cluster_policy='refuse'`` and a warning under
    ``'mask'`` (the builder then masks the lower-precedence labelled examples).
    """
    violations, warnings = [], []
    unknown = [k for k in disjoint if k not in DISJOINT_KEYS]
    if unknown:
        raise PreparationError(f"unknown disjointness key(s) {unknown}; choose from {list(DISJOINT_KEYS)}")
    if cluster_policy not in ("refuse", "mask"):
        raise PreparationError(f"cluster_policy {cluster_policy!r} must be 'refuse' or 'mask'")
    for r in rows:
        if r["split"] not in SPLITS:
            violations.append(f"genome {r['genome_id']}: split {r['split']!r} not in {SPLITS}")

    by_genome = defaultdict(list)
    for r in rows:
        by_genome[r["genome_id"]].append(r)
    for g, rs in sorted(by_genome.items()):
        if len(rs) == 1:
            continue
        conflicting = [c for c in ("group", "source", "split") + TAXON_KEYS if len({x.get(c, "") for x in rs}) > 1]
        if conflicting:
            violations.append(f"genome {g} has {len(rs)} rows with conflicting {conflicting}"
                              + (" (declared in several splits)" if "split" in conflicting else ""))
        else:
            warnings.append(f"genome {g} is listed {len(rs)} times with identical metadata")

    def spans(key):
        owner = defaultdict(set)
        for r in rows:
            if r.get(key):
                owner[r[key]].add(r["split"])
        return {k: sorted(v) for k, v in owner.items() if len(v) > 1}

    for key in TAXON_KEYS:
        if key in disjoint:
            missing = sorted({r["genome_id"] for r in rows if not r.get(key)})
            if missing:
                violations.append(f"disjointness on '{key}' requested but taxonomy is unavailable for "
                                  f"{len(missing)} genome(s), e.g. {missing[:3]}; refusing to assume it")
    for key in ("group",) + TAXON_KEYS:
        crossing = spans(key)
        if not crossing:
            continue
        msg = f"{len(crossing)} {key}(s) span several splits, e.g. {sorted(crossing.items())[:3]}"
        (violations if key in disjoint else warnings).append(msg)

    registered_ok = True
    for src in sorted({r["source"] for r in rows}):
        entry = registry.get(src)
        if entry is None:
            violations.append(f"source {src!r} is not in the cohort registry (unregistered provenance)")
            registered_ok = False
            continue
        used = {r["split"] for r in rows if r["source"] == src}
        illegal = sorted(used - set(entry["allowed_splits"]))
        if illegal:
            violations.append(f"source {src!r} (role {entry['role']}) is not allowed in split(s) {illegal}; "
                              f"allowed: {entry['allowed_splits']}")
        if entry.get("genome_level") is False:
            violations.append(f"source {src!r} is a protein-level cohort (no genome layout); it cannot appear in a genome split table")
        if entry["role"] == "examined" and "test" in used:
            warnings.append(f"source {src!r} was examined during earlier development; its use as test data "
                            f"supports a disjointness claim, not an untouched-test claim")

    context_sharing = None
    homology_checked = False
    if proteins is not None:
        proteins = _normalise_proteins(proteins)
        split_of = {g: rs[0]["split"] for g, rs in by_genome.items()}
        owner = {}
        shared_ctx_ids = set()
        for genome, protein, kind in proteins:
            if genome not in split_of:
                violations.append(f"protein list names genome {genome!r} that has no split")
                continue
            s = split_of[genome]
            if protein in owner and owner[protein][1] != s:
                if kind == "target" and owner[protein][2] == "target":
                    violations.append(f"protein id {protein!r} appears in {owner[protein][1]} and {s}")
                else:
                    shared_ctx_ids.add(protein)
            owner.setdefault(protein, (genome, s, kind))
        context_sharing = {"protein_ids_shared_across_splits": len(shared_ctx_ids)}
        if protein_clusters is not None:
            homology_checked = True
            tgt, ctx = defaultdict(set), defaultdict(set)
            for genome, p, kind in proteins:
                if p in protein_clusters and genome in split_of:
                    (tgt if kind == "target" else ctx)[protein_clusters[p]].add(split_of[genome])
            unmapped_t = sorted(p for p, (_, _, k) in owner.items() if k == "target" and p not in protein_clusters)
            unmapped_c = sorted(p for p, (_, _, k) in owner.items() if k == "context" and p not in protein_clusters)
            if unmapped_t:
                violations.append(f"homology map does not cover {len(unmapped_t)} target protein(s), e.g. {unmapped_t[:3]}")
            if unmapped_c:
                warnings.append(f"homology map does not cover {len(unmapped_c)} context protein(s); their family sharing is not assessed")
            span = {c: sorted(s) for c, s in tgt.items() if len(s) > 1}
            if span:
                msg = f"{len(span)} homology cluster(s) of labelled target proteins span several splits, e.g. {sorted(span.items())[:3]}"
                if cluster_policy == "refuse":
                    violations.append(msg)
                else:
                    warnings.append(msg + "; labelled examples will be masked by build-table (masked_examples.tsv)")
            context_sharing["clusters_shared_across_splits"] = sum(1 for s in ctx.values() if len(s) > 1)
            context_sharing["note"] = ("context proteins may share conserved families across splits; this is reported, "
                                       "not a violation (only labelled target examples must be family-disjoint)")
        else:
            warnings.append("no homology-cluster map declared: protein-level family disjointness NOT assessed")
    else:
        warnings.append("no protein list declared: protein-id overlap NOT assessed")

    counts = {s: {"genomes": len({r["genome_id"] for r in rows if r["split"] == s}),
                  "groups": len({r["group"] for r in rows if r["split"] == s}),
                  "sources": sorted({r["source"] for r in rows if r["split"] == s})} for s in SPLITS}
    for s in SPLITS:
        if counts[s]["genomes"] == 0:
            violations.append(f"split '{s}' has no genomes")
    n_test = counts["test"]["genomes"]
    if 0 < n_test < min_test_genomes:
        warnings.append(f"only {n_test} test genome(s) (< {min_test_genomes}): genome-level uncertainty will be wide")
    dev_in_test = any(registry.get(r["source"], {}).get("role") == "development" and r["split"] == "test" for r in rows)
    examined_in_test = any(registry.get(r["source"], {}).get("role") in ("development", "examined")
                           and r["split"] == "test" for r in rows)

    # Precisely which axes were actually looked at. An axis that was not checked is listed as such.
    taxon_checked = {k: (k in disjoint and not any("disjointness on '%s'" % k in v for v in violations)) for k in TAXON_KEYS}
    checked = ["genome_id"] + (["group"] if "group" in disjoint else []) + (["source_registry"] if registered_ok else [])
    checked += [f"taxonomy:{k}" for k, ok in taxon_checked.items() if ok]
    if proteins is not None:
        checked.append("protein_id_overlap(targets)")
    if homology_checked:
        checked.append("homology_clusters(targets)")
    unchecked = [f"taxonomy:{k}" for k, ok in taxon_checked.items() if not ok]
    if proteins is None:
        unchecked.append("protein_id_overlap")
    if not homology_checked:
        unchecked.append("homology_clusters")
    unchecked += ["pretrained_model_leakage", "untouched_test_use"]
    required = {"genome_id", "group", "source_registry", "protein_id_overlap(targets)", "homology_clusters(targets)"}
    independent = not violations and not dev_in_test and required <= set(checked)
    claims = {
        "genome_disjoint_splits": all(len({x["split"] for x in rs}) == 1 for rs in by_genome.values()),
        "group_disjoint_splits": not spans("group"),
        "species_disjoint_splits": (not spans("species")) if taxon_checked["species"] or any(r["species"] for r in rows) else "not_assessed",
        "protein_overlap_checked": proteins is not None,
        "homology_disjoint_checked": homology_checked,
        "homology_scope": "labelled_target_proteins_only" if homology_checked else "not_assessed",
        "test_contains_development_data": dev_in_test,
        "test_contains_examined_or_development_data": examined_in_test,
        "all_sources_registered": registered_ok,
        "pretrained_model_leakage": "not_assessed",
        # The registry records KNOWN past use only; the absence of a recorded use is not evidence
        # that a cohort is untouched, so this claim is never derived from it.
        "untouched_test_claim": "not_assessed",
        "checked_axes": checked,
        "unchecked_axes": unchecked,
        "independent_test_claim_supportable": independent,
    }
    claims["independent_test_claim_text"] = (
        ("test is disjoint from train/calibration on: " if independent else "NOT a supported independence claim; checked: ")
        + ", ".join(checked) + "; not checked: " + ", ".join(unchecked))
    return {"violations": violations, "warnings": warnings, "counts": counts, "claims": claims,
            "context_sharing": context_sharing, "disjoint_keys": list(disjoint)}


# ------------------------------------------------------------------------------------ GO

class GoView:
    """GO closure over ``v3_metrics.load_obo`` (``[Typedef]`` stanzas skipped, alt ids and obsolete terms honoured)."""

    def __init__(self, obo_path):
        self.onto = load_obo(str(obo_path))

    def primary(self, term):
        return self.onto.primary(term)

    def obsolete(self, term) -> bool:
        return self.primary(term) in self.onto.obsolete

    def closure(self, term) -> set:
        return set(self.onto.ancestors(term)) - GO_ROOTS

    def is_complex_term(self, term) -> bool:
        anc = self.onto.ancestors(term)
        return GO_COMPLEX in anc and GO_CC_ROOT in anc


# ------------------------------------------------------------------------------------ labels

def load_label_policy(path) -> dict:
    pol = json.loads(Path(path).read_text())
    for ds, per_kind in pol.get("datasets", {}).items():
        for kind, p in per_kind.items():
            if kind not in ("EC", "GO"):
                raise PreparationError(f"label policy {ds}: unsupported function type {kind}")
            neg = p.get("negative_policy", "none")
            if neg not in NEGATIVE_POLICIES:
                raise PreparationError(f"label policy {ds}/{kind}: negative_policy {neg!r} not in {NEGATIVE_POLICIES}")
            if neg != "none" and not str(p.get("justification", "")).strip():
                raise PreparationError(f"label policy {ds}/{kind}: a negative policy needs a written justification")
            ev = p.get("allowed_evidence", "experimental")
            if ev not in EV.POLICIES and not isinstance(ev, list):
                raise PreparationError(f"label policy {ds}/{kind}: allowed_evidence must be one of {EV.POLICIES} or a list")
            for role, src in (p.get("role_truth") or {}).items():
                if role not in CANDIDATE_ROLES or not str(src).strip():
                    raise PreparationError(f"label policy {ds}/{kind}: role_truth needs a known role and a source string")
            if "complex_membership_from_cc" in p:
                if kind != "GO" or not str(p["complex_membership_from_cc"]).strip():
                    raise PreparationError(f"label policy {ds}/{kind}: complex_membership_from_cc is a GO-only opt-in "
                                           "and needs a written source")
    return pol


def evidence_ok(code: str, allowed) -> bool:
    return EV.allowed(code, allowed)


def read_labels(path, split_rows, policy, go: GoView | None = None):
    """Explicit positive/negative truth, preserving polarity and assay references.

    Missing labels stay unknown. Negatives are exact assertions, never propagated
    to GO ancestors, and need both a policy source and a row-level reference.
    """
    source_of = {r["genome_id"]: r["source"] for r in split_rows}
    truth = defaultdict(lambda: defaultdict(dict))
    dropped = Counter()
    with open(path, newline="") as fh:
        header = csv.DictReader(fh, delimiter="\t").fieldnames or []
    has_role = "role" in header
    for r in read_tsv(path, ("genome_id", "protein_id", "function_id", "evidence")):
        g = r["genome_id"].strip()
        if g not in source_of:
            dropped["label_genome_without_split"] += 1
            continue
        fid = r["function_id"].strip().replace("_", ":")
        kind = "GO" if fid.upper().startswith("GO:") else "EC"
        fid = function_id(kind, fid)
        p = policy.get("datasets", {}).get(source_of[g], {}).get(kind)
        if p is None:
            dropped[f"no_label_policy_for_{kind}"] += 1
            continue
        code = r["evidence"].strip()
        if not evidence_ok(code, p.get("allowed_evidence", "experimental")):
            dropped["evidence_not_allowed"] += 1
            if EV.category(code) == "curated_inferential":
                dropped["evidence_curated_inferential_not_experimental"] += 1
            continue
        qualifier = (r.get('qualifier') or '').strip()
        if qualifier not in ('', 'NOT'):
            raise PreparationError(f'unsupported truth qualifier {qualifier!r}; resolve its semantics before training')
        label_text = (r.get('label') or ('0' if qualifier == 'NOT' else '1')).strip()
        if label_text not in ('0', '1') or (qualifier == 'NOT' and label_text != '0'):
            raise PreparationError('truth label must be exactly 0/1 and agree with the NOT qualifier')
        label = int(label_text)
        reference = (r.get('reference') or '').strip()
        if label == 0 and (not isinstance(p.get('explicit_negative_truth'), str) or
                           not p['explicit_negative_truth'].strip() or not reference):
            raise PreparationError('explicit negative truth requires explicit_negative_truth policy source and a row reference')
        protein = r["protein_id"].strip()
        explicit = (r.get("role") or "").strip() if has_role else ""
        if kind == "GO" and go is not None:
            if go.obsolete(fid):
                dropped["obsolete_go_term"] += 1
                continue
            fid = go.primary(fid)
            ids = go.closure(fid) if label == 1 else {fid}
        else:
            ids = {fid}
        if explicit and explicit != DEFAULT_ROLE_MAP[kind][0]:
            src = (p.get("role_truth") or {}).get(explicit)
            if src is None:
                dropped["explicit_role_not_declared_in_policy"] += 1
                continue
            _add(truth[(g, kind, explicit)], protein, {fid}, code, f"explicit_role_truth:{src}", label, reference)   # exact term, no propagation
            continue
        _add(truth[(g, kind, DEFAULT_ROLE_MAP[kind][0])], protein, ids, code, "default_role", label, reference)
        if label == 1 and kind == "GO" and p.get("complex_membership_from_cc") and go is not None and go.is_complex_term(fid):
            cids = {t for t in ids if go.is_complex_term(t)}
            _add(truth[(g, kind, COMPLEX_ROLE)], protein, cids, code,
                 f"cc_complex_membership:{p['complex_membership_from_cc']}")
    return truth, dropped


def _add(bucket, protein, ids, code, basis, label=1, reference=''):
    for t in ids:
        rec = bucket[protein].setdefault(t, {"codes": set(), "basis": basis, "label": label, "references": set()})
        if rec['label'] != label:
            raise PreparationError(f'conflicting positive/negative truth for {protein}/{t}')
        rec["codes"].add(code)
        if reference:
            rec['references'].add(reference)


# ------------------------------------------------------------------------------------ feature rows

def _finite(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def validate_export_dir(d: Path, contract: dict) -> int:
    """Check one ``<genome>/`` export against the production writer's contract
    (``gspa.integration.learned.FeatureExporter``); returns the number of rows."""
    for name in (EXPORT_ROWS_FILE, EXPORT_CONTRACT_FILE, EXPORT_MANIFEST_FILE):
        if not (d / name).exists():
            raise PreparationError(f"{d}: missing {name}")
    manifest = json.loads((d / EXPORT_MANIFEST_FILE).read_text())
    if manifest.get("schema") != EXPORT_MANIFEST_SCHEMA:
        raise PreparationError(f"{d}: export manifest schema {manifest.get('schema')!r} is not {EXPORT_MANIFEST_SCHEMA}")
    if manifest.get("labels_included") is not False:
        raise PreparationError(f"{d}: export manifest does not state labels_included=false")
    if manifest.get("rows_file") != EXPORT_ROWS_FILE or manifest.get("contract_file") != EXPORT_CONTRACT_FILE:
        raise PreparationError(f"{d}: export manifest names unexpected rows/contract files")
    if manifest.get("rows_sha256") != sha256_file(d / EXPORT_ROWS_FILE):
        raise PreparationError(f"{d}: candidate_features.jsonl does not match its export manifest hash")
    if manifest.get("contract_fingerprint") != contract["fingerprint"] or manifest.get("contract_id") != contract["id"]:
        raise PreparationError(f"{d}: exported with a different feature contract than the one supplied")
    shipped = json.loads((d / EXPORT_CONTRACT_FILE).read_text())
    if (shipped.get("id") != contract["id"] or shipped.get("fingerprint") != contract["fingerprint"]
            or [(f["name"], f["group"]) for f in shipped.get("features", [])] !=
            [(f["name"], f["group"]) for f in contract["features"]]):
        raise PreparationError(f"{d}: exported feature_contract.json differs from the supplied contract")
    n = 0
    with open(d / EXPORT_ROWS_FILE) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            where = f"{d}/{EXPORT_ROWS_FILE}:{ln}"
            if row.get("schema") != EXPORT_ROW_SCHEMA:
                raise PreparationError(f"{where}: row schema {row.get('schema')!r} is not {EXPORT_ROW_SCHEMA}")
            if row.get("contract") != {"id": contract["id"], "fingerprint": contract["fingerprint"]}:
                raise PreparationError(f"{where}: row contract id/fingerprint mismatch")
            for key in ("protein_id", "function_type", "function_id", "function_key", "role", "candidate_key"):
                if not isinstance(row.get(key), str) or not row[key]:
                    raise PreparationError(f"{where}: missing or empty '{key}'")
            if row["role"] not in CANDIDATE_ROLES:
                raise PreparationError(f"{where}: unknown role {row['role']!r}")
            if row.get("cohort_context_mode", "UNIVERSAL") not in ("UNIVERSAL", "PATHWAY", "ORGANISM"):
                raise PreparationError(f"{where}: invalid cohort context mode")
            fkey = f"{row['protein_id']}|{row['function_type']}|{row['function_id']}"
            if row["function_key"] != fkey or row["candidate_key"] != f"{fkey}|{row['role']}":
                raise PreparationError(f"{where}: function_key/candidate_key do not follow the writer's layout")
            if "genome_id" in row and row["genome_id"] != d.name:
                raise PreparationError(f"{where}: genome_id {row['genome_id']!r} differs from its export directory {d.name!r}")
            for forbidden in ("label", "split", "group", "example_id"):
                if forbidden in row:
                    raise PreparationError(f"{where}: exporter row already carries '{forbidden}'")
            feats = row.get("features")
            if not isinstance(feats, dict) or list(feats) != contract["names"]:
                raise PreparationError(f"{where}: features are not exactly the contract columns in order")
            bad = [k for k, v in feats.items() if not _finite(v)]
            if bad:
                raise PreparationError(f"{where}: non-finite feature value(s) {bad[:3]}")
            ctxs = row.get("available_contexts")
            if (not isinstance(ctxs, dict) or set(ctxs) not in
                    (set(EXPORT_CONTEXT_KEYS), set(EXPORT_CONTEXT_KEYS) | {'local_reaction'})
                    or not all(isinstance(v, bool) for v in ctxs.values())):
                raise PreparationError(f"{where}: available_contexts does not match the writer's layout")
            if 'local_reaction' in ctxs and ctxs['local_reaction'] != (feats.get('cohort.w.local_available', 0) > 0):
                raise PreparationError(f"{where}: available_contexts.local_reaction disagrees with the local availability feature")
            if not isinstance(row.get("notes"), list):
                raise PreparationError(f"{where}: notes must be a list")
            if row.get('quality_evidence'):
                raise PreparationError(f"{where}: gene-model quality observations are not supported for training; retain them and extend the model contract before joining labels")
            n += 1
    if manifest.get("rows") != n:
        raise PreparationError(f"{d}: export manifest says {manifest.get('rows')} rows, file has {n}")
    return n


def iter_feature_rows(feature_dir, contract):
    """Yield ``(genome_dir_name, row)`` after verifying every export against its manifest and the writer's layout."""
    feature_dir = Path(feature_dir)
    dirs = sorted(p for p in feature_dir.iterdir() if p.is_dir() and (p / EXPORT_ROWS_FILE).exists())
    if not dirs:
        raise PreparationError(f"{feature_dir}: no <genome>/candidate_features.jsonl found")
    for d in dirs:
        validate_export_dir(d, contract)
        with open(d / EXPORT_ROWS_FILE) as fh:
            for line in fh:
                if line.strip():
                    yield d.name, json.loads(line)


def load_locus_groups(path, policy_declares_complete):
    """candidate_key -> (locus_group_id) for groups whose truth is declared complete."""
    out = {}
    if path is None:
        return out
    seen, states = set(), {}
    for r in read_tsv(path, ("genome_id", "locus_group_id", "candidate_key", "truth_complete")):
        gid, group, candidate, complete = (r[k].strip() for k in ('genome_id', 'locus_group_id', 'candidate_key', 'truth_complete'))
        if not gid or not group or not candidate or complete.lower() not in ('true', 'false'):
            raise PreparationError('locus declarations need nonempty identities and truth_complete true/false')
        key, gkey = (gid, candidate), (gid, group)
        if key in seen:
            raise PreparationError(f'duplicate candidate or multiple locus groups for {key}')
        seen.add(key)
        state = complete.lower() == 'true'
        if gkey in states and states[gkey] != state:
            raise PreparationError(f'locus group {gkey} has conflicting truth_complete declarations')
        states[gkey] = state
        if state:
            if not isinstance(policy_declares_complete, str) or not policy_declares_complete.strip():
                raise PreparationError("locus groups declared complete but the label policy has no "
                                       "'locus_truth_complete' justification")
            out[key] = group
    return out


def certify_locus_groups(rows, declarations, candidate_universe, truth_source):
    """Certify complete production competition groups after all label/family filtering.

    Production compares all candidates of a function/role within one genome.
    Retain incomplete groups' marginal rows, but remove their locus supervision.
    """
    declared = defaultdict(set)
    for (genome, candidate), group in declarations.items():
        declared[(genome, group)].add(candidate)
    by_candidate = {(r['genome_id'], r['candidate_key']): r for r in rows}
    report = []
    for (genome, group), members in sorted(declared.items()):
        signatures = set()
        for candidate in members:
            parts = candidate.split('|')
            if len(parts) != 4:
                raise PreparationError(f'invalid locus candidate key {candidate!r}')
            signatures.add(tuple(parts[1:]))
        if len(signatures) != 1:
            raise PreparationError(f'locus {genome}/{group} mixes function identities or roles')
        signature = next(iter(signatures))
        universe = candidate_universe.get((genome, *signature), set())
        retained = {c for c in members if (genome, c) in by_candidate}
        reasons = []
        if members != universe:
            reasons.append('candidate_universe_incomplete')
        if retained != members:
            reasons.append('members_not_labelled_or_masked')
        if any(by_candidate[(genome, c)]['label_basis'] == 'closed_world_unconfirmed' for c in retained):
            reasons.append('unconfirmed_negative_labels')
        example_ids = sorted(by_candidate[(genome, c)]['example_id'] for c in retained)
        for c in retained:
            row = by_candidate[(genome, c)]
            row['locus_group_id'] = None if reasons else group
            if not reasons:
                row['locus_members'] = example_ids
                row['locus_truth_source'] = truth_source
        report.append({'genome_id': genome, 'locus_group_id': group, 'eligible': not reasons,
                       'reasons': reasons, 'declared_members': sorted(members), 'retained_members': sorted(retained),
                       'exported_members': sorted(universe)})
    return {'declared_groups': len(declared), 'eligible_groups': sum(r['eligible'] for r in report),
            'groups': report, 'scope': 'all exported candidates of one function/role in one genome; complete explicit truth only'}


# ------------------------------------------------------------------------------------ family masking

def read_protein_clusters(path) -> dict:
    """``protein_id -> cluster`` (global ids) or, with a ``genome_id`` column, ``(genome_id, protein_id) -> cluster``."""
    rows = read_tsv(path, ("protein_id", "cluster"))
    has_genome = "genome_id" in rows[0] if rows else False
    return {((r["genome_id"].strip(), r["protein_id"].strip()) if has_genome else r["protein_id"].strip()): r["cluster"].strip()
            for r in rows}


def _cluster_of(clusters, genome, protein):
    return clusters.get((genome, protein), clusters.get(protein))


def mask_homologous_examples(rows_out, clusters):
    """Mask labelled examples whose family also occurs among higher-precedence labelled examples.

    Returns ``(kept, masked, report)``. ``masked`` rows carry ``masked_reason``/``cluster``/
    ``clashes_with``. Raises when a labelled example has no cluster (never silently kept or dropped).
    """
    unmapped = sorted({(r["genome_id"], r["protein_id"]) for r in rows_out
                       if _cluster_of(clusters, r["genome_id"], r["protein_id"]) is None})
    if unmapped:
        raise PreparationError(f"homology map does not cover {len(unmapped)} labelled protein(s), e.g. {unmapped[:3]}")
    split_by_cluster = defaultdict(set)
    for r in rows_out:
        split_by_cluster[_cluster_of(clusters, r["genome_id"], r["protein_id"])].add(r["split"])
    kept, masked = [], []
    for r in rows_out:
        c = _cluster_of(clusters, r["genome_id"], r["protein_id"])
        clash = sorted(s for s in MASK_PRECEDENCE[r["split"]] if s in split_by_cluster[c])
        if clash:
            masked.append({**r, "cluster": c, "clashes_with": clash,
                           "masked_reason": "family_in_" + "+".join(clash) + "_labelled_examples"})
        else:
            kept.append(r)
    hard_test = sum(1 for r in rows_out if r["split"] == "test"
                    and "train" in split_by_cluster[_cluster_of(clusters, r["genome_id"], r["protein_id"])])
    report = {"labelled_examples_before": len(rows_out), "masked": len(masked), "kept": len(kept),
              "masked_by_split": {s: sum(1 for m in masked if m["split"] == s) for s in SPLITS},
              "masked_positives_by_split": {s: sum(m["label"] for m in masked if m["split"] == s) for s in SPLITS},
              "test_examples_with_family_in_train_retained": hard_test,
              "family_disjoint_verified": True}
    left = defaultdict(set)
    for r in kept:
        left[_cluster_of(clusters, r["genome_id"], r["protein_id"])].add(r["split"])
    if any(len(s) > 1 for s in left.values()):
        raise PreparationError("internal error: kept labelled examples still share families across splits")
    for s in SPLITS:
        if any(r["split"] == s for r in rows_out) and not any(r["split"] == s and r["label"] == 1 for r in kept):
            raise PreparationError(f"family masking left split '{s}' without positive labelled examples "
                                   f"({report['masked_by_split'][s]} masked); refusing to continue")
    return kept, masked, report


# ------------------------------------------------------------------------------------ build

def build_table(feature_dir, contract_path, labels_path, split_path, policy_path, out_dir, *,
                registry_path=REGISTRY_DEFAULT, go_obo=None, locus_groups=None, protein_list=None,
                protein_clusters=None, example_clusters=None, disjoint=DEFAULT_DISJOINT) -> dict:
    """``protein_list``/``protein_clusters`` feed the pre-build audit (kind=context|target per protein);
    ``example_clusters`` (``protein_id[,genome_id],cluster``) enables family masking of labelled examples."""
    out_dir = Path(out_dir)
    contract = json.loads(Path(contract_path).read_text())
    contract["names"] = [f["name"] for f in contract["features"]]
    if contract.get("fingerprint") not in (None, contract_fingerprint(contract)):
        raise PreparationError(f"{contract_path}: stated fingerprint does not match its feature list")
    contract["fingerprint"] = contract_fingerprint(contract)
    registry = load_registry(registry_path)
    split_rows = read_split_table(split_path)
    proteins = None
    if protein_list:
        proteins = [(r["genome_id"], r["protein_id"], r.get("kind") or "target")
                    for r in read_tsv(protein_list, ("genome_id", "protein_id"))]
    clusters = read_protein_clusters(protein_clusters) if protein_clusters else None
    audit = audit_splits(split_rows, registry, disjoint=disjoint, proteins=proteins,
                         protein_clusters=clusters if clusters and all(isinstance(k, str) for k in clusters) else None,
                         cluster_policy="mask" if example_clusters else "refuse")
    if audit["violations"]:
        raise PreparationError("split audit failed: " + "; ".join(audit["violations"]))
    policy = load_label_policy(policy_path)
    go = GoView(go_obo) if go_obo else None
    truth, label_drops = read_labels(labels_path, split_rows, policy, go)
    locus = load_locus_groups(locus_groups, policy.get("locus_truth_complete"))
    info = {r["genome_id"]: r for r in split_rows}

    rows_out, excluded = [], Counter()
    seen = set()
    positives_with_candidate = defaultdict(set)
    candidate_universe = defaultdict(set)
    for genome, row in iter_feature_rows(feature_dir, contract):
        gid = row.get("genome_id") or genome
        if gid not in info:
            excluded["genome_without_split"] += 1
            continue
        kind = row["function_type"]
        if kind not in ("EC", "GO"):
            excluded["unsupported_function_type"] += 1
            continue
        key = (gid, row["candidate_key"])
        if key in seen:
            raise PreparationError(f"duplicate candidate {key}")
        seen.add(key)
        candidate_universe[(gid, kind, row['function_id'], row['role'])].add(row['candidate_key'])
        fid = function_id(kind, row["function_id"])
        if kind == "GO" and go is not None:
            if go.obsolete(fid):
                excluded["obsolete_go_candidate"] += 1
                continue
            fid = go.primary(fid)
        role, default_role = row["role"], DEFAULT_ROLE_MAP[kind][0]
        info_g = info[gid]
        pol = policy["datasets"].get(info_g["source"], {}).get(kind)
        if pol is None:
            excluded["no_label_policy"] += 1
            continue
        protein = row["protein_id"]
        rec = truth.get((gid, kind, role), {}).get(protein, {}).get(fid)
        if rec is not None:
            label = rec['label']
            basis, codes, role_basis = ('positive_truth' if label else 'explicit_negative_truth'), sorted(rec["codes"]), rec["basis"]
            if label:
                positives_with_candidate[(gid, kind, role)].add((protein, fid))
        elif role != default_role:
            excluded["role_not_labelled"] += 1
            if role == COMPLEX_ROLE:
                excluded["complex_subunit_without_role_truth"] += 1
            continue
        else:
            labelled_protein = any(r['label'] == 1 for r in truth.get((gid, kind, default_role), {}).get(protein, {}).values())
            if labelled_protein and pol.get("negative_policy") == "closed_world_declared":
                label, basis, codes, role_basis = 0, "closed_world_unconfirmed", [], "default_role"
            else:
                excluded["no_truth_not_negative" if not labelled_protein else "unconfirmed_not_negative"] += 1
                continue
        rows_out.append({
            "example_id": f"{gid}|{row['candidate_key']}", "genome_id": gid, "group": info_g["group"],
            "split": info_g["split"], "label": label, "label_basis": basis, "dataset": info_g["source"],
            "label_evidence": codes, "label_evidence_categories": sorted({EV.category(c) for c in codes}),
            "label_role_basis": role_basis,
            "label_references": sorted(rec['references']) if rec else [],
            "candidate_key": row["candidate_key"], "protein_id": protein,
            "function_id": fid, "function_type": kind, "role": role,
            "locus_group_id": locus.get(key), "features": row["features"]})
        if 'sequence_sources' in row:
            rows_out[-1]['sequence_sources'] = row['sequence_sources']
        if 'sequence_contracts' in row:
            rows_out[-1]['sequence_contracts'] = row['sequence_contracts']
        if 'cohort_context_mode' in row:
            rows_out[-1]['cohort_context_mode'] = row['cohort_context_mode']

    masking = None
    masked = []
    if example_clusters:
        rows_out, masked, masking = mask_homologous_examples(rows_out, read_protein_clusters(example_clusters))
    locus_audit = certify_locus_groups(rows_out, locus, candidate_universe, policy.get('locus_truth_complete'))
    out_dir.mkdir(parents=True, exist_ok=True)
    table = out_dir / "training_table.jsonl"
    with open(table, "w") as fh:
        for r in rows_out:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
    masked_path = None
    if example_clusters:
        masked_path = out_dir / "masked_examples.tsv"
        with open(masked_path, "w") as fh:
            fh.write("example_id\tsplit\tlabel\tlabel_basis\tprotein_id\tfunction_id\trole\tcluster\tclashes_with\tmasked_reason\n")
            for m in masked:
                fh.write("\t".join(str(x) for x in (m["example_id"], m["split"], m["label"], m["label_basis"], m["protein_id"],
                                                    m["function_id"], m["role"], m["cluster"], "+".join(m["clashes_with"]),
                                                    m["masked_reason"])) + "\n")
    per_split = {}
    for s in SPLITS:
        sub = [r for r in rows_out if r["split"] == s]
        per_split[s] = {"rows": len(sub), "positives": sum(r["label"] for r in sub),
                        "closed_world_negatives": sum(1 for r in sub if r["label_basis"] == "closed_world_unconfirmed"),
                        "explicit_negatives": sum(1 for r in sub if r["label_basis"] == "explicit_negative_truth"),
                        "genomes": len({r["genome_id"] for r in sub}), "groups": len({r["group"] for r in sub}),
                        "locus_group_rows": sum(1 for r in sub if r["locus_group_id"]),
                        "role_basis": dict(Counter(r["label_role_basis"].split(":")[0] for r in sub if r["label"] == 1))}
    unreachable = {}
    for (gid, kind, role), by_protein in truth.items():
        total = sum(r['label'] == 1 for fs in by_protein.values() for r in fs.values())
        reached = len(positives_with_candidate.get((gid, kind, role), ()))
        unreachable[f"{gid}/{kind}/{role}"] = {"truth_pairs": total, "with_candidate_row": reached}
    manifest = {
        "schema": MANIFEST_SCHEMA, "table": str(table), "table_sha256": sha256_file(table),
        "inputs": {"feature_dir": str(feature_dir), "contract": str(contract_path),
                   "contract_fingerprint": contract["fingerprint"], "labels": str(labels_path),
                   "labels_sha256": sha256_file(labels_path), "splits": str(split_path),
                   "splits_sha256": sha256_file(split_path), "label_policy": str(policy_path),
                   "label_policy_sha256": sha256_file(policy_path), "registry": str(registry_path),
                   "locus_groups": str(locus_groups) if locus_groups else None,
                   "locus_groups_sha256": sha256_file(locus_groups) if locus_groups else None,
                   "go_obo": str(go_obo) if go_obo else None,
                   "go_obo_sha256": sha256_file(go_obo) if go_obo else None,
                   "go_propagation": "is_a+part_of" if go_obo else "exact_term",
                   "example_clusters": str(example_clusters) if example_clusters else None,
                   "example_clusters_sha256": sha256_file(example_clusters) if example_clusters else None},
        "locus_audit": locus_audit,
        "split_audit": audit, "per_split": per_split, "excluded_rows": dict(excluded),
        "label_records_dropped": dict(label_drops),
        "family_masking": masking if masking else {"applied": False,
                                                   "note": "no --example-clusters: labelled-example family disjointness NOT enforced"},
        "masked_examples_file": str(masked_path) if masked_path else None,
        "truth_recall_ceiling": unreachable,
        "label_basis_note": ("explicit_negative_truth rows carry polarity, assay policy and references; "
                             "closed_world_unconfirmed rows are unconfirmed, not verified-false, and cannot supervise locus null outcomes"),
        "dataset_status": "real" if all(registry[s]["role"] != "synthetic" for s in {r["source"] for r in split_rows}) else "synthetic",
    }
    (out_dir / "training_table_manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    return manifest


# ------------------------------------------------------------------------------------ CLI

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("audit-splits")
    a.add_argument("--splits", required=True)
    a.add_argument("--registry", default=str(REGISTRY_DEFAULT))
    a.add_argument("--proteins", help="TSV genome_id<TAB>protein_id[<TAB>kind: target|context]; protein_id must be a cross-genome identity")
    a.add_argument("--protein-clusters", help="TSV protein_id<TAB>cluster")
    a.add_argument("--cluster-policy", choices=("refuse", "mask"), default="refuse")
    a.add_argument("--disjoint", default=",".join(DEFAULT_DISJOINT), help=f"comma list from {list(DISJOINT_KEYS)}")
    a.add_argument("--out", help="write the report JSON here")
    b = sub.add_parser("build-table")
    b.add_argument("--feature-dir", required=True)
    b.add_argument("--contract", required=True)
    b.add_argument("--labels", required=True)
    b.add_argument("--splits", required=True)
    b.add_argument("--label-policy", required=True)
    b.add_argument("--registry", default=str(REGISTRY_DEFAULT))
    b.add_argument("--go-obo")
    b.add_argument("--locus-groups")
    b.add_argument("--proteins")
    b.add_argument("--protein-clusters")
    b.add_argument("--example-clusters", help="TSV protein_id[<TAB>genome_id]<TAB>cluster covering every labelled protein; "
                                              "enables family masking of labelled examples (masked_examples.tsv)")
    b.add_argument("--disjoint", default=",".join(DEFAULT_DISJOINT))
    b.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    try:
        disjoint = tuple(k.strip() for k in args.disjoint.split(",") if k.strip())
        if args.cmd == "audit-splits":
            proteins = clusters = None
            if args.proteins:
                proteins = [(r["genome_id"], r["protein_id"], r.get("kind") or "target")
                            for r in read_tsv(args.proteins, ("genome_id", "protein_id"))]
            if args.protein_clusters:
                clusters = read_protein_clusters(args.protein_clusters)
            report = audit_splits(read_split_table(args.splits), load_registry(args.registry), disjoint=disjoint,
                                  proteins=proteins, protein_clusters=clusters, cluster_policy=args.cluster_policy)
            text = json.dumps(report, indent=1, sort_keys=True)
            if args.out:
                Path(args.out).write_text(text + "\n")
            print(text)
            return 2 if report["violations"] else 0
        manifest = build_table(args.feature_dir, args.contract, args.labels, args.splits, args.label_policy, args.out_dir,
                               registry_path=args.registry, go_obo=args.go_obo, locus_groups=args.locus_groups,
                               protein_list=args.proteins, protein_clusters=args.protein_clusters,
                               example_clusters=args.example_clusters, disjoint=disjoint)
        print(json.dumps({"per_split": manifest["per_split"], "excluded_rows": manifest["excluded_rows"],
                          "family_masking": manifest["family_masking"]}, indent=1))
        return 0
    except PreparationError as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
