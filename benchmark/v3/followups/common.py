#!/usr/bin/env python3
"""Shared loaders and bookkeeping for the Stage 2B follow-up experiments."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import platform
import random
import sys
import tarfile
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from offgrowth.mask import (  # noqa: E402
    CandidateClaim,
    TransporterHit,
    TruthCase,
    base_seed_rxn_id,
    load_extra_cases,
    load_gapsmith_transporter_tbl,
    load_price148_cases,
    load_stage1_candidates_from_tsv,
    protein_aliases,
    strip_version,
)
from offgrowth.panel import normalize_ec, is_complete_ec  # noqa: E402

REPO_ROOT = V3_DIR.parents[1]
DESIGN_PATH = REPO_ROOT / "tasks" / "v3-stage2b-followups-design-2026-10-31.md"
DESIGN_SHA256 = "5a468ece2c892f3cc6b47622f8ef77f444d1e40f04bb11ca6b38ac5779ecb768"  # whole file at freezing time
DESIGN_BODY_SHA256 = "25a1a8676078e6e32d7c30120bb4024d7630c26509e18253ffcffdb612ef5f14"  # everything before "## 9. Deviations"
DEVIATIONS_MARKER = b"## 9. Deviations"

# Production constants of gspa.integration.fba.LocusResponsibilityUpdater / the sidecar (never tuned here).
BETA1 = 1.0
BETA3 = 1.5
TAU_NULL = 0.05
GAMMA_FBA = 2.2
SYNTENY_WINDOW = 10
NEIGHBOUR_MIN_SCORE = 0.25  # the candidate cut of the production reactions.tsv export
CANDIDATE_FLOOR = 1e-5      # the sidecar's original candidate floor
SHUFFLE_SEEDS: Tuple[int, ...] = (1, 2, 3, 4, 5)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def design_receipt(design_path: Optional[Path] = None) -> Dict[str, object]:
    """The frozen design's hash (the committed value and the file's current hash, which must agree)."""
    path = Path(design_path) if design_path is not None else DESIGN_PATH
    current = sha256_file(path) if path.exists() else None
    body = None
    if path.exists():
        raw = path.read_bytes()
        body = hashlib.sha256(raw[: raw.index(DEVIATIONS_MARKER)] if DEVIATIONS_MARKER in raw else raw).hexdigest()
    return {
        "design_path": str(path),
        "design_sha256_frozen": DESIGN_SHA256,
        "design_sha256_current": current,
        "design_body_sha256_frozen": DESIGN_BODY_SHA256,
        "design_body_sha256_current": body,
        "design_unchanged": current == DESIGN_SHA256,
        "design_body_unchanged": body == DESIGN_BODY_SHA256,
    }


def environment_receipt(inputs: Mapping[str, Path]) -> Dict[str, object]:
    hashes = {k: sha256_file(p) for k, p in sorted(inputs.items()) if Path(p).is_file()}
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "input_sha256": hashes,
        **design_receipt(),
    }


# ---------------------------------------------------------------------------
# Genome layout and synteny
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Locus:
    protein_id: str
    contig: str
    rank: int
    start: int
    strand: str


class GenomeLayout:
    """Proteins ordered by start coordinate per contig; neighbours are ranks within +/- window."""

    def __init__(self, loci: Sequence[Locus]):
        self.by_protein: Dict[str, Locus] = {}
        self.by_contig: Dict[str, List[Locus]] = {}
        for locus in loci:
            self.by_contig.setdefault(locus.contig, []).append(locus)
        for contig, items in self.by_contig.items():
            items.sort(key=lambda l: (l.start, l.protein_id))
            ranked = [Locus(l.protein_id, contig, i, l.start, l.strand) for i, l in enumerate(items)]
            self.by_contig[contig] = ranked
            for l in ranked:
                for alias in protein_aliases(l.protein_id):
                    self.by_protein.setdefault(alias, l)

    @classmethod
    def from_tsv(cls, path: Path) -> "GenomeLayout":
        loci: List[Locus] = []
        with Path(path).open("r", encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    start = int(row["start"])
                except (KeyError, ValueError):
                    continue
                loci.append(Locus(row["protein_id"], row["contig"], 0, start, row.get("strand", "")))
        return cls(loci)

    def locus(self, protein_id: str) -> Optional[Locus]:
        return self.by_protein.get(protein_id) or self.by_protein.get(strip_version(protein_id))

    def neighbours(self, protein_id: str, window: int = SYNTENY_WINDOW) -> List[str]:
        """Protein ids within +/- window loci of protein_id on its contig (protein_id itself excluded)."""
        locus = self.locus(protein_id)
        if locus is None:
            return []
        items = self.by_contig[locus.contig]
        lo = max(0, locus.rank - window)
        hi = min(len(items), locus.rank + window + 1)
        return [items[i].protein_id for i in range(lo, hi) if items[i].protein_id != locus.protein_id]

    def shuffled(self, seed: int) -> "GenomeLayout":
        """Control: the same positions with protein identities permuted (seeded)."""
        rng = random.Random(seed)
        all_loci = [l for items in self.by_contig.values() for l in items]
        ids = [l.protein_id for l in all_loci]
        rng.shuffle(ids)
        return GenomeLayout([Locus(pid, l.contig, 0, l.start, l.strand) for pid, l in zip(ids, all_loci)])

    def n_proteins(self) -> int:
        return sum(len(v) for v in self.by_contig.values())


def load_kegg_ec_modules(path: Path) -> Dict[str, Set[str]]:
    """{ec: {module ids}} from a KEGG module table with pathway_id/module_id and ec_number/ec columns."""
    ec_modules: Dict[str, Set[str]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mod = (row.get("pathway_id") or row.get("module_id") or "").strip()
            ec = normalize_ec(row.get("ec_number") or row.get("ec") or "")
            if mod and is_complete_ec(ec):
                ec_modules.setdefault(ec, set()).add(mod)
    return ec_modules


def reaction_to_ecs(ec_to_reactions: Mapping[str, Sequence[str]]) -> Dict[str, Set[str]]:
    rxn_ecs: Dict[str, Set[str]] = {}
    for ec, rxns in ec_to_reactions.items():
        for r in rxns:
            rxn_ecs.setdefault(base_seed_rxn_id(r), set()).add(normalize_ec(ec))
    return rxn_ecs


# ---------------------------------------------------------------------------
# Stored off-growth outputs
# ---------------------------------------------------------------------------


def read_genome_results(tar_path: Path) -> Dict[str, object]:
    """Read genome_results.json from a packed per-genome off-growth output tar."""
    with tarfile.open(tar_path, "r") as tf:
        for member in tf.getmembers():
            if member.name.endswith("genome_results.json"):
                fh = tf.extractfile(member)
                if fh is None:
                    break
                return json.loads(fh.read().decode("utf-8"))
    raise FileNotFoundError(f"genome_results.json not in {tar_path}")


def stored_tests_by_arm(arm_results: Mapping[str, Mapping[str, object]]) -> Dict[str, List[Dict[str, object]]]:
    """Per arm (A1, A2, A3): the stored per-test results (test_id, target_id, added_reactions)."""
    out: Dict[str, List[Dict[str, object]]] = {}
    for arm in ("A0", "A1", "A2", "A3"):
        ar = arm_results.get(arm)
        if not ar:
            continue
        out[arm] = [
            {
                "test_id": str(t["test_id"]),
                "target_id": str(t.get("target_id", "")),
                "added_reactions": tuple(str(x) for x in t.get("added_reactions", ())),
                "gapfill_ran": bool(t.get("gapfill_ran", False)),
            }
            for t in ar.get("test_results", ())
        ]
    return out


def added_reactions_with_tests(
    tests: Iterable[Dict[str, object]],
) -> Dict[str, Set[str]]:
    """{base reaction id: {test ids that added it}}."""
    out: Dict[str, Set[str]] = {}
    for t in tests:
        for r in t["added_reactions"]:
            out.setdefault(base_seed_rxn_id(str(r)), set()).add(str(t["test_id"]))
    return out


def load_known_experimental_ecs(
    genome_payload: Mapping[str, object],
    swissprot_tsv: Optional[Path],
    genome_id: str,
) -> Dict[str, Tuple[str, ...]]:
    """Price-148 truth (recorded in the stored payload) plus Swiss-Prot experimental labels for the genome."""
    known: Dict[str, Set[str]] = {}
    for prot, ecs in (genome_payload.get("known_experimental_ecs") or {}).items():
        known.setdefault(str(prot), set()).update(str(e) for e in ecs)
    if swissprot_tsv is not None and Path(swissprot_tsv).exists():
        with Path(swissprot_tsv).open("r", encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if (row.get("genome_id") or "").strip() != genome_id:
                    continue
                if (row.get("tier") or "experimental").strip().lower() != "experimental":
                    continue
                pid = (row.get("protein_id") or "").strip()
                ec = (row.get("function_id") or row.get("ec") or "").strip()
                if pid and ec:
                    known.setdefault(pid, set()).add(ec)
    return {p: tuple(sorted(v)) for p, v in known.items()}


def load_cases(manifest: Path, swissprot_tsv: Optional[Path], genome_id: str) -> List[TruthCase]:
    cases = [c for c in load_price148_cases(manifest) if c.genome_id == genome_id]
    if swissprot_tsv is not None and Path(swissprot_tsv).exists():
        cases.extend(load_extra_cases(swissprot_tsv, genome_id=genome_id))
    return cases


def load_candidates(path: Path) -> Tuple[CandidateClaim, ...]:
    return load_stage1_candidates_from_tsv([Path(path)])


def transporter_substrates(transporters: Sequence[TransporterHit]) -> Dict[str, Set[str]]:
    """{seed compound id: {protein ids with a transporter hit for it}}."""
    out: Dict[str, Set[str]] = {}
    for t in transporters:
        if t.seed_id:
            out.setdefault(t.seed_id, set()).add(t.protein_id)
    return out


def write_json(path: Path, obj: object) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], columns: Sequence[str]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
