#!/usr/bin/env python3
"""Per-genome label-free inputs shared by the re-scoring experiments (1, 2, 3, 6)."""

from __future__ import annotations

from dataclasses import dataclass, field
import io
from pathlib import Path
import tarfile
import tempfile
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from .common import (
    GenomeLayout,
    load_candidates,
    load_kegg_ec_modules,
    reaction_to_ecs,
    transporter_substrates,
)
from offgrowth.mask import CandidateClaim, TransporterHit, load_gapsmith_transporter_tbl
from offgrowth.panel import (
    PanelEntry,
    SeedMetabolite,
    SeedReaction,
    build_ec_to_reactions_map,
    load_panel_tsv,
    load_seed_metabolites,
    load_seed_reactions,
)


@dataclass
class ReferenceTables:
    reactions: Mapping[str, SeedReaction]
    metabolites: Mapping[str, SeedMetabolite]
    ec_to_reactions: Mapping[str, Sequence[str]]
    rxn_ecs: Mapping[str, Set[str]]
    ec_modules: Mapping[str, Set[str]]
    panel: Tuple[PanelEntry, ...]
    panel_by_id: Dict[str, PanelEntry] = field(default_factory=dict)

    @classmethod
    def load(
        cls,
        seed_reactions: Path,
        ec_aliases: Optional[Path],
        seed_metabolites: Optional[Path],
        kegg_modules: Path,
        panel_tsv: Path,
    ) -> "ReferenceTables":
        reactions = load_seed_reactions(seed_reactions, ec_aliases)
        ec_to_rxns = build_ec_to_reactions_map(reactions, usable_only=False)
        mets = load_seed_metabolites(seed_metabolites) if seed_metabolites is not None else {}
        panel = load_panel_tsv(panel_tsv)
        return cls(
            reactions=reactions,
            metabolites=mets,
            ec_to_reactions=ec_to_rxns,
            rxn_ecs=reaction_to_ecs(ec_to_rxns),
            ec_modules=load_kegg_ec_modules(kegg_modules),
            panel=panel,
            panel_by_id={p.panel_id: p for p in panel},
        )

    def substrate_by_test(self) -> Dict[str, str]:
        """A1 test id -> SEED compound of its sole source."""
        return {f"A1_{p.panel_id}": p.seed_id for p in self.panel}


def extract_model_member(models_tar: Path, genome_id: str, suffix: str, dest_dir: Path) -> Path:
    """Extract `<genome>/<genome><suffix>` from the benchmark models tar into dest_dir and return its path."""
    member_name = f"{genome_id}/{genome_id}{suffix}"
    dest = Path(dest_dir) / f"{genome_id}{suffix}"
    if dest.exists():
        return dest
    with tarfile.open(models_tar, "r") as tf:
        member = tf.getmember(member_name)
        fh = tf.extractfile(member)
        if fh is None:
            raise FileNotFoundError(member_name)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(fh.read())
    return dest


@dataclass
class GenomeInputs:
    genome_id: str
    layout: Optional[GenomeLayout]
    candidates: Tuple[CandidateClaim, ...]
    transporters: Tuple[TransporterHit, ...]
    transporter_proteins_by_cpd: Dict[str, Set[str]]

    @classmethod
    def load(
        cls,
        genome_id: str,
        layout_tsv: Optional[Path],
        candidates_tsv: Path,
        transporter_tbl: Optional[Path],
    ) -> "GenomeInputs":
        layout = GenomeLayout.from_tsv(layout_tsv) if layout_tsv is not None and Path(layout_tsv).exists() else None
        transporters = load_gapsmith_transporter_tbl(transporter_tbl) if transporter_tbl is not None else ()
        return cls(
            genome_id=genome_id,
            layout=layout,
            candidates=load_candidates(candidates_tsv),
            transporters=tuple(transporters),
            transporter_proteins_by_cpd=transporter_substrates(transporters),
        )
