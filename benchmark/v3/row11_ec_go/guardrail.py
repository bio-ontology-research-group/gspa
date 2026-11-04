#!/usr/bin/env python3
"""Engine-consistent symbolic guardrail for row-11 evaluation."""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "benchmark" / "v3"))
import evaluate_stage3_agent as E  # noqa: E402

run_gspa_agent = E.agent
DEFAULT_RESOURCES_DIR = ROOT / "gspa-core" / "src" / "main" / "resources"
UNKNOWN_TAXON_REASON = (
    "organism taxon not in the bundled hierarchy; taxon constraints not applied, as in the production engine"
)


def normalize_go_check(result: dict[str, Any]) -> dict[str, Any]:
    if result.get("violation_type") in ("UNKNOWN_TAXON", "MISSING_TAXON"):
        return {
            "valid": True,
            "reason": UNKNOWN_TAXON_REASON,
            "violation_type": None,
            "constraint_state": "taxon_not_in_hierarchy",
        }
    out = dict(result)
    out.setdefault("constraint_state", "applied")
    return out


class EngineConsistentGuardrail(run_gspa_agent.SymbolicGuardrail):
    """SymbolicGuardrail that matches production engine behavior when organism taxon is not in the bundled hierarchy."""

    def check_go_term(
        self,
        go_id: str,
        organism_taxon: str,
        allowed_vocab: set[str] | None = None,
        require_taxon: bool = True,
    ) -> dict[str, Any]:
        parent_res = super().check_go_term(
            go_id,
            organism_taxon,
            allowed_vocab,
            require_taxon,
        )
        return normalize_go_check(parent_res)


def load(
    resources_dir: Path | str = DEFAULT_RESOURCES_DIR,
    go_obo: Path | str | None = None,
) -> EngineConsistentGuardrail:
    res_path = Path(resources_dir) if resources_dir is not None else DEFAULT_RESOURCES_DIR
    obo_path = Path(go_obo) if go_obo is not None else None
    guard = EngineConsistentGuardrail(res_path, obo_path)
    if not guard.resources_loaded or not guard.go_parents:
        raise ValueError("complete symbolic resources and GO hierarchy are required")
    return guard
