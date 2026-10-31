"""Evidence categories shared by the v3 training-data tools.

Three distinct things were previously folded into one "experimental" set:

* ``EXPERIMENTAL`` - evidence from an experiment (GO: EXP/IDA/IPI/IMP/IGI/IEP and the
  high-throughput H* codes; Swiss-Prot EC/GO text: ECO:0000269).
* ``CURATED_INFERENTIAL`` - curator-asserted but not experimental (GO: TAS, IC; ECO:0000303
  author statement without traceable support, ECO:0000305 curator inference).
* everything else (IEA, ISS, ISO, ND, similarity/automatic ECO tags, untagged) - not truth.

``CAFA_ELIGIBLE`` = EXPERIMENTAL + CURATED_INFERENTIAL is the CAFA III/IV ground-truth code set;
it must be named as such, never called "experimental". The strict ``experimental`` policy is the
default everywhere.
"""
from __future__ import annotations

EXPERIMENTAL_GO = frozenset({"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"})
CURATED_INFERENTIAL_GO = frozenset({"TAS", "IC"})
EXPERIMENTAL_ECO = frozenset({"ECO:0000269"})
CURATED_INFERENTIAL_ECO = frozenset({"ECO:0000303", "ECO:0000305"})

EXPERIMENTAL = EXPERIMENTAL_GO | EXPERIMENTAL_ECO
CURATED_INFERENTIAL = CURATED_INFERENTIAL_GO | CURATED_INFERENTIAL_ECO
CAFA_ELIGIBLE = EXPERIMENTAL | CURATED_INFERENTIAL
POLICIES = ("experimental", "cafa_eligible", "all")


def category(code: str) -> str:
    """``experimental`` | ``curated_inferential`` | ``other`` for a GO evidence code or ECO id."""
    code = (code or "").strip()
    if code in EXPERIMENTAL:
        return "experimental"
    if code in CURATED_INFERENTIAL:
        return "curated_inferential"
    return "other"


def allowed(code: str, policy) -> bool:
    """``policy``: 'experimental' (strict), 'cafa_eligible', 'all', or an explicit list of codes."""
    if policy == "all":
        return True
    if policy == "experimental":
        return (code or "").strip() in EXPERIMENTAL
    if policy == "cafa_eligible":
        return (code or "").strip() in CAFA_ELIGIBLE
    if isinstance(policy, (list, tuple, set, frozenset)):
        return (code or "").strip() in set(policy)
    raise ValueError(f"unknown evidence policy {policy!r}")
