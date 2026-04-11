"""
Synthetic perturbations for the Phase 7.4 benchmark:

  ablate_claims_for_gapfill
      Remove 10% (configurable) of true GO annotations from the claims
      file. The claims file is the raw input to gspa-cli integrate; the
      removed (protein, aspect, go_term) pairs are returned as the
      ablated set so the evaluator can compute how many were recovered
      by the integrator's priors.

  inject_taxon_violations
      Sample from a curated list of "implausible for this organism"
      GO terms (e.g., photosynthesis in a heterotroph), add them to the
      claims file with calibrated weak-to-moderate raw scores drawn from
      the empirical low-confidence distribution, and return the
      (acc, aspect, go_term) set so the evaluator can score suppression.
"""

from __future__ import annotations

import json
import os
import random
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Set, Tuple


# ---------------------------------------------------------------------------
# Shared I/O
# ---------------------------------------------------------------------------


def _read_claims(path: str) -> List[dict]:
    out: List[dict] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            out.append(json.loads(line))
    return out


def _write_claims(path: str, claims: Iterable[dict]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or '.', exist_ok=True)
    with open(path, 'w') as fh:
        for c in claims:
            fh.write(json.dumps(c, separators=(',', ':')))
            fh.write('\n')


# ---------------------------------------------------------------------------
# Ablation (gap-fill test)
# ---------------------------------------------------------------------------


@dataclass
class AblationResult:
    ablated_claims_path: str
    removed_pairs: Dict[str, Set[Tuple[str, str]]] = field(default_factory=dict)
    n_removed: int = 0

    def as_dict(self) -> dict:
        return {
            'ablated_claims_path': self.ablated_claims_path,
            'n_removed': self.n_removed,
            'removed_pairs': {
                acc: sorted(list(pairs)) for acc, pairs in self.removed_pairs.items()
            },
        }


def ablate_claims_for_gapfill(
    claims_path: str,
    truth: Mapping[str, Mapping[str, Iterable[str]]],
    output_path: str,
    fraction: float = 0.1,
    seed: int = 42,
) -> AblationResult:
    """
    Remove `fraction` of each protein's true (aspect, go_term) pairs from
    the claims file so that we can check whether the Phase 7 priors
    recover them. Claims whose (protein, go_term) is in the ablated set
    are dropped regardless of source or aspect (the pair is considered
    globally "missing" for that protein).

    Note: we don't remove supporting claims for descendant GO terms of
    the ablated term — the test asks whether the integrator can recover
    the exact removed term, not a close relative.
    """
    rng = random.Random(seed)
    removed_pairs: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    for acc, aspect_map in truth.items():
        per_protein: List[Tuple[str, str]] = [
            (aspect, term)
            for aspect, terms in aspect_map.items()
            for term in terms
        ]
        if not per_protein:
            continue
        k = max(1, int(round(len(per_protein) * fraction)))
        removed = set(rng.sample(per_protein, k=min(k, len(per_protein))))
        removed_pairs[acc] = removed

    # Filter claims
    claims = _read_claims(claims_path)
    kept: List[dict] = []
    n_dropped = 0
    for c in claims:
        acc = c.get('protein_id')
        aspect = c.get('go_aspect')
        func = c.get('function_id')
        if acc in removed_pairs and (aspect, func) in removed_pairs[acc]:
            n_dropped += 1
            continue
        kept.append(c)
    _write_claims(output_path, kept)
    return AblationResult(
        ablated_claims_path=output_path,
        removed_pairs=removed_pairs,
        n_removed=n_dropped,
    )


# ---------------------------------------------------------------------------
# Taxon violation injection
# ---------------------------------------------------------------------------


# A curated list of "implausible for most heterotrophic bacteria" GO terms
# used as taxon-violation probes. These are chosen to trigger
# only_in_taxon / never_in_taxon constraints from go-plus.
DEFAULT_IMPLAUSIBLE_GO_TERMS: List[Tuple[str, str]] = [
    ('BP', 'GO:0015979'),   # photosynthesis
    ('BP', 'GO:0009773'),   # photosynthetic electron transport in photosystem I
    ('CC', 'GO:0009523'),   # photosystem II
    ('CC', 'GO:0009654'),   # photosystem II oxygen-evolving complex
    ('BP', 'GO:0042710'),   # biofilm formation (flagged for many strict anaerobes)
    ('BP', 'GO:0006952'),   # defense response (eukaryote-biased)
    ('CC', 'GO:0005737'),   # cytoplasm — valid GO but eukaryote-preferring
    ('BP', 'GO:0009835'),   # fruit ripening
    ('BP', 'GO:0048856'),   # anatomical structure development
]


@dataclass
class InjectionResult:
    injected_claims_path: str
    injected_pairs: Dict[str, Set[Tuple[str, str]]] = field(default_factory=dict)
    input_scores: Dict[str, Dict[Tuple[str, str], float]] = field(default_factory=dict)
    n_injected: int = 0

    def as_dict(self) -> dict:
        return {
            'injected_claims_path': self.injected_claims_path,
            'n_injected': self.n_injected,
            'injected_pairs': {
                acc: sorted(list(pairs)) for acc, pairs in self.injected_pairs.items()
            },
        }


def inject_taxon_violations(
    claims_path: str,
    truth: Mapping[str, Mapping[str, Iterable[str]]],
    output_path: str,
    n_injections: int = 50,
    raw_score_range: Tuple[float, float] = (0.3, 0.65),
    source: str = 'diamond',
    implausible_terms: Optional[List[Tuple[str, str]]] = None,
    seed: int = 4242,
) -> InjectionResult:
    """
    Inject n_injections taxon-violating GO claims into the claims file.
    The target proteins are randomly sampled from the truth set; the
    terms are drawn from a curated implausible-for-bacteria list by
    default.
    """
    if implausible_terms is None:
        implausible_terms = DEFAULT_IMPLAUSIBLE_GO_TERMS

    rng = random.Random(seed)
    proteins = list(truth.keys())
    if not proteins or not implausible_terms:
        return InjectionResult(injected_claims_path=output_path)

    claims = _read_claims(claims_path)
    injected_pairs: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    input_scores: Dict[str, Dict[Tuple[str, str], float]] = defaultdict(dict)
    n_injected = 0

    for _ in range(n_injections):
        acc = rng.choice(proteins)
        aspect, term = rng.choice(implausible_terms)
        if (aspect, term) in injected_pairs[acc]:
            continue
        raw_score = rng.uniform(*raw_score_range)
        claims.append({
            'protein_id': acc,
            'function_type': 'GO',
            'function_id': term,
            'go_aspect': aspect,
            'source': source,
            'raw_score': raw_score,
            'metadata': {'injected_taxon_violation': True},
        })
        injected_pairs[acc].add((aspect, term))
        input_scores[acc][(aspect, term)] = raw_score
        n_injected += 1

    _write_claims(output_path, claims)
    return InjectionResult(
        injected_claims_path=output_path,
        injected_pairs=injected_pairs,
        input_scores=input_scores,
        n_injected=n_injected,
    )
