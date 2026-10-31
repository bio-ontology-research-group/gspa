#!/usr/bin/env python3
"""Train and export the learned context-integration model (GSPA v3, items 4/7).

Reads a feature table whose columns follow ``gspa.context-features/1`` (the contract JSON
written by ``FeatureContract.toJson()``), fits an interpretable regularised logistic
regression and/or shallow gradient-boosted trees, calibrates on a SEPARATE split, and writes
a frozen JSON artifact that ``gspa.integration.learned.LearnedContextModel`` loads.

Protocol (every step is enforced in code, not by convention)
-----------------------------------------------------------
* Input rows MUST declare ``example_id``, ``genome_id``, ``group`` (family / species cluster)
  and ``split`` in {train, calibration, test}. Missing provenance is rejected; nothing is
  ever labelled "out-of-fold" on trust.
* The splits must be disjoint on the declared keys (default ``group`` and ``genome_id``).
* ``train()`` receives train + calibration rows ONLY. Test rows are touched by
  ``evaluate_frozen()`` after the artifact has been serialised and hashed, so test labels
  cannot influence weights, normalisation, hyper-parameters (chosen by group-wise CV inside
  the train split), the Platt calibration or the display threshold (chosen on calibration).
* Sample weights (``--weighting group`` by default) stop large genomes / families dominating.
* The artifact carries the source file hash, split identity (ids + groups hash per split),
  feature list + contract fingerprint, normalisation, weights/trees, calibration, threshold,
  hyper-parameters and ``frozen_before_test: true``. It has no timestamps, so identical
  inputs give byte-identical artifacts.
* The optional locus layer (explicit "no supplied candidate is correct" outcome) is trained
  only if the TRAIN split has enough fully-negative groups and enough single-positive groups;
  otherwise the artifact states ``trained: false`` with the reason and the JVM refuses to
  emit any null probability. There is no default null weight.

Nothing here is evidence of accuracy by itself: results on synthetic data only verify the
mechanism. Real claims need the datasets listed in docs/concepts/learned-context.md.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import optimize
from scipy.stats import rankdata

TRAINER_VERSION = "4"
ARTIFACT_SCHEMA = "gspa.learned-context-model/1"
LOCUS_GROUP_ID = "gspa.locus-group-features/1"
SPLITS = ("train", "calibration", "test")
LOCUS_GROUP_FEATURES = [
    "grp.max_margin", "grp.second_margin_available", "grp.second_margin",
    "grp.log1p_n_candidates", "grp.any_seq_support", "grp.asm_fragmented",
]


class ProvenanceError(ValueError):
    """The input violates the split / provenance protocol; training must not proceed."""


# ----------------------------------------------------------------------------- contract / data

def contract_fingerprint(cid: str, features: list[dict]) -> str:
    md = hashlib.sha256()
    md.update(cid.encode("utf-8"))
    for f in features:
        md.update(("\n" + f["name"] + "\t" + f["group"]).encode("utf-8"))
    return md.hexdigest()


def locus_group_fingerprint() -> str:
    # Mirrors FeatureContract.LOCUS_GROUP (all columns are group CONTEXT_NEIGHBORHOOD).
    return contract_fingerprint(LOCUS_GROUP_ID, [{"name": n, "group": "CONTEXT_NEIGHBORHOOD"} for n in LOCUS_GROUP_FEATURES])


def load_contract(path: Path) -> dict:
    c = json.loads(Path(path).read_text())
    fp = contract_fingerprint(c["id"], c["features"])
    if c.get("fingerprint") not in (None, fp):
        raise ProvenanceError(f"contract {c['id']}: stated fingerprint does not match its feature list")
    c["fingerprint"] = fp
    c["names"] = [f["name"] for f in c["features"]]
    return c


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_table(path: Path, contract: dict) -> list[dict]:
    """Read JSONL (``features`` object) or TSV (one column per feature) into validated rows."""
    path = Path(path)
    names = contract["names"]
    raw: list[dict] = []
    if path.suffix == ".jsonl":
        with open(path) as fh:
            for ln, line in enumerate(fh, 1):
                if line.strip():
                    try:
                        raw.append(json.loads(line))
                    except json.JSONDecodeError as e:
                        raise ProvenanceError(f"{path}:{ln}: invalid JSON ({e})")
    elif path.suffix in (".tsv", ".txt"):
        with open(path, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                row = {k: rec.get(k) for k in ("example_id", "genome_id", "group", "split", "label",
                                              "locus_group_id", "role", "function_type", "function_id", "candidate_key",
                                              "label_basis", "label_references", "locus_members", "locus_truth_source", "sequence_sources", "sequence_contracts", "quality_evidence")}
                for key in ('label_references', 'locus_members', 'sequence_sources', 'sequence_contracts', 'quality_evidence'):
                    if row.get(key):
                        row[key] = json.loads(row[key])
                row["features"] = {n: rec.get(n) for n in names}
                if 'cohort_context_mode' in rec:
                    row['cohort_context_mode'] = rec['cohort_context_mode']
                raw.append(row)
    else:
        raise ProvenanceError(f"unsupported table type {path.suffix}; use .jsonl or .tsv")
    rows = []
    seen = set()
    for i, r in enumerate(raw, 1):
        where = f"{path.name} row {i}"
        if r.get('quality_evidence'):
            raise ProvenanceError(f"{where}: gene-model quality observations are not supported by the numeric feature contract; do not discard them to train or score this row")
        for req in ("example_id", "genome_id", "group", "split", "label", "features"):
            if r.get(req) in (None, ""):
                raise ProvenanceError(f"{where}: missing required field '{req}' (provenance is mandatory)")
        if r["split"] not in SPLITS:
            raise ProvenanceError(f"{where}: split must be one of {SPLITS}, got {r['split']!r}")
        try:
            label = int(r["label"])
        except (TypeError, ValueError):
            raise ProvenanceError(f"{where}: label must be 0/1, got {r['label']!r}")
        if isinstance(r['label'], bool) or r['label'] not in (0, 1, '0', '1'):
            raise ProvenanceError(f"{where}: label must be 0/1, got {r['label']!r}")
        eid = str(r["example_id"])
        if eid in seen:
            raise ProvenanceError(f"{where}: duplicate example_id {eid!r}")
        seen.add(eid)
        feats = r["features"]
        missing = [n for n in names if n not in feats]
        extra = [k for k in feats if k not in set(names)]
        if missing or extra:
            raise ProvenanceError(f"{where}: features do not match contract {contract['id']} "
                                  f"(missing {missing[:3]}, unexpected {extra[:3]})")
        try:
            x = np.array([float(feats[n]) for n in names], dtype=float)
        except (TypeError, ValueError):
            raise ProvenanceError(f"{where}: non-numeric feature value")
        if not np.all(np.isfinite(x)):
            raise ProvenanceError(f"{where}: non-finite feature value (encode unavailable as 0 + its availability flag)")
        lg = r.get("locus_group_id")
        rows.append(dict(example_id=eid, genome_id=str(r["genome_id"]), group=str(r["group"]), split=r["split"],
                         label=label, x=x, locus_group_id=(str(lg) if lg not in (None, "") else None),
                         role=r.get("role"), **{k: r[k] for k in ('function_type', 'function_id', 'candidate_key',
                         'label_basis', 'label_references', 'locus_members', 'locus_truth_source', 'sequence_sources',
                         'cohort_context_mode', 'sequence_contracts') if k in r}))
        validate_context_mode(r.get('cohort_context_mode', 'UNIVERSAL'))
        if 'sequence_sources' in r and r['sequence_sources'] is not None:
            validate_sequence_sources(r['sequence_sources'])
        if r.get('sequence_contracts') is not None:
            validate_sequence_contracts(r['sequence_contracts'], r.get('sequence_sources'))
    validate_locus_truth(rows, required=False)
    return rows


def validate_sequence_sources(sources):
    if (not isinstance(sources, list) or any(not isinstance(s, str) or not s or s != s.strip().lower() for s in sources)
            or sources != sorted(set(sources))):
        raise ProvenanceError('sequence_sources must be sorted unique normalized source names')


def support_features(names):
    return [n for n in names if n.startswith(('fn.type_', 'fn.go_aspect_', 'role.'))
            or n.endswith('.reported') or n.endswith('.available') or n.endswith('_available')]


def score_support_features(names):
    return [n for n in names if n.startswith('seq.') and n not in support_features(names)]


def validate_sequence_contracts(entries, sources):
    if entries is None:
        return False
    keys = {'schema', 'tool', 'tool_version', 'tool_sha256', 'reference_sha256',
            'parameters_sha256', 'adapter_sha256', 'score_model'}
    if not isinstance(entries, list) or sources is None:
        raise ProvenanceError('invalid sequence_contracts identity')
    validate_sequence_sources(sources)
    for entry in entries:
        if not isinstance(entry, dict) or set(entry) != {'source', 'contract'} or entry['source'] not in sources:
            raise ProvenanceError('invalid sequence_contracts source')
        c = entry['contract']
        if c is not None and (not isinstance(c, dict) or set(c) != keys or c['schema'] != 'gspa.predictor-contract/1'
                or any(not isinstance(v, str) or not v or v != v.strip() for v in c.values())
                or any(not re.fullmatch('[0-9a-f]{64}', c[k]) for k in keys if k.endswith('_sha256'))):
            raise ProvenanceError('invalid predictor contract')
    canonical = [json.dumps(e, sort_keys=True, separators=(',', ':'), ensure_ascii=True) for e in entries]
    if canonical != sorted(set(canonical)) or sorted({e['source'] for e in entries}) != sources:
        raise ProvenanceError('sequence_contracts must be sorted unique and cover the sources')
    return all(e['contract'] is not None for e in entries)


def validate_context_mode(mode):
    if mode not in ('UNIVERSAL', 'PATHWAY', 'ORGANISM'):
        raise ProvenanceError('cohort_context_mode must be UNIVERSAL, PATHWAY or ORGANISM')
    return mode


def build_input_support(train_rows, calib_rows, names):
    """Necessary input-coverage checks, not proof of calibration/generalization.

    Profiles preserve actual source identities even when numerical columns pool
    them as 'other'. Only train/calibration evidence determines coverage.
    """
    if any(r.get('sequence_sources') is None for r in train_rows + calib_rows):
        return {'schema': 'gspa.model-input-support/3', 'available': False,
                'reason': 'training rows lack exact sequence_sources; re-export production features'}
    if any(not validate_sequence_contracts(r.get('sequence_contracts'), r['sequence_sources']) for r in train_rows + calib_rows):
        return {'schema': 'gspa.model-input-support/3', 'available': False,
                'reason': 'training rows lack complete predictor identities; re-export verified inputs'}
    flags = support_features(names)
    if not flags:
        return {'schema': 'gspa.model-input-support/3', 'available': False, 'reason': 'contract has no source/function/role indicators'}
    positions = [names.index(n) for n in flags]
    score_names = score_support_features(names)
    score_positions = [names.index(n) for n in score_names]
    profiles = {}
    for split, rows in (('train', train_rows), ('calibration', calib_rows)):
        for r in rows:
            sources = r['sequence_sources']
            validate_sequence_sources(sources)
            values = [float(r['x'][i]) for i in positions]
            if any(v not in (0, 1) for v in values):
                raise ProvenanceError('input-support indicators must be binary')
            values = [int(v) for v in values]
            mode = validate_context_mode(r.get('cohort_context_mode', 'UNIVERSAL'))
            contracts = r['sequence_contracts']
            identity = [sources, values, mode, contracts]
            key = json.dumps(identity, separators=(',', ':'), sort_keys=True)
            scores = [float(r['x'][i]) for i in score_positions]
            if any(not math.isfinite(v) for v in scores):
                raise ProvenanceError('sequence score ranges require finite features')
            p = profiles.setdefault(key, {'sources': sources, 'flags': values,
                'train_labels': [0, 0], 'calibration_labels': [0, 0], 'cohort_context_mode': mode,
                'sequence_contracts': contracts, 'lower': scores[:], 'upper': scores[:]})
            p['lower'] = [min(a, b) for a, b in zip(p['lower'], scores)]
            p['upper'] = [max(a, b) for a, b in zip(p['upper'], scores)]
            p[split + '_labels'][r['label']] += 1
    return json.loads(json.dumps({'schema': 'gspa.model-input-support/3', 'available': True,
            'features': flags, 'score_features': score_names, 'profiles': [profiles[k] for k in sorted(profiles)],
            'interpretation': 'Exact predictor/source/function/role/context identities with both classes in train and calibration; sequence features within their observed profile ranges. Necessary coverage only; no joint-distribution, taxonomic-transfer or calibration guarantee.'}))


def validate_contract_support(support, names):
    if support.get('features') != support_features(names) or support.get('score_features') != score_support_features(names):
        raise ProvenanceError('input-support feature lists do not match the contract')
    if not isinstance(support.get('profiles'), list):
        raise ProvenanceError('invalid input-support profiles')
    seen = set()
    for p in support['profiles']:
        if not isinstance(p, dict):
            raise ProvenanceError('invalid input-support profile')
        sources = p.get('sources')
        validate_sequence_sources(sources)
        if not validate_sequence_contracts(p.get('sequence_contracts'), sources):
            raise ProvenanceError('incomplete predictor profile')
        validate_context_mode(p.get('cohort_context_mode'))
        flags = p.get('flags')
        if not isinstance(flags, list) or len(flags) != len(support['features']) or any(type(v) not in (int, float) or v not in (0, 1) for v in flags):
            raise ProvenanceError('invalid profile flags')
        for k in ('train_labels', 'calibration_labels'):
            counts = p.get(k)
            if not isinstance(counts, list) or len(counts) != 2 or any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 or v != int(v) for v in counts):
                raise ProvenanceError('invalid profile class counts')
        for k in ('lower', 'upper'):
            bounds = p.get(k)
            if not isinstance(bounds, list) or len(bounds) != len(support['score_features']) or any(type(v) not in (int, float) or not math.isfinite(v) for v in bounds):
                raise ProvenanceError('invalid sequence score ranges')
        if any(lo > hi for lo, hi in zip(p['lower'], p['upper'])):
            raise ProvenanceError('invalid sequence score ranges')
        identity = json.dumps([sources, [int(v) for v in flags], p['cohort_context_mode'], p['sequence_contracts']], sort_keys=True)
        if identity in seen:
            raise ProvenanceError('duplicate input-support profile')
        seen.add(identity)


def unsupported_input_reason(art, x, sequence_sources=None, quality_evidence=None, cohort_context_mode='UNIVERSAL', sequence_contracts=None):
    validate_context_mode(cohort_context_mode)
    if quality_evidence:
        return 'gene_model_quality_not_supported_by_model'
    if art.get('intended_use') == 'mechanism_test':
        return None
    support = art.get('input_support')
    if support is not None and (not isinstance(support, dict) or type(support.get('available')) is not bool or
            support.get('schema') not in ('gspa.model-input-support/1', 'gspa.model-input-support/2', 'gspa.model-input-support/3')):
        raise ProvenanceError('invalid input-support schema/availability')
    if support is None or not support['available']:
        return 'model_input_support_unavailable'
    if support.get('schema') not in ('gspa.model-input-support/1', 'gspa.model-input-support/2', 'gspa.model-input-support/3'):
        raise ProvenanceError('invalid input-support schema')
    names = art['contract']['features']
    if len(x) != len(names) or not np.all(np.isfinite(x)):
        raise ProvenanceError('invalid feature vector')
    if support['schema'] == 'gspa.model-input-support/3':
        validate_contract_support(support, names)
    mode_aware = support['schema'] != 'gspa.model-input-support/1'
    modes = {validate_context_mode(p.get('cohort_context_mode')) for p in support['profiles']} if mode_aware else {'UNIVERSAL'}
    if cohort_context_mode not in modes:
        return 'unseen_reaction_context_mode'
    if sequence_sources is None:
        return 'candidate_source_identity_unavailable'
    validate_sequence_sources(sequence_sources)
    if '(unnamed)' in sequence_sources:
        return 'unnamed_sequence_source'
    if support['schema'] != 'gspa.model-input-support/3':
        return 'model_predictor_contract_support_unavailable'
    if not validate_sequence_contracts(sequence_contracts, sequence_sources):
        return 'candidate_predictor_contract_unavailable'
    values = [float(x[names.index(n)]) for n in support['features']]
    if any(v not in (0, 1) for v in values):
        raise ProvenanceError('input-support indicators must be binary')
    known_base = False
    for p in support['profiles']:
        if mode_aware and p['cohort_context_mode'] != cohort_context_mode:
            continue
        if p['sources'] == sequence_sources and p['flags'] == values:
            known_base = True
        if p['sources'] == sequence_sources and p['flags'] == values and p['sequence_contracts'] == sequence_contracts:
            if not all(v > 0 for v in p['train_labels'] + p['calibration_labels']):
                return 'profile_lacks_both_classes_in_train_and_calibration'
            score_values = [float(x[names.index(n)]) for n in support['score_features']]
            if any(v < lo or v > hi for v, lo, hi in zip(score_values, p['lower'], p['upper'])):
                return 'sequence_score_outside_training_calibration_range'
            return None
    return 'unseen_predictor_contract' if known_base else 'unseen_source_function_role_or_availability_profile'


def validate_locus_truth(rows, required=True):
    """Preserve complete declared membership through filtering/serialization.

    Production null supervision needs explicit negative assays, never unconfirmed
    closed-world rows. The source statement is auditable provenance, not proof of
    biological completeness by itself; the builder verifies the exported universe.
    """
    groups = {}
    for r in rows:
        if r.get('locus_group_id') is not None:
            groups.setdefault((r['genome_id'], r['locus_group_id']), []).append(r)
    for key, members in groups.items():
        actual = {r['example_id'] for r in members}
        if len(actual) != len(members):
            raise ProvenanceError(f'locus {key}: duplicate example identity')
        if len({(r['split'], r['group']) for r in members}) != 1:
            raise ProvenanceError(f'locus {key}: split or ancestry group differs among members')
        present = any(r.get('locus_members') is not None or r.get('locus_truth_source') is not None for r in members)
        if not required and not present:
            continue  # legacy mechanism fixtures only; production train() requires the certificate
        for r in members:
            expected = r.get('locus_members')
            if (not isinstance(expected, list) or not expected or not all(isinstance(x, str) for x in expected)
                    or len(expected) != len(set(expected)) or set(expected) != actual):
                raise ProvenanceError(f'locus {key}: missing or incomplete locus_members certificate')
            if not isinstance(r.get('locus_truth_source'), str) or not r['locus_truth_source'].strip():
                raise ProvenanceError(f'locus {key}: missing locus_truth_source')
        if len({r['locus_truth_source'] for r in members}) != 1:
            raise ProvenanceError(f'locus {key}: inconsistent truth source')
        if required:
            identity = {(r.get('function_type'), r.get('function_id'), r.get('role')) for r in members}
            if len(identity) != 1 or any(not v for v in next(iter(identity))):
                raise ProvenanceError(f'locus {key}: members need the same function/role identity')
            for r in members:
                expected_basis = 'positive_truth' if r['label'] else 'explicit_negative_truth'
                if r.get('label_basis') != expected_basis or (not r['label'] and not r.get('label_references')):
                    raise ProvenanceError(f'locus {key}: needs explicit positive/negative truth, not unconfirmed labels')


def check_disjoint(rows: list[dict], keys: list[str]) -> None:
    """Splits must not share any value of the declared keys."""
    for key in keys:
        owner: dict[str, str] = {}
        for r in rows:
            v = r[key]
            if v in owner and owner[v] != r["split"]:
                raise ProvenanceError(
                    f"{key} {v!r} appears in both '{owner[v]}' and '{r['split']}': splits must be "
                    f"disjoint by {key}; refusing to claim held-out evaluation")
            owner[v] = r["split"]
    present = {r["split"] for r in rows}
    for s in SPLITS:
        if s not in present:
            raise ProvenanceError(f"no rows declared for split '{s}'")


def split_identity(rows: list[dict], keys: list[str]) -> dict:
    out = {}
    for s in SPLITS:
        sub = [r for r in rows if r["split"] == s]
        ids = sorted(r["example_id"] for r in sub)
        out[s] = {
            "n": len(sub),
            "example_ids_sha256": hashlib.sha256("\n".join(ids).encode()).hexdigest(),
            **{f"n_{k}s": len({r[k] for r in sub}) for k in ("group", "genome_id")},
            **{f"{k}s_sha256": hashlib.sha256("\n".join(sorted({r[k] for r in sub})).encode()).hexdigest()
               for k in ("group", "genome_id")},
        }
    out["disjoint_keys_enforced"] = list(keys)
    return out


def training_data_digest(train_rows: list[dict], calib_rows: list[dict]) -> str:
    """Hash of the train + calibration content only (ids, provenance, labels, features).

    The artifact records this instead of the whole input file's hash: the file also holds the
    test split, and an artifact that depended on test labels (even through a hash) would not be
    provably independent of them. The whole-file hash goes into the test report.
    """
    md = hashlib.sha256()
    for r in sorted(train_rows + calib_rows, key=lambda r: r["example_id"]):
        md.update(json.dumps([r["example_id"], r["genome_id"], r["group"], r["split"], r["label"],
                              r["locus_group_id"], [float(v) for v in r["x"]]]).encode("utf-8"))
        md.update(b"\n")
        truth = {k: r[k] for k in ('locus_members', 'locus_truth_source', 'label_basis', 'label_references',
                                  'function_type', 'function_id', 'role') if r.get(k) is not None}
        if any(r.get(k) is not None for k in ('locus_members', 'locus_truth_source', 'label_basis', 'label_references')):
            md.update(json.dumps(truth, sort_keys=True).encode('utf-8'))
            md.update(b'\n')
        if r.get('sequence_sources') is not None:
            md.update(json.dumps(r['sequence_sources']).encode('utf-8'))
            md.update(b'\n')
        if r.get('sequence_contracts') is not None:
            validate_sequence_contracts(r['sequence_contracts'], r.get('sequence_sources'))
            md.update(json.dumps(r['sequence_contracts'], sort_keys=True).encode('utf-8'))
            md.update(b'\n')
        mode = validate_context_mode(r.get('cohort_context_mode', 'UNIVERSAL'))
        if mode != 'UNIVERSAL':  # Preserve historical universal-only training identities.
            md.update(json.dumps({'cohort_context_mode': mode}, sort_keys=True).encode('utf-8'))
            md.update(b'\n')
    return md.hexdigest()


def sample_weights(rows: list[dict], mode: str) -> np.ndarray:
    """Weights that stop replication dominating; normalised to mean 1 within the given rows."""
    if mode == "none":
        return np.ones(len(rows))
    keyf = {"genome": lambda r: r["genome_id"], "group": lambda r: r["group"],
            "genome_group": lambda r: (r["genome_id"], r["group"])}.get(mode)
    if keyf is None:
        raise ProvenanceError(f"unknown weighting {mode}")
    counts: dict = {}
    for r in rows:
        counts[keyf(r)] = counts.get(keyf(r), 0) + 1
    w = np.array([1.0 / counts[keyf(r)] for r in rows])
    return w * (len(w) / w.sum())


def stack(rows: list[dict]):
    X = np.vstack([r["x"] for r in rows])
    y = np.array([r["label"] for r in rows], dtype=float)
    return X, y


# ----------------------------------------------------------------------------- models

def _sigmoid(z):
    return 0.5 * (1.0 + np.tanh(0.5 * z))


class LogisticModel:
    kind = "logistic"

    def __init__(self, l2: float):
        self.l2 = l2

    def fit(self, X, y, w):
        self.mean = X.mean(axis=0)
        sd = X.std(axis=0)
        self.scale = np.where(sd < 1e-12, 1.0, sd)
        Z = (X - self.mean) / self.scale
        d = Z.shape[1]
        l2 = self.l2

        def fun(theta):
            b0, beta = theta[0], theta[1:]
            z = b0 + Z @ beta
            loss = float(np.sum(w * (np.logaddexp(0.0, z) - y * z)) + 0.5 * l2 * beta @ beta)
            r = w * (_sigmoid(z) - y)
            grad = np.concatenate(([r.sum()], Z.T @ r + l2 * beta))
            return loss, grad

        res = optimize.minimize(fun, np.zeros(d + 1), jac=True, method="L-BFGS-B",
                                options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-8})
        self.intercept = float(res.x[0])
        self.weights = res.x[1:].copy()
        return self

    def margin(self, X):
        return self.intercept + ((X - self.mean) / self.scale) @ self.weights

    def artifact_part(self):
        return {
            "kind": "logistic",
            "normalization": {"mean": self.mean.tolist(), "scale": self.scale.tolist()},
            "logistic": {"weights": self.weights.tolist(), "intercept": self.intercept},
        }

    def params(self):
        return {"l2": self.l2}


class GbdtModel:
    """Newton-boosted shallow regression trees on the weighted logistic loss (numpy only)."""
    kind = "gbdt"

    def __init__(self, n_trees: int, depth: int, lr: float, min_leaf: int = 5, lam: float = 1.0):
        self.n_trees, self.depth, self.lr, self.min_leaf, self.lam = n_trees, depth, lr, min_leaf, lam

    def fit(self, X, y, w):
        base = float(np.clip(np.sum(w * y) / np.sum(w), 1e-6, 1 - 1e-6))
        self.init = math.log(base / (1 - base))
        z = np.full(len(y), self.init)
        self.trees = []
        for _ in range(self.n_trees):
            p = _sigmoid(z)
            g = w * (p - y)
            h = np.maximum(w * p * (1 - p), 1e-12)
            tree = self._grow(X, g, h)
            self.trees.append(tree)
            z = z + self.lr * self._eval_tree(tree, X)
        return self

    def _grow(self, X, g, h):
        nodes = {"feature": [], "threshold": [], "left": [], "right": [], "value": []}

        def new_node():
            for k in nodes:
                nodes[k].append(0)
            return len(nodes["feature"]) - 1

        def leaf(i, idx):
            nodes["feature"][i] = -1
            nodes["threshold"][i] = 0.0
            nodes["left"][i] = 0
            nodes["right"][i] = 0
            nodes["value"][i] = float(-g[idx].sum() / (h[idx].sum() + self.lam))

        def split(i, idx, depth):
            if depth >= self.depth or len(idx) < 2 * self.min_leaf:
                return leaf(i, idx)
            G, H = g[idx].sum(), h[idx].sum()
            parent = G * G / (H + self.lam)
            best = (1e-9, None, None)
            for f in range(X.shape[1]):
                col = X[idx, f]
                order = np.argsort(col, kind="stable")
                cs = col[order]
                if cs[0] == cs[-1]:
                    continue
                GL = np.cumsum(g[idx][order])[:-1]
                HL = np.cumsum(h[idx][order])[:-1]
                nL = np.arange(1, len(idx))
                ok = (cs[1:] > cs[:-1]) & (nL >= self.min_leaf) & (len(idx) - nL >= self.min_leaf)
                if not ok.any():
                    continue
                gain = GL ** 2 / (HL + self.lam) + (G - GL) ** 2 / (H - HL + self.lam) - parent
                gain = np.where(ok, gain, -np.inf)
                k = int(np.argmax(gain))
                if gain[k] > best[0]:
                    best = (float(gain[k]), f, 0.5 * (cs[k] + cs[k + 1]), )
            if best[1] is None:
                return leaf(i, idx)
            _, f, thr = best
            mask = X[idx, f] <= thr
            nodes["feature"][i] = int(f)
            nodes["threshold"][i] = float(thr)
            nodes["value"][i] = 0.0
            li = new_node(); nodes["left"][i] = li
            ri = new_node(); nodes["right"][i] = ri
            split(li, idx[mask], depth + 1)
            split(ri, idx[~mask], depth + 1)

        root = new_node()
        split(root, np.arange(X.shape[0]), 0)
        return nodes

    @staticmethod
    def _eval_tree(t, X):
        out = np.empty(X.shape[0])
        feat, thr, left, right, val = (np.array(t[k]) for k in ("feature", "threshold", "left", "right", "value"))
        for i in range(X.shape[0]):
            n = 0
            while feat[n] >= 0:
                n = left[n] if X[i, feat[n]] <= thr[n] else right[n]
            out[i] = val[n]
        return out

    def margin(self, X):
        z = np.full(X.shape[0], self.init)
        for t in self.trees:
            z = z + self.lr * self._eval_tree(t, X)
        return z

    def artifact_part(self):
        return {"kind": "gbdt", "normalization": None,
                "gbdt": {"init": self.init, "learning_rate": self.lr, "trees": self.trees}}

    def params(self):
        return {"n_trees": self.n_trees, "depth": self.depth, "learning_rate": self.lr,
                "min_leaf": self.min_leaf, "l2_leaf": self.lam}


# ----------------------------------------------------------------------------- selection / calibration

def weighted_logloss(y, p, w=None):
    p = np.clip(p, 1e-15, 1 - 1e-15)
    ll = -(y * np.log(p) + (1 - y) * np.log(1 - p))
    return float(np.average(ll, weights=w))


def grouped_folds(groups: list[str], k: int, seed: int) -> list[int]:
    """Deterministic group -> fold assignment; a group is never split across folds."""
    uniq = sorted(set(groups), key=lambda g: hashlib.sha256(f"{seed}|{g}".encode()).hexdigest())
    fold_of = {g: i % k for i, g in enumerate(uniq)}
    return [fold_of[g] for g in groups]


def cv_margins(make_model, X, y, w, groups, k, seed):
    """Out-of-fold margins within the training split (group-wise folds)."""
    folds = np.array(grouped_folds(groups, k, seed))
    oof = np.zeros(len(y))
    for f in range(k):
        tr, te = folds != f, folds == f
        if te.sum() == 0 or tr.sum() == 0 or len(set(y[tr])) < 2:
            raise ProvenanceError("a CV training fold has a single class; use fewer folds or more groups")
        m = make_model().fit(X[tr], y[tr], w[tr])
        oof[te] = m.margin(X[te])
    return oof


def select_hyperparameters(candidates, X, y, w, groups, k, seed):
    """Pick the configuration with the lowest grouped-CV weighted log-loss on TRAIN only."""
    scored = []
    for params, make in candidates:
        oof = cv_margins(make, X, y, w, groups, k, seed)
        scored.append((weighted_logloss(y, _sigmoid(oof), w), params, make))
    scored.sort(key=lambda t: t[0])
    return scored[0], [(s, p) for s, p, _ in scored]


def fit_platt(margin, y, w):
    if len(set(y.tolist())) < 2:
        raise ProvenanceError("calibration split has a single class; cannot fit a calibration")

    def fun(ab):
        z = ab[0] * margin + ab[1]
        loss = float(np.sum(w * (np.logaddexp(0.0, z) - y * z)))
        r = w * (_sigmoid(z) - y)
        return loss, np.array([np.sum(r * margin), np.sum(r)])

    res = optimize.minimize(fun, np.array([1.0, 0.0]), jac=True, method="L-BFGS-B",
                            bounds=[(1e-3, 100.0), (-30.0, 30.0)], options={"ftol": 1e-14, "gtol": 1e-10})
    return float(res.x[0]), float(res.x[1])


def choose_threshold(p, y, rule):
    """Advisory display threshold from CALIBRATION predictions only."""
    # Candidate cut points are the unique predicted values; "predict positive" means p >= cut.
    order = np.argsort(-p, kind="stable")
    ps, ys = p[order], y[order]
    last = np.concatenate((ps[1:] != ps[:-1], [True]))          # last index of each unique value
    tp = np.cumsum(ys)[last]
    fp = np.cumsum(1.0 - ys)[last]
    cuts = ps[last]
    npos = float(ys.sum())
    if rule == "max_f1":
        f1 = np.where(tp > 0, 2 * tp / (2 * tp + fp + (npos - tp)), 0.0)
        k = int(np.argmax(f1))                                   # ties: the higher (more conservative) cut
        return {"rule": rule, "value": float(cuts[k]), "selected_on": "calibration", "advisory_only": True,
                "calibration_f1": float(f1[k])}
    if rule.startswith("precision_at:"):
        target = float(rule.split(":", 1)[1])
        prec = tp / (tp + fp)
        ok = np.nonzero(prec >= target)[0]
        if len(ok):
            return {"rule": rule, "value": float(cuts[ok[-1]]), "selected_on": "calibration",
                    "advisory_only": True, "achievable": True}
        return {"rule": rule, "value": 1.0, "selected_on": "calibration", "advisory_only": True,
                "achievable": False}
    raise ProvenanceError(f"unknown threshold rule {rule}")


# ----------------------------------------------------------------------------- locus (null) layer

def _group_features(margins, feats_rows, idx):
    """Mirror of LearnedContextModel.groupFeatures; ``feats_rows`` are candidate feature vectors."""
    ms = sorted(margins, reverse=True)
    g = np.zeros(len(LOCUS_GROUP_FEATURES))
    g[0] = ms[0]
    if len(ms) > 1:
        g[1], g[2] = 1.0, ms[1]
    g[3] = math.log1p(len(ms))
    g[4] = 1.0 if any(f[idx["seq.n_groups_reporting"]] > 0 for f in feats_rows) else 0.0
    g[5] = 1.0 if any(f[idx["asm.available"]] > 0 and f[idx["asm.n_contigs"]] > math.log(2.0) + 1e-9
                      for f in feats_rows) else 0.0
    return g


def build_locus_groups(rows, margins, names):
    """Group within a genome; identical local names in other genomes are independent."""
    idx = {n: i for i, n in enumerate(names)}
    by: dict[str, list[int]] = {}
    for i, r in enumerate(rows):
        if r["locus_group_id"] is not None:
            by.setdefault((r['genome_id'], r["locus_group_id"]), []).append(i)
    out, multi = [], 0
    for gid, members in by.items():
        pos = [j for j, i in enumerate(members) if rows[i]["label"] == 1]
        if len(pos) > 1:
            multi += 1       # several correct loci (isoenzymes): not expressible as one-of-n; excluded, not forced
            continue
        m = [margins[i] for i in members]
        f = [rows[i]["x"] for i in members]
        out.append(dict(gid=gid, margins=np.array(m), g=_group_features(m, f, idx),
                        target=(pos[0] + 1 if pos else 0)))
    return out, multi


def locus_nll(theta, groups, l2):
    s, b0, w = math.exp(theta[0]), theta[1], theta[2:]
    nll = 0.0
    for gr in groups:
        logits = np.concatenate(([b0 + float(w @ gr["g"])], s * gr["margins"]))
        mx = logits.max()
        lse = mx + math.log(np.exp(logits - mx).sum())
        nll += lse - logits[gr["target"]]
    return nll / max(len(groups), 1) + 0.5 * l2 * float(w @ w)


def fit_locus(groups, l2):
    theta0 = np.zeros(2 + len(LOCUS_GROUP_FEATURES))
    res = optimize.minimize(lambda t: locus_nll(t, groups, l2), theta0, method="L-BFGS-B",
                            bounds=[(-5, 5)] + [(-30, 30)] * (1 + len(LOCUS_GROUP_FEATURES)))
    return res.x


def locus_probabilities(theta, margins, g):
    s, b0, w = math.exp(theta[0]), theta[1], theta[2:]
    logits = np.concatenate(([b0 + float(w @ g)], s * np.asarray(margins)))
    e = np.exp(logits - logits.max())
    p = e / e.sum()
    return float(p[0]), p[1:]


# ----------------------------------------------------------------------------- train / freeze / evaluate

def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, indent=1, allow_nan=False) + "\n"


def default_candidates(kind: str):
    if kind == "logistic":
        return [({"l2": l2}, (lambda l2=l2: LogisticModel(l2))) for l2 in (0.1, 1.0, 10.0, 100.0)]
    if kind == "gbdt":
        return [({"n_trees": 40, "depth": d, "learning_rate": 0.15},
                 (lambda d=d: GbdtModel(40, d, 0.15))) for d in (2, 3)]
    raise ProvenanceError(f"unknown model kind {kind}")


def train(train_rows, calib_rows, contract, kind, *, weighting="group", folds=4, seed=0,
          threshold_rule="max_f1", min_null_groups=20, min_pos_groups=20,
          identity, intended_use="production", synthetic=False):
    """Fit on train, choose hyper-parameters by grouped CV on train, calibrate on calibration.

    Deliberately has no access to test rows. Returns the artifact dict (not yet serialised).
    """
    if any(r.get('quality_evidence') for r in train_rows + calib_rows):
        raise ProvenanceError('gene-model quality observations require a supported feature contract before training')
    validate_locus_truth(train_rows + calib_rows, required=not synthetic)
    names = contract["names"]
    data_sha = training_data_digest(train_rows, calib_rows)
    Xtr, ytr = stack(train_rows)
    Xca, yca = stack(calib_rows)
    wtr = sample_weights(train_rows, weighting)
    wca = sample_weights(calib_rows, weighting)
    gtr = [r["group"] for r in train_rows]

    (cv_loss, params, make), all_scores = select_hyperparameters(default_candidates(kind), Xtr, ytr, wtr, gtr, folds, seed)
    model = make().fit(Xtr, ytr, wtr)

    cal_margin = model.margin(Xca)
    a, b = fit_platt(cal_margin, yca, wca)
    p_cal = _sigmoid(a * cal_margin + b)
    threshold = choose_threshold(p_cal, yca, threshold_rule)

    # ---- locus layer (explicit null outcome), trained on OOF margins of the train split
    locus = {"trained": False, "reason": "not attempted"}
    tr_has_loc = [r for r in train_rows if r["locus_group_id"] is not None]
    if tr_has_loc:
        oof = cv_margins(make, Xtr, ytr, wtr, gtr, folds, seed)
        groups, multi = build_locus_groups(train_rows, oof, names)
        n_null = sum(1 for g in groups if g["target"] == 0)
        n_pos = len(groups) - n_null
        cal_groups, _ = build_locus_groups(calib_rows, cal_margin, names)
        cn_null = sum(1 for g in cal_groups if g["target"] == 0)
        cn_pos = len(cal_groups) - cn_null
        if n_null < min_null_groups or n_pos < min_pos_groups:
            locus = {"trained": False, "reason":
                     f"train split has {n_null} no-correct-candidate groups and {n_pos} single-positive groups; "
                     f"need >= {min_null_groups} and >= {min_pos_groups} to learn an honest null outcome"}
        elif cn_null < 5 or cn_pos < 5:
            locus = {"trained": False, "reason":
                     f"calibration split has {cn_null} null and {cn_pos} positive groups (< 5 each): cannot validate the null layer"}
        else:
            best = None
            for l2 in (0.1, 1.0, 10.0):
                th = fit_locus(groups, l2)
                score = locus_nll(th, cal_groups, 0.0)
                if best is None or score < best[0]:
                    best = (score, l2, th)
            _, l2, th = best
            locus = {"trained": True,
                     "contract": {"id": LOCUS_GROUP_ID, "fingerprint": locus_group_fingerprint(),
                                  "features": LOCUS_GROUP_FEATURES},
                     "scale": float(math.exp(th[0])), "null_bias": float(th[1]),
                     "null_weights": [float(v) for v in th[2:]],
                     "hyperparameters": {"l2": l2, "selected_on": "calibration_nll"},
                     "n_train_groups": len(groups), "n_train_null_groups": n_null,
                     "n_multi_positive_groups_excluded": multi}
    else:
        locus = {"trained": False, "reason": "no rows carry locus_group_id; no grouped null-outcome data"}

    artifact = {
        "schema": ARTIFACT_SCHEMA,
        "model_id": f"context-{kind}-{data_sha[:12]}-{contract['fingerprint'][:8]}",
        "intended_use": intended_use,
        "contract": {"id": contract["id"], "fingerprint": contract["fingerprint"], "features": names},
        **model.artifact_part(),
        "calibration": {"method": "platt", "a": a, "b": b, "fitted_on": "calibration",
                        "n": int(len(yca)), "logloss_before": weighted_logloss(yca, _sigmoid(cal_margin), wca),
                        "logloss_after": weighted_logloss(yca, p_cal, wca)},
        "display_threshold": threshold,
        "locus_model": locus,
        "input_support": build_input_support(train_rows, calib_rows, names),
        "training": {
            "trainer_version": TRAINER_VERSION,
            "frozen_before_test": True,
            "synthetic": bool(synthetic),
            "training_data_sha256": data_sha,        # train + calibration rows only; never the test split
            "split_identity": identity,
            "weighting": weighting,
            "cv": {"folds": folds, "seed": seed, "unit": "group", "logloss": cv_loss,
                   "candidates": [{"params": p, "cv_logloss": s} for s, p in all_scores]},
            "hyperparameters": params,
            "n_train": int(len(ytr)), "n_calibration": int(len(yca)),
            "positive_rate_train": float(ytr.mean()),
            "truth": {
                "label_basis_counts": {s: dict(sorted(Counter(r.get('label_basis', 'unspecified') for r in rs).items()))
                                       for s, rs in (('train', train_rows), ('calibration', calib_rows))},
                "probability_target": ('synthetic_mechanism' if synthetic else
                    'reference_annotation_support_under_declared_closed_world_policy'
                    if any(r.get('label_basis') == 'closed_world_unconfirmed' for r in train_rows + calib_rows)
                    else 'supplied_explicit_binary_truth'
                    if all(r.get('label_basis') in ('positive_truth', 'explicit_negative_truth') for r in train_rows + calib_rows)
                    else 'unspecified_binary_label_target'),
                "locus_supervision": 'complete per-genome function/role groups with explicit truth certificates' if not synthetic else 'synthetic mechanism groups',
            },
        },
    }
    return artifact


def write_frozen(artifact, path: Path) -> str:
    text = canonical_json(artifact)
    Path(path).write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


# ---- reference predictor (the JVM class must agree with this to ~1e-12)

def load_artifact(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def predict_margin(art: dict, x: np.ndarray) -> float:
    if art["kind"] == "logistic":
        n = art["normalization"]
        xs = (x - np.array(n["mean"])) / np.array(n["scale"])
        return float(art["logistic"]["intercept"] + xs @ np.array(art["logistic"]["weights"]))
    g = art["gbdt"]
    z = g["init"]
    for t in g["trees"]:
        k = 0
        while t["feature"][k] >= 0:
            k = t["left"][k] if x[t["feature"][k]] <= t["threshold"][k] else t["right"][k]
        z += g["learning_rate"] * t["value"][k]
    return float(z)


def predict_proba(art: dict, x: np.ndarray, sequence_sources=None, quality_evidence=None, cohort_context_mode='UNIVERSAL', sequence_contracts=None) -> float:
    reason = unsupported_input_reason(art, x, sequence_sources, quality_evidence, cohort_context_mode, sequence_contracts)
    if reason:
        raise ProvenanceError('unsupported learned probability: ' + reason)
    c = art["calibration"]
    z = c["a"] * predict_margin(art, x) + c["b"]
    return float(_sigmoid(np.array(z)))


def assign_locus(art: dict, candidate_vectors: list[np.ndarray], sequence_sources=None, quality_evidence=None, cohort_context_modes=None, sequence_contracts=None):
    for i, x in enumerate(candidate_vectors):
        reason = unsupported_input_reason(art, x, sequence_sources[i] if sequence_sources is not None else None,
                                          quality_evidence[i] if quality_evidence is not None else None,
                                          cohort_context_modes[i] if cohort_context_modes is not None else 'UNIVERSAL',
                                          sequence_contracts[i] if sequence_contracts is not None else None)
        if reason:
            raise ProvenanceError('unsupported locus member: ' + reason)
    lm = art["locus_model"]
    if not lm.get("trained"):
        raise ValueError("no trained null-outcome layer: " + lm.get("reason", ""))
    names = art["contract"]["features"]
    idx = {n: i for i, n in enumerate(names)}
    margins = [predict_margin(art, x) for x in candidate_vectors]
    g = _group_features(margins, candidate_vectors, idx)
    theta = np.concatenate(([math.log(lm["scale"]), lm["null_bias"]], lm["null_weights"]))
    p0, pc = locus_probabilities(theta, margins, g)
    return p0, [float(v) for v in pc]


def evaluate_frozen(artifact_path: Path, contract: dict, test_rows, model_sha256: str,
                    source_file_sha256: str | None = None) -> dict:
    """Test-split metrics of the already-serialised artifact. Called once, after freezing."""
    art = load_artifact(artifact_path)
    if hashlib.sha256(Path(artifact_path).read_bytes()).hexdigest() != model_sha256:
        raise ProvenanceError("artifact changed between freezing and evaluation")
    validate_locus_truth(test_rows, required=not art['training'].get('synthetic', False))
    original_n = len(test_rows)
    reasons = [unsupported_input_reason(art, r['x'], r.get('sequence_sources'), r.get('quality_evidence'),
                                       r.get('cohort_context_mode', 'UNIVERSAL'), r.get('sequence_contracts')) for r in test_rows]
    unsupported_groups = {(r['genome_id'], r['locus_group_id']) for r, reason in zip(test_rows, reasons)
                          if reason and r.get('locus_group_id') is not None}
    coverage = {'total_rows': original_n, 'supported_rows': reasons.count(None),
                'unsupported_reasons': dict(Counter(r for r in reasons if r)),
                'locus_groups_with_unsupported_members': len(unsupported_groups)}
    all_positive_rows = sum(r['label'] for r in test_rows)
    test_rows = [dict(r, locus_group_id=None) if (r['genome_id'], r['locus_group_id']) in unsupported_groups else r
                 for r, reason in zip(test_rows, reasons) if reason is None]
    if not test_rows:
        return {'model_sha256': model_sha256, 'source_file_sha256': source_file_sha256,
                'evaluated_on': 'test after freezing', 'input_coverage': coverage, 'n': 0,
                'status': 'no_supported_test_rows', 'auroc': None, 'logloss': None, 'brier': None,
                'ece_10bin': None, 'locus': {'trained': False, 'reason': 'no supported test rows'}}
    X, y = stack(test_rows)
    margins = np.array([predict_margin(art, x) for x in X])
    c = art["calibration"]
    p = _sigmoid(c["a"] * margins + c["b"])
    thr = art["display_threshold"]["value"]
    pred = p >= thr
    tp = float(np.sum(pred & (y == 1)))
    ranks = rankdata(p)
    npos, nneg = float(y.sum()), float(len(y) - y.sum())
    auroc = float((ranks[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)) if npos and nneg else None
    bins = np.minimum((p * 10).astype(int), 9)
    ece = float(sum(abs(p[bins == b].mean() - y[bins == b].mean()) * np.mean(bins == b)
                    for b in range(10) if np.any(bins == b)))
    report = {
        "model_sha256": model_sha256,
        "input_coverage": coverage,
        "metric_scope": 'supported rows only; input_coverage reports abstentions over the full test set',
        "source_file_sha256": source_file_sha256,
        "evaluated_on": "test (labels used only now, after the artifact was frozen)",
        "probability_target": art['training'].get('truth', {}).get('probability_target', 'unspecified'),
        "n": int(len(y)), "positive_rate": float(y.mean()),
        "logloss": weighted_logloss(y, p), "brier": float(np.mean((p - y) ** 2)),
        "auroc": auroc, "ece_10bin": ece,
        "at_display_threshold": {"threshold": thr, "n_predicted_positive": int(pred.sum()),
                                 "precision": tp / pred.sum() if pred.sum() else None,
                                 "recall": tp / npos if npos else None,
                                 "recall_all_test_positives": tp / all_positive_rows if all_positive_rows else None},
    }
    lm = art["locus_model"]
    if lm.get("trained"):
        groups, multi = build_locus_groups(test_rows, margins, art["contract"]["features"])
        theta = np.concatenate(([math.log(lm["scale"]), lm["null_bias"]], lm["null_weights"]))
        nll = locus_nll(theta, groups, 0.0) if groups else None
        null_gr = [g for g in groups if g["target"] == 0]
        null_p = [locus_probabilities(theta, g["margins"], g["g"])[0] for g in null_gr]
        report["locus"] = {"n_groups": len(groups), "nll": nll, "n_null_groups": len(null_gr),
                           "mean_null_probability_on_null_groups": float(np.mean(null_p)) if null_p else None,
                           "n_multi_positive_groups_excluded": multi}
    else:
        report["locus"] = {"trained": False, "reason": lm.get("reason")}
    return report


def run(args) -> dict:
    contract = load_contract(Path(args.contract))
    table = Path(args.table)
    rows = read_table(table, contract)
    keys = [k.strip() for k in args.disjoint_keys.split(",") if k.strip()]
    for k in keys:
        if k not in ("group", "genome_id"):
            raise ProvenanceError(f"cannot enforce disjointness on unknown key {k}")
    check_disjoint(rows, keys)
    identity = split_identity(rows, keys)
    tr = [r for r in rows if r["split"] == "train"]
    ca = [r for r in rows if r["split"] == "calibration"]
    te = [r for r in rows if r["split"] == "test"]
    src = sha256_file(table)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    kinds = ["logistic", "gbdt"] if args.model == "both" else [args.model]
    summary = {}
    for kind in kinds:
        art = train(tr, ca, contract, kind, weighting=args.weighting, folds=args.folds, seed=args.seed,
                    threshold_rule=args.threshold_rule, min_null_groups=args.min_null_groups,
                    min_pos_groups=args.min_pos_groups, identity=identity,
                    intended_use=args.intended_use, synthetic=args.synthetic)
        path = out / f"model_{kind}.json"
        digest = write_frozen(art, path)                       # frozen + hashed BEFORE any test label is read
        report = evaluate_frozen(path, contract, te, digest, src)
        (out / f"test_report_{kind}.json").write_text(canonical_json(report))
        summary[kind] = {"model": str(path), "sha256": digest, "cv_logloss": art["training"]["cv"]["logloss"],
                         "test_logloss": report["logloss"]}
    # Primary choice uses the train-only CV score; the test report is never consulted.
    summary["primary_by_train_cv"] = min(kinds, key=lambda k: summary[k]["cv_logloss"])
    (out / "training_summary.json").write_text(canonical_json(summary))
    return summary


# ----------------------------------------------------------------------------- synthetic mechanism fixture

def synthetic_rows(contract, seed=7, n_groups=240):
    """Rows with PLANTED structure for mechanism tests; never a production dataset."""
    rng = np.random.default_rng(seed)
    names = contract["names"]
    idx = {n: i for i, n in enumerate(names)}
    rows = []
    n_split = {"train": 0.6, "calibration": 0.2, "test": 0.2}
    cuts = np.cumsum([n_split[s] for s in SPLITS])
    for gi in range(n_groups):
        u = gi / n_groups
        split = SPLITS[int(np.searchsorted(cuts, u, side="right"))] if u < cuts[-1] else "test"
        genome, fam = f"genome{gi % 40:03d}_{split}", f"fam{gi:04d}"
        n_cand = int(rng.integers(1, 5))
        has_correct = rng.random() > 0.3
        correct_at = int(rng.integers(0, n_cand)) if has_correct else -1
        for ci in range(n_cand):
            x = np.zeros(len(names))
            correct = ci == correct_at
            hom = float(np.clip(rng.normal(0.75 if correct else 0.35, 0.18), 0, 1))
            x[idx["fn.type_ec"]] = 1
            x[idx["role.catalytic"]] = 1
            x[idx["seq.group.homology.max_raw"]] = hom
            x[idx["seq.group.homology.reported"]] = 1
            x[idx["seq.src.diamond.max_raw"]] = hom
            x[idx["seq.src.diamond.reported"]] = 1
            x[idx["seq.n_sources_reporting"]] = 1
            x[idx["seq.n_groups_reporting"]] = 1
            x[idx["seq.max_raw_any"]] = hom
            if rng.random() < 0.7:     # structured neighbourhood supplied
                x[idx["ctx.nbr.available"]] = 1
                x[idx["ctx.nbr.support_max"]] = float(np.clip(rng.normal(0.7 if correct else 0.3, 0.2), 0, 1))
                x[idx["ctx.nbr.conflict_max"]] = float(np.clip(rng.normal(0.1 if correct else 0.4, 0.2), 0, 1))
            if rng.random() < 0.5:
                x[idx["ctx.competitor.max_raw"]] = float(np.clip(rng.normal(0.3 if correct else 0.6, 0.2), 0, 1))
                x[idx["ctx.competitor.n_with_claim"]] = float(rng.integers(0, 3))
            if rng.random() < 0.6:
                x[idx["asm.available"]] = 1
                x[idx["asm.log_contig_length"]] = float(rng.normal(12, 1))
                x[idx["asm.log_distance_to_contig_end"]] = float(rng.normal(9, 1.5))
                x[idx["asm.n_contigs"]] = float(math.log(1 + rng.integers(1, 40)))
            rows.append(dict(example_id=f"{fam}_{ci}", genome_id=genome, group=fam, split=split,
                             label=int(correct), x=x, locus_group_id=f"{fam}_loc", role="CATALYTIC"))
    return rows


def write_rows_jsonl(rows, names, path: Path):
    with open(path, "w") as fh:
        for r in rows:
            fh.write(json.dumps({
                "example_id": r["example_id"], "genome_id": r["genome_id"], "group": r["group"],
                "split": r["split"], "label": r["label"], "locus_group_id": r["locus_group_id"],
                "role": r["role"], "features": {n: float(v) for n, v in zip(names, r["x"])},
                **{k: r[k] for k in ('locus_members', 'locus_truth_source', 'label_basis', 'label_references',
                                     'function_type', 'function_id', 'candidate_key', 'sequence_sources', 'quality_evidence',
                                     'cohort_context_mode', 'sequence_contracts') if k in r},
            }, sort_keys=True) + "\n")


def make_mechanism_fixture(contract_path: Path, out_dir: Path) -> dict:
    """Write a synthetic-data fixture used ONLY to test the JVM <-> Python mechanism."""
    contract = load_contract(contract_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = synthetic_rows(contract)
    # The rows are regenerable from the seed; artifacts record only the train+calibration digest.
    check_disjoint(rows, ["group", "genome_id"])
    identity = split_identity(rows, ["group", "genome_id"])
    tr = [r for r in rows if r["split"] == "train"]
    ca = [r for r in rows if r["split"] == "calibration"]
    names = contract["names"]
    idx = {n: i for i, n in enumerate(names)}

    def vec(**kw):
        x = np.zeros(len(names))
        for k, v in kw.items():
            x[idx[k.replace("__", ".")]] = v
        return x

    base = dict(fn__type_ec=1, role__catalytic=1)
    ev = dict(seq__group__homology__max_raw=0.8, seq__group__homology__reported=1, seq__src__diamond__max_raw=0.8,
              seq__src__diamond__reported=1, seq__n_sources_reporting=1, seq__n_groups_reporting=1, seq__max_raw_any=0.8)
    vectors = {
        "empty": vec(**base),
        "strong_sequence": vec(**base, **ev),
        "strong_sequence_favourable_context": vec(**base, **ev, ctx__nbr__available=1, ctx__nbr__support_max=0.9,
                                                  ctx__nbr__conflict_max=0.0),
        "strong_sequence_unfavourable_context": vec(**base, **ev, ctx__nbr__available=1, ctx__nbr__support_max=0.05,
                                                    ctx__nbr__conflict_max=0.9, ctx__competitor__max_raw=0.9,
                                                    ctx__competitor__n_with_claim=2),
        "weak_sequence": vec(**base, seq__group__homology__max_raw=0.15, seq__group__homology__reported=1,
                             seq__src__diamond__max_raw=0.15, seq__src__diamond__reported=1,
                             seq__n_sources_reporting=1, seq__n_groups_reporting=1, seq__max_raw_any=0.15),
    }
    result = {"synthetic_mechanism_fixture": True, "contract_id": contract["id"],
              "note": "Planted-structure synthetic data. Verifies mechanisms only; carries no accuracy claim.",
              "vectors": {k: {names[i]: float(v) for i, v in enumerate(x) if v != 0.0} for k, x in vectors.items()},
              "expected": {}, "locus_cases": []}
    for kind in ("logistic", "gbdt"):
        art = train(tr, ca, contract, kind, weighting="group", folds=3, seed=1, identity=identity,
                    intended_use="mechanism_test", synthetic=True, min_null_groups=10, min_pos_groups=10)
        path = out_dir / f"mechanism_fixture_{kind}.json"
        write_frozen(art, path)
        loaded = load_artifact(path)
        result["expected"][kind] = {k: {"margin": predict_margin(loaded, x), "probability": predict_proba(loaded, x)}
                                    for k, x in vectors.items()}
        if loaded["locus_model"]["trained"]:
            cases = {
                "singleton_strong": ["strong_sequence"],
                "singleton_weak": ["weak_sequence"],
                "pair_strong_weak": ["strong_sequence", "weak_sequence"],
                "pair_two_weak": ["weak_sequence", "weak_sequence_b"],
            }
            vectors["weak_sequence_b"] = vectors["weak_sequence"].copy()
            result["vectors"]["weak_sequence_b"] = result["vectors"]["weak_sequence"]
            for cname, members in cases.items():
                p0, pc = assign_locus(loaded, [vectors[m] for m in members])
                result["locus_cases"].append({"model": kind, "case": cname, "members": members,
                                              "null_probability": p0, "candidate_probabilities": pc})
    (out_dir / "mechanism_fixture_vectors.json").write_text(canonical_json(result))
    (out_dir / "feature_contract_v1.json").write_text(Path(contract_path).read_text())
    return result


# ----------------------------------------------------------------------------- CLI

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("train", help="train, freeze, then evaluate once on the test split")
    t.add_argument("--table", required=True, help="training table (.jsonl or .tsv)")
    t.add_argument("--contract", required=True, help="feature contract JSON (FeatureContract.toJson)")
    t.add_argument("--out-dir", required=True)
    t.add_argument("--model", choices=["logistic", "gbdt", "both"], default="both")
    t.add_argument("--weighting", choices=["none", "genome", "group", "genome_group"], default="group")
    t.add_argument("--disjoint-keys", default="group,genome_id")
    t.add_argument("--folds", type=int, default=5)
    t.add_argument("--seed", type=int, default=0)
    t.add_argument("--threshold-rule", default="max_f1")
    t.add_argument("--min-null-groups", type=int, default=20)
    t.add_argument("--min-pos-groups", type=int, default=20)
    t.add_argument("--intended-use", choices=["production", "mechanism_test"], default="production")
    t.add_argument("--synthetic", action="store_true", help="mark the training data as synthetic")
    f = sub.add_parser("make-mechanism-fixture", help="write the synthetic JVM<->Python agreement fixture")
    f.add_argument("--contract", required=True)
    f.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    try:
        if args.cmd == "train":
            if args.synthetic and args.intended_use == "production":
                raise ProvenanceError("synthetic data cannot produce a production artifact")
            print(json.dumps(run(args), indent=1))
        else:
            make_mechanism_fixture(Path(args.contract), Path(args.out_dir))
    except ProvenanceError as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
