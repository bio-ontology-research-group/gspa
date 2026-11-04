#!/usr/bin/env python3
"""Positive-unlabelled (PU) trainer, PU reliability calibration, and PU marginal-likelihood null layer.

Reuses ``scripts/train_context_model.py`` without modifying it. Unconfirmed candidate rows are treated
strictly as unlabelled (never as negatives or false positives); explicit negatives are excluded from
every fit and from prior/label-frequency estimation and kept for evaluation only.
"""
from __future__ import annotations

import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import optimize

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts"))
import train_context_model as T  # noqa: E402

PU_TRAINER_VERSION = T.TRAINER_VERSION + "+pu1"
PU_PROBABILITY_TARGET = "experimental_support_under_positive_unlabelled_policy"
PU_LOCUS_PROBABILITY_TARGET = (
    "none of the proposed candidates is experimentally supported, beyond what the label frequency c explains"
)
PU_DISPLAY_THRESHOLD_NOTE = "F1 against labelled positives; unlabelled rows are not negatives"
PU_TREE_OBJECTIVE_DOC = (
    "newton_boosted_unbiased_pu_logistic_loss: init=logit(pi), g_pos=-pi*w/W_P, h_pos=0, "
    "g_unl=w*sigmoid(z)/W_U, h_unl=w*sigmoid(z)*(1-sigmoid(z))/W_U, leaf_denominator_clipped>=1e-12"
)
UNLABELLED_BASES = (None, "", "unspecified", "closed_world_unconfirmed", "no_truth_not_negative")
POSITIVE_BASES = (None, "", "unspecified", "positive_truth")


# ----------------------------------------------------------------------------- row semantics

def classify_pu_rows(rows: list[dict]) -> tuple[list[dict], np.ndarray, np.ndarray, int]:
    """Filter out explicit_negative_truth rows and return (kept_rows, is_pos, is_unl, n_explicit_neg).

    Unconfirmed rows are never called negatives anywhere in this module.
    """
    kept: list[dict] = []
    pos_flags: list[bool] = []
    unl_flags: list[bool] = []
    n_explicit_neg = 0
    for r in rows:
        basis = r.get("label_basis")
        label = int(r["label"])
        if basis == "explicit_negative_truth":
            if label != 0:
                raise T.ProvenanceError(
                    f"row {r.get('example_id')!r} has label={label} with explicit_negative_truth"
                )
            n_explicit_neg += 1
            continue
        if label == 1:
            if basis not in POSITIVE_BASES:
                raise T.ProvenanceError(
                    f"positive row {r.get('example_id')!r} has unexpected label_basis {basis!r}"
                )
            kept.append(r)
            pos_flags.append(True)
            unl_flags.append(False)
        elif label == 0:
            if basis not in UNLABELLED_BASES:
                raise T.ProvenanceError(
                    f"unlabelled row {r.get('example_id')!r} has unexpected label_basis {basis!r}"
                )
            kept.append(r)
            pos_flags.append(False)
            unl_flags.append(True)
        else:
            raise T.ProvenanceError(f"row {r.get('example_id')!r}: label must be 0 or 1, got {label!r}")
    return kept, np.array(pos_flags, dtype=bool), np.array(unl_flags, dtype=bool), n_explicit_neg


def _infer_masks(y: np.ndarray, is_pos: np.ndarray | None, is_unl: np.ndarray | None) -> tuple[np.ndarray, np.ndarray]:
    if is_pos is None:
        is_pos = np.asarray(y, dtype=float) == 1.0
    else:
        is_pos = np.asarray(is_pos, dtype=bool)
    if is_unl is None:
        is_unl = ~is_pos
    else:
        is_unl = np.asarray(is_unl, dtype=bool)
    if not np.any(is_pos) or not np.any(is_unl):
        raise T.ProvenanceError("PU training requires both positive and unlabelled rows")
    return is_pos, is_unl


def row_block_key(row: dict) -> tuple[str, str]:
    """Block key (genome_id, function_type) for pairwise PU-GO ranking."""
    gid = str(row["genome_id"])
    ft = row.get("function_type")
    if not ft:
        role = row.get("role")
        ft = "EC" if role == "CATALYTIC" else "GO"
    return (gid, str(ft))


def competition_key(row: dict) -> tuple | None:
    """Competition key (genome_id, function_type, function_id, role) with synthetic locus_group_id fallback."""
    gid = str(row["genome_id"])
    ft = row.get("function_type")
    fid = row.get("function_id")
    role = row.get("role")
    if ft not in (None, "") and fid not in (None, ""):
        norm_fid = str(fid).strip().replace("_", ":")
        if str(ft) == "EC" and not norm_fid.upper().startswith("EC:"):
            norm_fid = "EC:" + norm_fid
        return (gid, str(ft), norm_fid, str(role or ""))
    lg = row.get("locus_group_id")
    if lg not in (None, ""):
        return (gid, str(ft or "LOCUS"), str(lg), str(role or ""))
    return None


# ----------------------------------------------------------------------------- PU risk & objectives

def pu_negative_risk(margins: np.ndarray, w: np.ndarray, is_pos: np.ndarray, is_unl: np.ndarray, pi: float) -> float:
    """Estimated unlabelled-minus-positive risk component R_neg(theta) = R_U^-(theta) - pi * R_P^-(theta)."""
    z = np.asarray(margins, dtype=float)
    w = np.asarray(w, dtype=float)
    pos, unl = _infer_masks(z, is_pos, is_unl)
    wp = float(np.sum(w[pos]))
    wu = float(np.sum(w[unl]))
    ell_0 = np.logaddexp(0.0, z)
    r_u_neg = float(np.sum(w[unl] * ell_0[unl]) / wu)
    r_p_neg = float(np.sum(w[pos] * ell_0[pos]) / wp)
    return float(r_u_neg - float(pi) * r_p_neg)


def pu_unbiased_risk(margins: np.ndarray, w: np.ndarray, is_pos: np.ndarray, is_unl: np.ndarray, pi: float) -> float:
    """Unbiased PU risk (du Plessis 2015) without L2 penalty:
    pi * R_P^+(theta) + R_U^-(theta) - pi * R_P^-(theta).
    """
    z = np.asarray(margins, dtype=float)
    w = np.asarray(w, dtype=float)
    pos, unl = _infer_masks(z, is_pos, is_unl)
    wp = float(np.sum(w[pos]))
    wu = float(np.sum(w[unl]))
    ell_1_pos = np.logaddexp(0.0, -z[pos])
    ell_0_pos = np.logaddexp(0.0, z[pos])
    ell_0_unl = np.logaddexp(0.0, z[unl])
    r_p_pos = float(np.sum(w[pos] * ell_1_pos) / wp)
    r_p_neg = float(np.sum(w[pos] * ell_0_pos) / wp)
    r_u_neg = float(np.sum(w[unl] * ell_0_unl) / wu)
    return float(float(pi) * r_p_pos + r_u_neg - float(pi) * r_p_neg)


def pn_risk(margins: np.ndarray, y_true: np.ndarray, w: np.ndarray | None = None) -> float:
    """Weighted logistic risk on binary labels (for unbiasedness verification against true labels)."""
    z = np.asarray(margins, dtype=float)
    y = np.asarray(y_true, dtype=float)
    w_arr = np.ones(len(z), dtype=float) if w is None else np.asarray(w, dtype=float)
    loss = np.logaddexp(0.0, z) - y * z
    return float(np.sum(w_arr * loss) / np.sum(w_arr))


def pu_unbiased_objective_and_grad(
    theta: np.ndarray,
    Z: np.ndarray,
    w: np.ndarray,
    is_pos: np.ndarray,
    is_unl: np.ndarray,
    pi: float,
    l2: float,
) -> tuple[float, np.ndarray]:
    """Convex unbiased PU logistic loss and analytic gradient on standardised features Z."""
    b0, beta = float(theta[0]), theta[1:]
    z = b0 + Z @ beta
    wp = float(np.sum(w[is_pos]))
    wu = float(np.sum(w[is_unl]))
    # Using ell(z, 1) - ell(z, 0) = -z on positives:
    pos_term = -float(pi) * float(np.sum(w[is_pos] * z[is_pos]) / wp)
    unl_term = float(np.sum(w[is_unl] * np.logaddexp(0.0, z[is_unl])) / wu)
    reg = 0.5 * float(l2) * float(beta @ beta)
    loss = pos_term + unl_term + reg

    r = np.zeros(len(z), dtype=float)
    r[is_pos] = -float(pi) * (w[is_pos] / wp)
    r[is_unl] = (w[is_unl] / wu) * T._sigmoid(z[is_unl])
    grad = np.concatenate(([float(r.sum())], Z.T @ r + float(l2) * beta))
    return float(loss), grad


def pu_nnpu_objective_and_grad(
    theta: np.ndarray,
    Z: np.ndarray,
    w: np.ndarray,
    is_pos: np.ndarray,
    is_unl: np.ndarray,
    pi: float,
    l2: float,
    k: float = 100.0,
    penalty: float = 1e6,
) -> tuple[float, np.ndarray]:
    """Smooth non-negative PU objective:
    pi * R_P^+(theta) + softplus_k(R_neg(theta)) + (penalty/2) * max(0, -R_neg(theta))^2 + (l2/2)||beta||^2.
    """
    b0, beta = float(theta[0]), theta[1:]
    z = b0 + Z @ beta
    wp = float(np.sum(w[is_pos]))
    wu = float(np.sum(w[is_unl]))

    ell_1_pos = np.logaddexp(0.0, -z[is_pos])
    ell_0_pos = np.logaddexp(0.0, z[is_pos])
    ell_0_unl = np.logaddexp(0.0, z[is_unl])

    r_p_pos = float(np.sum(w[is_pos] * ell_1_pos) / wp)
    r_p_neg = float(np.sum(w[is_pos] * ell_0_pos) / wp)
    r_u_neg = float(np.sum(w[is_unl] * ell_0_unl) / wu)
    r_neg = r_u_neg - float(pi) * r_p_neg

    sp = float(np.logaddexp(0.0, k * r_neg) / k)
    sp_deriv = float(T._sigmoid(np.array(k * r_neg)))
    viol = max(0.0, -r_neg)
    quad_pen = 0.5 * penalty * (viol ** 2)
    # Derivative of softplus_k(r_neg) + 0.5 * penalty * max(0, -r_neg)^2 w.r.t. r_neg:
    coeff_r_neg = sp_deriv - (penalty * viol if r_neg < 0.0 else 0.0)

    sig_z = T._sigmoid(z)
    r = np.zeros(len(z), dtype=float)
    # Gradient of pi * R_P^+ + coeff_r_neg * (R_U^- - pi * R_P^-) w.r.t. z:
    r[is_pos] = float(pi) * (w[is_pos] / wp) * ((sig_z[is_pos] - 1.0) - coeff_r_neg * sig_z[is_pos])
    r[is_unl] = coeff_r_neg * (w[is_unl] / wu) * sig_z[is_unl]

    reg = 0.5 * float(l2) * float(beta @ beta)
    loss = float(pi) * r_p_pos + sp + quad_pen + reg
    grad = np.concatenate(([float(r.sum())], Z.T @ r + float(l2) * beta))
    return float(loss), grad


def _prepare_rank_blocks(is_pos: np.ndarray, is_unl: np.ndarray, blocks: list) -> list[tuple[np.ndarray, np.ndarray]]:
    by_block_pos: dict = {}
    by_block_unl: dict = {}
    for idx, b in enumerate(blocks):
        if is_pos[idx]:
            by_block_pos.setdefault(b, []).append(idx)
        elif is_unl[idx]:
            by_block_unl.setdefault(b, []).append(idx)
    pairs = []
    for b in sorted(set(by_block_pos) & set(by_block_unl)):
        p_idx = np.array(by_block_pos[b], dtype=int)
        u_idx = np.array(by_block_unl[b], dtype=int)
        if len(p_idx) > 0 and len(u_idx) > 0:
            pairs.append((p_idx, u_idx))
    return pairs


def pu_rank_risk(
    margins: np.ndarray,
    w: np.ndarray,
    is_pos: np.ndarray,
    is_unl: np.ndarray,
    pi: float,
    blocks: list,
) -> float:
    """PU-GO ranking risk without L2: pi * R_P^+ + pairwise ranking loss inside (genome, function_type) blocks."""
    z = np.asarray(margins, dtype=float)
    w = np.asarray(w, dtype=float)
    pos, unl = _infer_masks(z, is_pos, is_unl)
    wp = float(np.sum(w[pos]))
    r_p_pos = float(np.sum(w[pos] * np.logaddexp(0.0, -z[pos])) / wp)

    block_pairs = _prepare_rank_blocks(pos, unl, blocks)
    total_pair_w = 0.0
    rank_num = 0.0
    for p_idx, u_idx in block_pairs:
        wp_b = w[p_idx][:, None]
        wu_b = w[u_idx][None, :]
        w_mat = wp_b * wu_b
        diff = z[p_idx][:, None] - z[u_idx][None, :]
        total_pair_w += float(np.sum(w_mat))
        rank_num += float(np.sum(w_mat * np.logaddexp(0.0, -diff)))
    r_rank = (rank_num / total_pair_w) if total_pair_w > 0.0 else 0.0
    return float(float(pi) * r_p_pos + r_rank)


def pu_rank_objective_and_grad(
    theta: np.ndarray,
    Z: np.ndarray,
    w: np.ndarray,
    is_pos: np.ndarray,
    is_unl: np.ndarray,
    pi: float,
    blocks: list,
    l2: float,
    _block_pairs: list[tuple[np.ndarray, np.ndarray]] | None = None,
) -> tuple[float, np.ndarray]:
    """PU-GO ranking loss and analytic gradient on standardised features Z."""
    b0, beta = float(theta[0]), theta[1:]
    z = b0 + Z @ beta
    wp = float(np.sum(w[is_pos]))
    r_p_pos = float(np.sum(w[is_pos] * np.logaddexp(0.0, -z[is_pos])) / wp)

    r = np.zeros(len(z), dtype=float)
    r[is_pos] = float(pi) * (w[is_pos] / wp) * (T._sigmoid(z[is_pos]) - 1.0)

    if _block_pairs is None:
        _block_pairs = _prepare_rank_blocks(is_pos, is_unl, blocks)

    total_pair_w = sum(float(np.sum(w[p_idx]) * np.sum(w[u_idx])) for p_idx, u_idx in _block_pairs)
    rank_num = 0.0
    if total_pair_w > 0.0:
        for p_idx, u_idx in _block_pairs:
            wp_b = w[p_idx][:, None]
            wu_b = w[u_idx][None, :]
            w_mat = wp_b * wu_b
            diff = z[p_idx][:, None] - z[u_idx][None, :]
            rank_num += float(np.sum(w_mat * np.logaddexp(0.0, -diff)))
            # d/d(diff) log(1 + exp(-diff)) = -sigmoid(-diff)
            s_mat = (w_mat / total_pair_w) * T._sigmoid(-diff)
            r[p_idx] -= np.sum(s_mat, axis=1)
            r[u_idx] += np.sum(s_mat, axis=0)
    r_rank = (rank_num / total_pair_w) if total_pair_w > 0.0 else 0.0
    reg = 0.5 * float(l2) * float(beta @ beta)
    loss = float(pi) * r_p_pos + r_rank + reg
    grad = np.concatenate(([float(r.sum())], Z.T @ r + float(l2) * beta))
    return float(loss), grad


# ----------------------------------------------------------------------------- model classes

class PuLogisticModel(T.LogisticModel):
    """Unbiased PU logistic regression with Kiryo et al. (2017) non-negative risk correction when R_neg < 0."""

    def __init__(self, l2: float):
        super().__init__(l2)
        self.clamp_applied = False
        self.r_neg_before_clamp = 0.0
        self.r_neg = 0.0

    def fit(self, X, y, w, *, is_pos=None, is_unl=None, pi: float = 0.5, blocks=None, force_clamp: bool = False):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        w = np.asarray(w, dtype=float)
        pos, unl = _infer_masks(y, is_pos, is_unl)

        self.mean = X.mean(axis=0)
        sd = X.std(axis=0)
        self.scale = np.where(sd < 1e-12, 1.0, sd)
        Z = (X - self.mean) / self.scale
        d = Z.shape[1]

        res = optimize.minimize(
            lambda th: pu_unbiased_objective_and_grad(th, Z, w, pos, unl, pi, self.l2),
            np.zeros(d + 1),
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-8},
        )
        theta = res.x.copy()
        z = theta[0] + Z @ theta[1:]
        r_neg_initial = pu_negative_risk(z, w, pos, unl, pi)
        self.r_neg_before_clamp = float(r_neg_initial)

        if r_neg_initial < 0.0 or force_clamp:
            self.clamp_applied = True
            theta = self._fit_clamped(Z, w, pos, unl, pi, d)
            z = theta[0] + Z @ theta[1:]
            self.r_neg = float(pu_negative_risk(z, w, pos, unl, pi))
        else:
            self.clamp_applied = False
            self.r_neg = float(r_neg_initial)

        self.intercept = float(theta[0])
        self.weights = theta[1:].copy()
        return self

    def _fit_clamped(self, Z: np.ndarray, w: np.ndarray, pos: np.ndarray, unl: np.ndarray, pi: float, d: int) -> np.ndarray:
        """Fit under the non-negative PU risk (softplus_k + quadratic barrier + Kiryo gradient ascent fallback)."""
        res = optimize.minimize(
            lambda th: pu_nnpu_objective_and_grad(th, Z, w, pos, unl, pi, self.l2, k=100.0, penalty=1e6),
            np.zeros(d + 1),
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-8},
        )
        theta = res.x.copy()
        wp = float(np.sum(w[pos]))
        wu = float(np.sum(w[unl]))

        # Kiryo et al. (2017) outer alternation: if R_neg(theta) < 0, ascend R_neg(theta).
        step = 0.05
        for _ in range(200):
            z = theta[0] + Z @ theta[1:]
            r_neg = pu_negative_risk(z, w, pos, unl, pi)
            if r_neg >= 0.0:
                break
            sig_z = T._sigmoid(z)
            # Gradient of R_neg w.r.t. z:
            r_vec = np.zeros(len(z), dtype=float)
            r_vec[unl] = (w[unl] / wu) * sig_z[unl]
            r_vec[pos] = -float(pi) * (w[pos] / wp) * sig_z[pos]
            grad_r_neg = np.concatenate(([float(r_vec.sum())], Z.T @ r_vec))
            norm = float(np.linalg.norm(grad_r_neg))
            if norm < 1e-12:
                break
            theta = theta + step * (grad_r_neg / norm)

        # At theta = 0, R_neg(0) = (1 - pi) * log(2) > 0; bisect if still below 0.
        z = theta[0] + Z @ theta[1:]
        if pu_negative_risk(z, w, pos, unl, pi) < 0.0:
            lo, hi = 0.0, 1.0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                cand = (1.0 - mid) * theta
                if pu_negative_risk(cand[0] + Z @ cand[1:], w, pos, unl, pi) >= 0.0:
                    hi = mid
                else:
                    lo = mid
            theta = (1.0 - hi) * theta
        return theta


class PuRankLogisticModel(T.LogisticModel):
    """PU-GO pairwise ranking logistic regression inside (genome_id, function_type) blocks."""

    def __init__(self, l2: float):
        super().__init__(l2)
        self.clamp_applied = False
        self.r_neg_before_clamp = 0.0
        self.r_neg = 0.0

    def fit(self, X, y, w, *, is_pos=None, is_unl=None, pi: float = 0.5, blocks=None):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        w = np.asarray(w, dtype=float)
        pos, unl = _infer_masks(y, is_pos, is_unl)
        if blocks is None:
            blocks = [("default", "default")] * len(y)

        self.mean = X.mean(axis=0)
        sd = X.std(axis=0)
        self.scale = np.where(sd < 1e-12, 1.0, sd)
        Z = (X - self.mean) / self.scale
        d = Z.shape[1]
        block_pairs = _prepare_rank_blocks(pos, unl, blocks)

        res = optimize.minimize(
            lambda th: pu_rank_objective_and_grad(th, Z, w, pos, unl, pi, blocks, self.l2, _block_pairs=block_pairs),
            np.zeros(d + 1),
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-8},
        )
        self.intercept = float(res.x[0])
        self.weights = res.x[1:].copy()
        z = self.intercept + Z @ self.weights
        self.r_neg_before_clamp = float(pu_negative_risk(z, w, pos, unl, pi))
        self.r_neg = self.r_neg_before_clamp
        self.clamp_applied = False
        return self


def pu_gbdt_grad_hess(
    z: np.ndarray,
    w: np.ndarray,
    is_pos: np.ndarray,
    is_unl: np.ndarray,
    pi: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-sample gradient and Hessian for PU GBDT boosting rounds.

    Positives: g_i = -pi * w_i / W_P, h_i = 0.
    Unlabelled: g_i = w_i * sigmoid(z_i) / W_U, h_i = w_i * sigmoid(z_i) * (1 - sigmoid(z_i)) / W_U.
    """
    z = np.asarray(z, dtype=float)
    w = np.asarray(w, dtype=float)
    pos, unl = _infer_masks(z, is_pos, is_unl)
    wp = float(np.sum(w[pos]))
    wu = float(np.sum(w[unl]))
    p = T._sigmoid(z)

    g = np.zeros(len(z), dtype=float)
    h = np.zeros(len(z), dtype=float)
    g[pos] = -float(pi) * w[pos] / wp
    g[unl] = w[unl] * p[unl] / wu
    h[pos] = 0.0
    h[unl] = w[unl] * p[unl] * (1.0 - p[unl]) / wu
    return g, h


class PuGbdtModel(T.GbdtModel):
    """Newton-boosted shallow regression trees on the unbiased PU logistic loss."""

    def __init__(self, n_trees: int, depth: int, lr: float, min_leaf: int = 5, lam: float = 1.0):
        super().__init__(n_trees, depth, lr, min_leaf=min_leaf, lam=max(float(lam), 1e-12))
        self.clamp_applied = False
        self.r_neg_before_clamp = 0.0
        self.r_neg = 0.0

    def fit(self, X, y, w, *, is_pos=None, is_unl=None, pi: float = 0.5, blocks=None, force_clamp: bool = False):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        w = np.asarray(w, dtype=float)
        pos, unl = _infer_masks(y, is_pos, is_unl)

        pi_c = float(np.clip(pi, 1e-6, 1.0 - 1e-6))
        self.init = math.log(pi_c / (1.0 - pi_c))
        z = np.full(len(y), self.init, dtype=float)
        self.trees = []
        for _ in range(self.n_trees):
            g, h = pu_gbdt_grad_hess(z, w, pos, unl, pi_c)
            tree = self._grow(X, g, h)
            self.trees.append(tree)
            z = z + self.lr * self._eval_tree(tree, X)

        r_neg_initial = pu_negative_risk(z, w, pos, unl, pi_c)
        self.r_neg_before_clamp = float(r_neg_initial)
        if r_neg_initial < 0.0 or force_clamp:
            self.clamp_applied = True
            delta = z - self.init
            lo, hi = 0.0, 1.0
            # Find largest scale factor alpha in [0, 1] on tree predictions with R_neg >= 0; if init itself < 0, shift init.
            if pu_negative_risk(np.full(len(y), min(self.init, 0.0)), w, pos, unl, pi_c) >= 0.0:
                if pu_negative_risk(np.full(len(y), self.init), w, pos, unl, pi_c) < 0.0:
                    self.init = 0.0
                for _ in range(50):
                    mid = 0.5 * (lo + hi)
                    if pu_negative_risk(self.init + mid * delta, w, pos, unl, pi_c) >= 0.0:
                        lo = mid
                    else:
                        hi = mid
                for t in self.trees:
                    t["value"] = [float(lo * v) for v in t["value"]]
            z = self.margin(X)
            self.r_neg = float(pu_negative_risk(z, w, pos, unl, pi_c))
        else:
            self.clamp_applied = False
            self.r_neg = float(r_neg_initial)
        return self


# ----------------------------------------------------------------------------- prior, CV & calibration

def estimate_prior(oof_margin: np.ndarray, y_labelled: np.ndarray, w: np.ndarray) -> dict:
    """Elkan-Noto e1 estimator on train out-of-fold margins of a PN logistic model.

    c = sum_P w * sigmoid(z) / W_P
    labelled_rate r = W_P / (W_P + W_U)
    pi_hat = min(max(r / c, r), 0.95)
    """
    z = np.asarray(oof_margin, dtype=float)
    y = np.asarray(y_labelled, dtype=float)
    w = np.asarray(w, dtype=float)
    pos, unl = _infer_masks(y, y == 1.0, y == 0.0)
    wp = float(np.sum(w[pos]))
    wu = float(np.sum(w[unl]))
    p = T._sigmoid(z)
    c = float(np.sum(w[pos] * p[pos]) / wp)
    r = float(wp / (wp + wu))
    pi_hat = float(min(max(r / max(c, 1e-12), r), 0.95))
    return {
        "pi": pi_hat,
        "c": c,
        "labelled_rate": r,
        "method": "elkan_noto_e1_train_oof",
    }


def pu_default_candidates(kind: str, objective: str):
    if kind == "logistic":
        if objective == "pu":
            return [({"l2": l2}, (lambda l2=l2: PuLogisticModel(l2))) for l2 in (0.1, 1.0, 10.0, 100.0)]
        if objective == "pu_rank":
            return [({"l2": l2}, (lambda l2=l2: PuRankLogisticModel(l2))) for l2 in (0.1, 1.0, 10.0, 100.0)]
        raise T.ProvenanceError(f"unknown logistic objective {objective!r}")
    if kind == "gbdt":
        if objective != "pu":
            raise T.ProvenanceError(f"gbdt only supports objective 'pu', got {objective!r}")
        return [
            ({"n_trees": 40, "depth": d, "learning_rate": 0.15}, (lambda d=d: PuGbdtModel(40, d, 0.15)))
            for d in (2, 3)
        ]
    raise T.ProvenanceError(f"unknown model kind {kind!r}")


def pu_cv_margins(
    make_model,
    X: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    groups: list[str],
    k: int,
    seed: int,
    *,
    is_pos: np.ndarray,
    is_unl: np.ndarray,
    pi: float,
    blocks: list | None = None,
) -> np.ndarray:
    """Group-wise out-of-fold margins for a PU model."""
    folds = np.array(T.grouped_folds(groups, k, seed))
    oof = np.zeros(len(y), dtype=float)
    for f in range(k):
        tr, te = folds != f, folds == f
        if te.sum() == 0 or tr.sum() == 0 or not np.any(is_pos[tr]) or not np.any(is_unl[tr]):
            raise T.ProvenanceError("a CV training fold lacks positive or unlabelled examples; use fewer folds or more groups")
        tr_blocks = [blocks[i] for i in np.nonzero(tr)[0]] if blocks is not None else None
        m = make_model().fit(X[tr], y[tr], w[tr], is_pos=is_pos[tr], is_unl=is_unl[tr], pi=pi, blocks=tr_blocks)
        oof[te] = m.margin(X[te])
    return oof


def select_pu_hyperparameters(
    candidates,
    X: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    groups: list[str],
    k: int,
    seed: int,
    *,
    objective: str,
    is_pos: np.ndarray,
    is_unl: np.ndarray,
    pi: float,
    blocks: list | None = None,
):
    """Select hyper-parameters by grouped CV on train minimising the PU risk R at the chosen pi."""
    scored = []
    for params, make in candidates:
        oof = pu_cv_margins(make, X, y, w, groups, k, seed, is_pos=is_pos, is_unl=is_unl, pi=pi, blocks=blocks)
        if objective == "pu":
            risk = pu_unbiased_risk(oof, w, is_pos, is_unl, pi)
        elif objective == "pu_rank":
            risk = pu_rank_risk(oof, w, is_pos, is_unl, pi, blocks or [("default", "default")] * len(y))
        else:
            raise T.ProvenanceError(f"unknown objective {objective!r}")
        scored.append((float(risk), params, make, oof))
    scored.sort(key=lambda t: t[0])
    return scored[0], [(s, p) for s, p, _, _ in scored]


def pu_reliability_table(
    p: np.ndarray,
    y_labelled: np.ndarray,
    c: float,
    n_bins: int = 10,
) -> tuple[list[dict], float, float]:
    """Ten equal-width bins of p on [0, 1]; returns (bins_table, pu_ece, pn_ece_reference)."""
    p = np.asarray(p, dtype=float)
    y = np.asarray(y_labelled, dtype=float)
    n = len(p)
    if n == 0:
        return [], 0.0, 0.0
    bins = np.clip((p * n_bins).astype(int), 0, n_bins - 1)
    table = []
    pu_ece = 0.0
    pn_ece = 0.0
    for b in range(n_bins):
        mask = bins == b
        nb = int(np.sum(mask))
        if nb > 0:
            mean_p = float(np.mean(p[mask]))
            obs_rate = float(np.mean(y[mask]))
            exp_rate = float(c * mean_p)
            weight = nb / n
            pu_ece += weight * abs(obs_rate - exp_rate)
            pn_ece += weight * abs(obs_rate - mean_p)
        else:
            mean_p = None
            obs_rate = None
            exp_rate = None
        table.append(
            {
                "bin": b,
                "lower": round(b / n_bins, 4),
                "upper": round((b + 1) / n_bins, 4),
                "n": nb,
                "mean_p": mean_p,
                "observed_labelled_rate": obs_rate,
                "expected_labelled_rate": exp_rate,
            }
        )
    return table, float(pu_ece), float(pn_ece)


def fit_pu_temperature(cal_margin: np.ndarray, y_cal: np.ndarray, c: float) -> dict:
    """Choose temperature a on a 61-point log-spaced grid in [0.25, 4.0] minimising PU-ECE (ties -> closest to 1)."""
    cal_margin = np.asarray(cal_margin, dtype=float)
    y_cal = np.asarray(y_cal, dtype=float)
    grid = np.geomspace(0.25, 4.0, 61)
    grid[30] = 1.0  # exact 2^0

    _, ece_before, _ = pu_reliability_table(T._sigmoid(cal_margin), y_cal, c)
    best_a = 1.0
    best_ece = float("inf")
    for a in grid:
        a_val = float(a)
        p = T._sigmoid(a_val * cal_margin)
        _, ece, _ = pu_reliability_table(p, y_cal, c)
        if ece < best_ece - 1e-15 or (abs(ece - best_ece) <= 1e-15 and abs(a_val - 1.0) < abs(best_a - 1.0)):
            best_ece = ece
            best_a = a_val

    return {
        "method": "platt",
        "a": float(best_a),
        "b": 0.0,
        "fitted_on": "calibration",
        "n": int(len(y_cal)),
        "pu_temperature": True,
        "pu_ece_before": float(ece_before),
        "pu_ece_after": float(best_ece),
        "c": float(c),
        "selected_on": "pu_ece",
    }


# ----------------------------------------------------------------------------- null layer (PU marginal likelihood)

def build_pu_groups(rows: list[dict], margins, names: list[str], *, min_members: int = 2) -> tuple[list[dict], int]:
    """Build PU locus competitions keyed by (genome_id, function_type, function_id, role).

    Excludes groups with >= 2 labelled positives (counted in ``multi``).
    Single-positive group -> target = index + 1.
    No-label group -> target = 0 ("no label observed", NOT "null is true").
    """
    idx = {n: i for i, n in enumerate(names)}
    by: dict[tuple, list[int]] = {}
    for i, r in enumerate(rows):
        if r.get("label_basis") == "explicit_negative_truth":
            continue
        key = competition_key(r)
        if key is not None:
            by.setdefault(key, []).append(i)
    out: list[dict] = []
    multi = 0
    for gid, members in by.items():
        if len(members) < min_members:
            continue
        pos = [j for j, i in enumerate(members) if rows[i]["label"] == 1]
        if len(pos) > 1:
            multi += 1
            continue
        m = [float(margins[i]) for i in members]
        f = [rows[i]["x"] for i in members]
        out.append(
            dict(
                gid=gid,
                margins=np.array(m, dtype=float),
                g=T._group_features(m, f, idx),
                target=(pos[0] + 1 if pos else 0),
                member_indices=list(members),
            )
        )
    return out, multi


def pu_locus_nll(theta: np.ndarray, groups: list[dict], c: float, l2: float) -> float:
    """PU marginal negative log-likelihood for the null layer.

    theta = [log s, b0, w].
    Labelled group (target > 0): -log softmax[target].
    No-label group (target == 0): -log(p0 + (1 - c) * sum_{j>=1} p_j), computed stably in log space.
    """
    s, b0, w = math.exp(float(theta[0])), float(theta[1]), np.asarray(theta[2:], dtype=float)
    c_val = float(np.clip(c, 0.0, 1.0))
    nll = 0.0
    for gr in groups:
        logits = np.concatenate(([b0 + float(w @ gr["g"])], s * gr["margins"]))
        mx = float(logits.max())
        lse = mx + math.log(float(np.exp(logits - mx).sum()))
        target = int(gr["target"])
        if target > 0:
            nll += lse - float(logits[target])
        else:
            if c_val >= 1.0 - 1e-15:
                nll += lse - float(logits[0])
            else:
                eff_logits = logits.copy()
                eff_logits[1:] += math.log(1.0 - c_val)
                mx_eff = float(eff_logits.max())
                lse_eff = mx_eff + math.log(float(np.exp(eff_logits - mx_eff).sum()))
                nll += lse - lse_eff
    return nll / max(len(groups), 1) + 0.5 * float(l2) * float(w @ w)


def uniform_pu_locus_nll(groups: list[dict], c: float) -> float:
    """Closed-form marginal NLL of the uniform model (p0 = p_j = 1 / (n + 1)) under the same c."""
    if not groups:
        return 0.0
    c_val = float(np.clip(c, 0.0, 1.0))
    total = 0.0
    for gr in groups:
        n = len(gr["margins"])
        if int(gr["target"]) > 0:
            total += math.log(n + 1.0)
        else:
            total += math.log(n + 1.0) - math.log(1.0 + (1.0 - c_val) * n)
    return float(total / len(groups))


def plain_uniform_nll(groups: list[dict]) -> float:
    """Plain mean log(n + 1) baseline across groups."""
    if not groups:
        return 0.0
    return float(np.mean([math.log(len(gr["margins"]) + 1.0) for gr in groups]))


def fit_pu_locus(groups: list[dict], c: float, l2: float) -> np.ndarray:
    """Fit null-layer parameters [log s, b0, w] by L-BFGS-B with the same bounds as T.fit_locus."""
    theta0 = np.zeros(2 + len(T.LOCUS_GROUP_FEATURES), dtype=float)
    bounds = [(-5.0, 5.0)] + [(-30.0, 30.0)] * (1 + len(T.LOCUS_GROUP_FEATURES))
    res = optimize.minimize(
        lambda t: pu_locus_nll(t, groups, c, l2),
        theta0,
        method="L-BFGS-B",
        bounds=bounds,
    )
    return res.x


# ----------------------------------------------------------------------------- train_pu

def _safe_label_basis_counts(rows: list[dict]) -> dict[str, int]:
    """Count label_basis values without ever placing unlabelled rows under a key containing 'negative'."""
    counts: Counter = Counter()
    for r in rows:
        b = r.get("label_basis", "unspecified")
        if b in (None, ""):
            b = "unspecified"
        elif b == "no_truth_not_negative":
            b = "unlabelled_export_no_truth"
        counts[str(b)] += 1
    return dict(sorted(counts.items()))


def train_pu(
    train_rows: list[dict],
    calib_rows: list[dict],
    contract: dict,
    kind: str,
    objective: str = "pu",
    *,
    pi: float | None = None,
    weighting: str = "group",
    folds: int = 5,
    seed: int = 42,
    threshold_rule: str = "max_f1",
    min_no_label_groups: int = 20,
    min_pos_groups: int = 20,
    min_cal_no_label_groups: int = 5,
    min_cal_pos_groups: int = 5,
    identity: dict,
    intended_use: str = "production",
    synthetic: bool = False,
) -> dict:
    """Fit a PU candidate-layer model and PU marginal-likelihood null layer; return a JVM-compatible artifact."""
    if any(r.get("split") != "train" for r in train_rows) or any(r.get("split") != "calibration" for r in calib_rows):
        raise T.ProvenanceError("train_pu requires split in {'train', 'calibration'}; test rows must never reach training")
    if any(r.get("quality_evidence") for r in train_rows + calib_rows):
        raise T.ProvenanceError("gene-model quality observations require a supported feature contract before training")
    if objective not in ("pu", "pu_rank"):
        raise T.ProvenanceError(f"unsupported PU objective {objective!r}; expected 'pu' or 'pu_rank'")

    names = contract["names"]
    data_sha = T.training_data_digest(train_rows, calib_rows)

    pu_train_rows, is_pos_tr, is_unl_tr, n_exp_neg_tr = classify_pu_rows(train_rows)
    pu_calib_rows, is_pos_ca, is_unl_ca, n_exp_neg_ca = classify_pu_rows(calib_rows)
    n_explicit_negatives_excluded = n_exp_neg_tr + n_exp_neg_ca

    if not np.any(is_pos_tr) or not np.any(is_unl_tr):
        raise T.ProvenanceError("train split needs both positive and unlabelled rows after excluding explicit negatives")
    if not np.any(is_pos_ca) or not np.any(is_unl_ca):
        raise T.ProvenanceError("calibration split needs both positive and unlabelled rows after excluding explicit negatives")

    Xtr, ytr = T.stack(pu_train_rows)
    Xca, yca = T.stack(pu_calib_rows)
    wtr = T.sample_weights(pu_train_rows, weighting)
    wca = T.sample_weights(pu_calib_rows, weighting)
    gtr = [r["group"] for r in pu_train_rows]
    blocks_tr = [row_block_key(r) for r in pu_train_rows]

    # Step 1: Estimate c and pi on train OOF margins of a PN logistic fit.
    (_, _, pn_make), _ = T.select_hyperparameters(T.default_candidates("logistic"), Xtr, ytr, wtr, gtr, folds, seed)
    pn_oof = T.cv_margins(pn_make, Xtr, ytr, wtr, gtr, folds, seed)
    prior_info = estimate_prior(pn_oof, ytr, wtr)
    c_val = float(prior_info["c"])
    labelled_rate = float(prior_info["labelled_rate"])
    if pi is None:
        pi_val = float(prior_info["pi"])
        pi_source = str(prior_info["method"])
    else:
        pi_val = float(pi)
        if not (0.0 < pi_val < 1.0):
            raise T.ProvenanceError(f"override pi must be in (0, 1), got {pi_val}")
        pi_source = f"override:{pi_val}"

    # Step 2: Select PU hyper-parameters by grouped CV on train minimising the PU risk at pi_val.
    candidates = pu_default_candidates(kind, objective)
    (cv_risk, params, make, oof), all_scores = select_pu_hyperparameters(
        candidates,
        Xtr,
        ytr,
        wtr,
        gtr,
        folds,
        seed,
        objective=objective,
        is_pos=is_pos_tr,
        is_unl=is_unl_tr,
        pi=pi_val,
        blocks=blocks_tr,
    )

    model = make().fit(Xtr, ytr, wtr, is_pos=is_pos_tr, is_unl=is_unl_tr, pi=pi_val, blocks=blocks_tr)
    r_neg_train_oof = float(pu_negative_risk(oof, wtr, is_pos_tr, is_unl_tr, pi_val))
    r_neg_train_in = float(pu_negative_risk(model.margin(Xtr), wtr, is_pos_tr, is_unl_tr, pi_val))

    if objective == "pu" and (r_neg_train_in < 0.0 or r_neg_train_oof < 0.0) and not model.clamp_applied:
        model = make().fit(
            Xtr, ytr, wtr, is_pos=is_pos_tr, is_unl=is_unl_tr, pi=pi_val, blocks=blocks_tr, force_clamp=True
        )
        r_neg_train_in = float(pu_negative_risk(model.margin(Xtr), wtr, is_pos_tr, is_unl_tr, pi_val))

    clamp_applied = bool(getattr(model, "clamp_applied", False))

    # Step 3: PU reliability calibration on the calibration split (temperature only, b = 0.0).
    cal_margin = model.margin(Xca)
    calibration = fit_pu_temperature(cal_margin, yca, c_val)
    p_cal = T._sigmoid(calibration["a"] * cal_margin + calibration["b"])
    threshold = T.choose_threshold(p_cal, yca, threshold_rule)
    threshold["note"] = PU_DISPLAY_THRESHOLD_NOTE

    # Step 4: Null layer (PU marginal likelihood) on competitions with >= 2 members.
    train_groups, train_multi = build_pu_groups(pu_train_rows, oof, names, min_members=2)
    cal_groups, cal_multi = build_pu_groups(pu_calib_rows, cal_margin, names, min_members=2)
    n_tr_no_label = sum(1 for g in train_groups if g["target"] == 0)
    n_tr_pos = len(train_groups) - n_tr_no_label
    n_ca_no_label = sum(1 for g in cal_groups if g["target"] == 0)
    n_ca_pos = len(cal_groups) - n_ca_no_label

    if not train_groups and not cal_groups:
        locus = {
            "trained": False,
            "reason": "no competitions with >= 2 members found in train or calibration",
            "n_train_groups": 0,
            "n_train_no_label_groups": 0,
            "n_train_single_positive_groups": 0,
            "n_multi_positive_groups_excluded": int(train_multi),
            "n_calibration_groups": 0,
            "n_calibration_no_label_groups": 0,
            "n_calibration_single_positive_groups": 0,
        }
    elif n_tr_no_label < min_no_label_groups or n_tr_pos < min_pos_groups:
        locus = {
            "trained": False,
            "reason": (
                f"train split has {n_tr_no_label} no-label groups and {n_tr_pos} single-positive groups; "
                f"need >= {min_no_label_groups} and >= {min_pos_groups} to learn a PU null outcome"
            ),
            "c": c_val,
            "n_train_groups": len(train_groups),
            "n_train_no_label_groups": n_tr_no_label,
            "n_train_single_positive_groups": n_tr_pos,
            "n_multi_positive_groups_excluded": int(train_multi),
            "n_calibration_groups": len(cal_groups),
            "n_calibration_no_label_groups": n_ca_no_label,
            "n_calibration_single_positive_groups": n_ca_pos,
            "calibration_uniform_marginal_nll": float(uniform_pu_locus_nll(cal_groups, c_val)),
            "calibration_mean_log_n_plus_1": float(plain_uniform_nll(cal_groups)),
        }
    elif n_ca_no_label < min_cal_no_label_groups or n_ca_pos < min_cal_pos_groups:
        locus = {
            "trained": False,
            "reason": (
                f"calibration split has {n_ca_no_label} no-label and {n_ca_pos} single-positive groups "
                f"(< {min_cal_no_label_groups} and < {min_cal_pos_groups} each): cannot validate the null layer"
            ),
            "c": c_val,
            "n_train_groups": len(train_groups),
            "n_train_no_label_groups": n_tr_no_label,
            "n_train_single_positive_groups": n_tr_pos,
            "n_multi_positive_groups_excluded": int(train_multi),
            "n_calibration_groups": len(cal_groups),
            "n_calibration_no_label_groups": n_ca_no_label,
            "n_calibration_single_positive_groups": n_ca_pos,
            "calibration_uniform_marginal_nll": float(uniform_pu_locus_nll(cal_groups, c_val)),
            "calibration_mean_log_n_plus_1": float(plain_uniform_nll(cal_groups)),
        }
    else:
        best = None
        for l2_loc in (0.1, 1.0, 10.0):
            th = fit_pu_locus(train_groups, c_val, l2_loc)
            score = pu_locus_nll(th, cal_groups, c_val, 0.0)
            if best is None or score < best[0]:
                best = (score, l2_loc, th)
        best_score, best_l2, best_th = best
        locus = {
            "trained": True,
            "contract": {
                "id": T.LOCUS_GROUP_ID,
                "fingerprint": T.locus_group_fingerprint(),
                "features": T.LOCUS_GROUP_FEATURES,
            },
            "scale": float(math.exp(float(best_th[0]))),
            "null_bias": float(best_th[1]),
            "null_weights": [float(v) for v in best_th[2:]],
            "hyperparameters": {"l2": float(best_l2), "selected_on": "calibration_marginal_nll"},
            "objective": "pu_marginal_likelihood",
            "c": c_val,
            "n_train_groups": len(train_groups),
            "n_train_no_label_groups": n_tr_no_label,
            "n_train_single_positive_groups": n_tr_pos,
            "n_multi_positive_groups_excluded": int(train_multi),
            "n_calibration_groups": len(cal_groups),
            "n_calibration_no_label_groups": n_ca_no_label,
            "n_calibration_single_positive_groups": n_ca_pos,
            "n_calibration_multi_positive_groups_excluded": int(cal_multi),
            "calibration_marginal_nll": float(best_score),
            "calibration_uniform_marginal_nll": float(uniform_pu_locus_nll(cal_groups, c_val)),
            "calibration_mean_log_n_plus_1": float(plain_uniform_nll(cal_groups)),
            "probability_target": PU_LOCUS_PROBABILITY_TARGET,
        }

    pu_meta = {
        "objective": objective,
        "pi": pi_val,
        "pi_source": pi_source,
        "c": c_val,
        "labelled_rate": labelled_rate,
        "r_neg_train": r_neg_train_in,
        "r_neg_train_unclamped": float(getattr(model, "r_neg_before_clamp", r_neg_train_in)),
        "r_neg_train_oof": r_neg_train_oof,
        "clamp_applied": clamp_applied,
        "explicit_negatives_excluded": int(n_explicit_negatives_excluded),
        "n_positive": int(np.sum(is_pos_tr)),
        "n_unlabelled": int(np.sum(is_unl_tr)),
        "n_calibration_positive": int(np.sum(is_pos_ca)),
        "n_calibration_unlabelled": int(np.sum(is_unl_ca)),
    }
    if kind == "gbdt":
        pu_meta["tree_objective"] = PU_TREE_OBJECTIVE_DOC

    artifact = {
        "schema": T.ARTIFACT_SCHEMA,
        "model_id": f"context-{kind}-{objective}-{data_sha[:12]}-{contract['fingerprint'][:8]}",
        "intended_use": intended_use,
        "contract": {"id": contract["id"], "fingerprint": contract["fingerprint"], "features": names},
        **model.artifact_part(),
        "calibration": calibration,
        "display_threshold": threshold,
        "locus_model": locus,
        "input_support": T.build_input_support(pu_train_rows, pu_calib_rows, names),
        "training": {
            "trainer_version": PU_TRAINER_VERSION,
            "frozen_before_test": True,
            "synthetic": bool(synthetic),
            "training_data_sha256": data_sha,
            "split_identity": identity,
            "weighting": weighting,
            "cv": {
                "folds": folds,
                "seed": seed,
                "unit": "group",
                "objective": objective,
                "pu_risk": float(cv_risk),
                "candidates": [{"params": p, "cv_pu_risk": float(s)} for s, p in all_scores],
            },
            "hyperparameters": params,
            "n_train": int(len(ytr)),
            "n_calibration": int(len(yca)),
            "labelled_rate_train": float(ytr.mean()),
            "pu": pu_meta,
            "truth": {
                "label_basis_counts": {
                    "train": _safe_label_basis_counts(pu_train_rows),
                    "calibration": _safe_label_basis_counts(pu_calib_rows),
                },
                "probability_target": PU_PROBABILITY_TARGET,
                "locus_supervision": (
                    "per-genome function/role competitions under PU marginal likelihood"
                    if not synthetic
                    else "synthetic mechanism competitions under PU marginal likelihood"
                ),
            },
        },
    }
    return artifact
