"""Figure 5: Stage 3 selective guarded reasoning (benchmark v2).

a  results/pillar4_v2_evaluation.tsv (configuration rows of tab:agent) + pillar4_v2_multiseed_summary.tsv MEAN_SD
b  results/pillar4_v2_evaluation.tsv violations_accepted (guarded agent: max over the three seed rows)
c  results/pillar4_v2_pg_tp_hit1.tsv hit1_pct (20 true pathway gaps);
   results/pillar4_v2_ml_baselines.tsv lro50_hit1_pct; pillar4_v2_multiseed_summary.tsv lro_hit1_pct
d  results/pillar4_v2_trigger_rates.tsv genome=SUMMARY_<domain>; pct / pct_uncapped (pooled=...)
"""

from __future__ import annotations

import numpy as np

import common as C

EV = "pillar4_v2_evaluation.tsv"
SEED = "pillar4_v2_multiseed_summary.tsv"
PG = "pillar4_v2_pg_tp_hit1.tsv"
ML = "pillar4_v2_ml_baselines.tsv"
TR = "pillar4_v2_trigger_rates.tsv"
SEEDS = [f"Stage 3 Agent — Qwen 27B + Guardrail (seed={s})" for s in (42, 43, 44)]
GUARDED = "Stage 3 multi-turn + guardrail"


def best_supervised(features):
    """Same rule as tab:agent: highest F1, then accuracy, among '<features>) + Guardrail' rows."""
    cands = [r for r in C.load(EV) if r["configuration"].endswith(f"{features}) + Guardrail")
             and not r["configuration"].startswith("Single Feature")]
    return max(cands, key=lambda r: (float(r["f1"]), float(r["accuracy"])))["configuration"]


def configs():
    rows = [
        ("Stage 1 + Stage 2 Deterministic Only", "Deterministic Stages 1 + 2 only"),
        ("Direct Zero-Shot LLM (No Tools, No Guardrail)", "Direct zero-shot LLM (no tools, no guardrail)"),
        ("Posterior Threshold (>=0.35) + Guardrail", "Posterior threshold (≥ 0.35) + guardrail"),
    ]
    for feats, flab in [("3 Pre-Tool Scalars", "3 pre-tool scalars"), ("8 Post-Tool Features", "8 post-tool features")]:
        conf = best_supervised(feats)
        name = conf.split(" (")[0]
        rows.append((conf, f"Best supervised, {flab} + guardrail ({name})"))
    rows += [("Stage 3 Agent — Fast 1-Step Mode", "Stage 3 fast 1-step"),
             ("Stage 3 Agent — Qwen 27B Without Guardrail", "Stage 3 multi-turn without guardrail"),
             (None, GUARDED + " (3-seed mean ± s.d.)")]
    return rows


def colour(label):
    if label.startswith(GUARDED):
        return C.BLUE4
    if label.startswith("Stage 3"):
        return C.BLUE1
    return C.GREY


def panel_a(sf):
    C.panel("fig5", "a")
    metrics = [("precision", "Precision"), ("recall", "Recall"), ("f1", "$F_1$"), ("accuracy", "Accuracy")]
    rows = configs()
    axes = sf.subplots(1, len(metrics), sharey=True, gridspec_kw={"wspace": 0.10})
    y = np.arange(len(rows))[::-1]
    ms = C.one(SEED, {"seed_label": "MEAN_SD"})
    for ax, (col, lab) in zip(axes, metrics):
        for yi, (conf, rlab) in zip(y, rows):
            if conf is None:
                v, e = C.pm(C.get_raw(SEED, {"seed_label": "MEAN_SD"}, col, "mean +/- population s.d. over seeds 42-44"))
            else:
                v, e = C.get(EV, {"configuration": conf}, col), None
            ax.barh(yi, v, height=0.68, color=colour(rlab), zorder=2, xerr=e,
                    error_kw=dict(ecolor=C.INK, elinewidth=1.3, capsize=3))
            if v == 0:
                ax.text(0.02, yi, "0", va="center", fontsize=C.BASE - 2.5, color=C.MUTED)
        ax.set_xlim(0, 1.0)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels(["0", "0.5", "1"])
        ax.set_xlabel(lab)
        ax.set_title(" ", pad=10)
        C.hgrid(ax, "x")
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([r[1] for r in rows], fontsize=C.BASE - 2)
    n = int(C.get(EV, {"configuration": rows[0][0]}, "n_anomalies"))
    sf.suptitle(f"Anomaly benchmark ({n} de-identified anomalies)", x=0.045, y=0.995, ha="left", va="top",
                fontsize=C.BASE + 0.5, fontweight="bold")
    return [r[1] for r in rows]


def panel_b(ax):
    C.panel("fig5", "b")
    rows = configs()
    y = np.arange(len(rows))[::-1]
    for yi, (conf, rlab) in zip(y, rows):
        if conf is None:
            vs = [C.get(EV, {"configuration": c}, "violations_accepted") for c in SEEDS]
            v = max(vs)
            C.note(EV, "configuration in the three seed rows", "violations_accepted", "max over seeds", v)
        else:
            v = C.get(EV, {"configuration": conf}, "violations_accepted")
        ax.barh(yi, v, height=0.68, color=C.ORANGE if v > 0 else colour(rlab), zorder=2)
        ax.text(v + 0.04, yi, f"{int(v)}", va="center", fontsize=C.BASE - 1.5, color=C.INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 2.2)
    ax.set_xticks([0, 1, 2])
    ax.set_xlabel("Forbidden terms")
    ax.set_title(" ", pad=10)
    ax.set_ylim(y[-1] - 0.6, y[0] + 0.6)
    sf = ax.figure
    sf.suptitle("Violations", x=0.22, y=0.995, ha="left", va="top", fontsize=C.BASE + 0.5, fontweight="bold")


def panel_c(sf):
    C.panel("fig5", "c")
    axl, axr = sf.subplots(1, 2, sharey=True, gridspec_kw={"wspace": 0.08})
    # (label, pg_tp_hit1 method, ml_baselines method); None = not evaluated on that benchmark;
    # "SEEDS" = guarded agent, three seeds
    rows = [
        ("Uniform random candidate", "uniform random candidate (expected)", None),
        ("Operon-context heuristic",
         "argmax operon context (neighbours with present pathway steps + in_pathway_operon)", None),
        ("Stage 1 suggestion score", "argmax suggestion_score",
         "Single Feature: RLGC Suggestion Score (Rank Candidates by suggestion_score)"),
        ("Stage 1 final log-odds", "argmax final_logodds",
         "Single Feature: Stage 1 Final Log-Odds (Rank Candidates by final_logodds)"),
        ("ESM-2 reaction-centroid similarity", None,
         "Single Feature: ESM-2 Reaction-Centroid Cosine (Rank Candidates by esm2_centroid calibrated_prob)"),
        ("Direct zero-shot LLM", "Direct Zero-Shot LLM (No Tools, No Guardrail)",
         "Direct Zero-Shot LLM (No Tools, No Guardrail)"),
        ("Stage 3 fast 1-step", "Stage 3 Agent — Fast 1-Step Mode", "Stage 3 Agent — Fast 1-Step Mode"),
        ("Stage 3 multi-turn without guardrail", "Stage 3 Agent — Qwen 27B Without Guardrail",
         "Stage 3 Agent — Qwen 27B Without Guardrail"),
        (GUARDED + " (3 seeds)", "SEEDS", "SEEDS"),
    ]
    left, right = [], []
    for lab, pg, ml in rows:
        if pg == "SEEDS":
            vs = [C.get(PG, {"method": s}, "hit1_pct") for s in SEEDS]
            mu = sum(vs) / len(vs)
            sd = (sum((v - mu) ** 2 for v in vs) / len(vs)) ** 0.5
            C.note(PG, "method in the three seed rows", "hit1_pct", "mean and population s.d. over seeds", mu)
            left.append((mu / 100, sd / 100, lab))
        else:
            left.append((C.get(PG, {"method": pg}, "hit1_pct", "/100") / 100 if pg else None, None, lab))
        if ml == "SEEDS":
            mu, sd = C.pm(C.get_raw(SEED, {"seed_label": "MEAN_SD"}, "lro_hit1_pct", "mean +/- s.d.; /100"))
            right.append((mu / 100, sd / 100, lab))
        else:
            right.append((C.get(ML, {"method": ml}, "lro50_hit1_pct", "/100") / 100 if ml else None, None, lab))
    n_tp = int(C.get(PG, {"method": SEEDS[0]}, "n_tp_items"))
    _hbars(axl, left, f"{n_tp} true pathway gaps:\ncorrect gene promoted")
    _hbars(axr, right, "50 ambiguous LRO gaps:\nHit@1")
    axl.set_yticklabels([v[2] for v in left], fontsize=C.BASE - 1.5)
    axr.tick_params(axis="y", length=0)
    C.title(axl, "Choosing the correct gene", pad=8)


def _hbars(ax, vals, xlabel):
    y = np.arange(len(vals))[::-1]
    for yi, (v, e, lab) in zip(y, vals):
        if v is None:
            ax.text(0.02, yi, "not evaluated", va="center", fontsize=C.BASE - 2.5, color=C.MUTED, style="italic")
            continue
        col = C.BLUE4 if lab.startswith(GUARDED) else (C.BLUE1 if lab.startswith("Stage 3") else C.GREY)
        ax.barh(yi, v, height=0.68, color=col, zorder=2, xerr=e, error_kw=dict(ecolor=C.INK, elinewidth=1.3, capsize=3))
        ax.text(v + (e or 0) + 0.02, yi, f"{v:.2f}", va="center", fontsize=C.BASE - 2.5, color=C.INK2)
    ax.set_yticks(y)
    ax.set_xlim(0, 1.2)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_ylim(y[-1] - 0.6, y[0] + 0.6)
    ax.set_xlabel(xlabel, fontsize=C.BASE - 1)
    C.hgrid(ax, "x")


def panel_d(ax):
    C.panel("fig5", "d")
    doms = [("SUMMARY_bacteria", "Bacteria"), ("SUMMARY_archaea", "Archaea"), ("SUMMARY_eukaryote", "Eukaryota"),
            ("SUMMARY_virus", "Viruses"), ("SUMMARY_ALL", "All")]
    x = np.arange(len(doms))
    series = [("pct", "Capped (150 anomalies per proteome)", C.BLUE3), ("pct_uncapped", "Uncapped", C.ORANGE)]
    capped = [r for r in C.load(TR) if not r["genome"].startswith("SUMMARY") and r.get("cap_reached") == "True"]
    cap = int(max(C.num(r["n_anomalies"]) for r in capped))
    C.note(TR, "cap_reached=True (per-genome rows)", "n_anomalies", "max = production anomaly cap", cap)
    for k, (col, lab, color) in enumerate(series):
        vals = [C.kv(C.get_raw(TR, {"genome": g}, col, "pooled= value"), "pooled") for g, _ in doms]
        xs = x + (k - 0.5) * 0.38
        ax.bar(xs, vals, width=0.36, color=color, zorder=2,
               label=lab.replace("150", str(cap)) if k == 0 else lab)
        for xi, v in zip(xs, vals):
            ax.text(xi, v + 0.25, f"{v:.1f}", ha="center", va="bottom", fontsize=C.BASE - 2.5, color=C.INK2)
    labs = []
    for g, lab in doms:
        n = C.get_raw(TR, {"genome": g}, "taxon").split("=")[-1]
        labs.append(f"{lab}\nn = {n}")
    ax.set_xticks(x)
    ax.set_xticklabels(labs, fontsize=C.BASE - 1.5)
    ax.set_xlabel("Domain (n = proteomes)", fontsize=C.BASE - 1)
    ax.set_ylabel("Proteins triggering Stage 3 (%)")
    C.hgrid(ax)
    C.top_legend(ax, ncol=2)
    C.title(ax, "Stage 3 trigger rate by domain", pad=30)


def generate():
    import matplotlib.pyplot as plt
    C.apply_style()
    fig = plt.figure(figsize=C.FIGSIZE, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.08, h_pad=0.08, wspace=0.04, hspace=0.06)
    top, bot = fig.subfigures(2, 1, hspace=0.05, height_ratios=[1, 1.05])
    sfa, sfb = top.subfigures(1, 2, width_ratios=[4.3, 1], wspace=0.0)
    sfc, sfd = bot.subfigures(1, 2, wspace=0.04, width_ratios=[1.15, 0.85])
    panel_a(sfa)
    panel_b(sfb.subplots())
    panel_c(sfc)
    panel_d(sfd.subplots())
    C.letters([sfa, sfb, sfc, sfd])
    C.save(fig, "fig5_pillar4_agent")
