"""Shared style, palette and row selection for the manuscript figures.

Every plotted value is read from ``results/*.tsv`` through :func:`one` / :func:`rows`
with an explicit column=value selector (the same selectors as
``paper-reproducibility/tables/table_sources.tsv``).  Each selection is logged and
``SOURCES`` is written next to the figures by ``generate_figures.py`` so a reader can
trace every mark back to a file, row and column.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "results"
FIG_DIR = REPO / "paper" / "figures"

# ---------------------------------------------------------------------------
# Palette (dataviz reference palette, light mode, validated on white with
# scripts/validate_palette.js: categorical slots 1-4 pass adjacent CVD / normal-vision
# floors; the blue ordinal ramp passes the one-hue ordinal checks).
#
# Visual grammar used in every figure:
#   * GSPA refinement stages: one-hue blue ramp, light -> dark with more refinement
#   * unconstrained / no-prior / reference baselines: neutral grey
#   * external tools and controls: non-blue categorical slots (orange, aqua, yellow, magenta)
#   * paired experimental conditions (design v1/v2, decoy designs, capped/uncapped,
#     truth sets): blue (primary) and orange (stress / alternative)
# ---------------------------------------------------------------------------
BLUE1 = "#86b6ef"   # ramp 250 - first refinement step
BLUE2 = "#5598e7"   # ramp 350
BLUE3 = "#2a78d6"   # ramp 450 - categorical slot 1
BLUE4 = "#184f95"   # ramp 600 - full GSPA
GREY = "#9a9893"    # neutral baseline
GREY_D = "#6f6e69"  # darker neutral for markers on grey
ORANGE = "#eb6834"  # slot 2
AQUA = "#1baf7a"    # slot 3
YELLOW = "#eda100"  # slot 4
MAGENTA = "#e87ba4" # slot 5
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

SEQ_ORANGE = LinearSegmentedColormap.from_list("seq_orange", ["#fdf1ea", "#f6b897", ORANGE, "#a23d14"])

# Figures are placed at 0.96\textwidth (6.86 in) in a two-column article; the canvas
# keeps the previous figure size (about 11.8 x 8.3 in), so 12 pt here prints near 7 pt.
FIGSIZE = (11.8, 8.3)
BASE = 12

STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": BASE,
    "axes.titlesize": BASE + 0.5,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.titlepad": 8,
    "axes.labelsize": BASE,
    "axes.labelcolor": INK,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 1.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.labelsize": BASE - 1,
    "ytick.labelsize": BASE - 1,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "legend.fontsize": BASE - 1.5,
    "legend.frameon": False,
    "legend.handlelength": 1.2,
    "legend.handletextpad": 0.5,
    "legend.borderaxespad": 0.2,
    "legend.columnspacing": 1.2,
    "lines.linewidth": 2.0,
    "lines.markersize": 8,
    "patch.linewidth": 0,
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.default": "regular",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Sans",
    "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold",
}


def apply_style():
    plt.rcParams.update(STYLE)


# ---------------------------------------------------------------------------
# Row selection with provenance
# ---------------------------------------------------------------------------
_CACHE: dict = {}
SOURCES: list = []   # (figure, panel, file, selector, column, value)
_CURRENT = {"fig": "", "panel": ""}


def panel(fig: str, p: str):
    _CURRENT["fig"], _CURRENT["panel"] = fig, p


def load(name: str) -> list:
    if name not in _CACHE:
        with open(RESULTS / name, newline="") as fh:
            _CACHE[name] = list(csv.DictReader(fh, delimiter="\t"))
    return _CACHE[name]


def prefix(p: str):
    """Predicate: method label equals p or is p followed by a ' (' / ' [' suffix."""
    f = lambda v: v is not None and (v == p or v.startswith(p + " (") or v.startswith(p + " ["))
    f.desc = f"'{p}*'"
    return f


def startswith(p: str):
    f = lambda v: v is not None and v.startswith(p)
    f.desc = f"'{p}*'"
    return f


def real_input(v):
    return v is not None and (v == "real" or v.startswith("jar:"))


real_input.desc = "'real' or 'jar:*'"


def _match(r, crit):
    for k, v in crit.items():
        rv = r.get(k)
        if callable(v):
            if not v(rv):
                return False
        elif rv != v:
            return False
    return True


def crit_str(crit: dict) -> str:
    return "; ".join(f"{k}~{getattr(v, 'desc', 'predicate')}" if callable(v) else f"{k}={v}"
                     for k, v in crit.items())


def rows(name: str, crit: dict) -> list:
    return [r for r in load(name) if _match(r, crit)]


def one(name: str, crit: dict) -> dict:
    rs = rows(name, crit)
    if len(rs) != 1:
        raise SystemExit(f"[figures] {name}: {len(rs)} rows match {crit_str(crit)} (need exactly 1)")
    return rs[0]


def num(s):
    if s is None:
        return None
    s = str(s).strip()
    if s in ("", "n/a", "NA", "nan"):
        return None
    mo = re.match(r"^\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", s)
    return float(mo.group(1)) if mo else None


def pm(s):
    """'0.6177+/-0.0866' -> (0.6177, 0.0866)."""
    mo = re.match(r"^\s*([+-]?\d+(?:\.\d+)?)\s*\+/-\s*(\d+(?:\.\d+)?)", s or "")
    return (float(mo.group(1)), float(mo.group(2))) if mo else (num(s), None)


def kv(s, key):
    """'pooled=5.178;median=5.702' -> value of key."""
    for part in (s or "").split(";"):
        if part.strip().startswith(key + "="):
            return float(part.split("=", 1)[1])
    return None


def get(name: str, crit: dict, col: str, derivation: str = ""):
    r = one(name, crit)
    v = num(r.get(col))
    SOURCES.append((_CURRENT["fig"], _CURRENT["panel"], f"results/{name}", crit_str(crit), col, derivation,
                    "" if v is None else v))
    return v


def get_raw(name: str, crit: dict, col: str, derivation: str = ""):
    r = one(name, crit)
    v = r.get(col)
    SOURCES.append((_CURRENT["fig"], _CURRENT["panel"], f"results/{name}", crit_str(crit), col, derivation, v))
    return v


def note(name: str, crit: str, col: str, derivation: str, value):
    SOURCES.append((_CURRENT["fig"], _CURRENT["panel"], f"results/{name}", crit, col, derivation, value))


def write_sources(path: Path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["figure", "panel", "source_file", "row_selector", "source_column", "derivation", "value"])
        for r in SOURCES:
            w.writerow(r)


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def label(ax, letter, x=None, y=None, fig=None):
    """Panel letter, lower-case bold, placed at the top-left of the axes' tight box."""
    fig = fig or ax.figure
    bb = ax.get_position()
    fig.text(bb.x0 - 0.055 if x is None else x, bb.y1 + 0.03 if y is None else y, letter,
             fontsize=BASE + 5, fontweight="bold", color=INK, va="bottom", ha="left")


def figure():
    """2 x 2 grid of sub-figures with constrained layout; returns (fig, [a, b, c, d])."""
    fig = plt.figure(figsize=FIGSIZE, layout="constrained")
    fig.get_layout_engine().set(w_pad=0.08, h_pad=0.08, wspace=0.04, hspace=0.06)
    sf = fig.subfigures(2, 2, wspace=0.04, hspace=0.05)
    return fig, [sf[0, 0], sf[0, 1], sf[1, 0], sf[1, 1]]


def letters(sfs, letters_="abcd"):
    for s, l in zip(sfs, letters_):
        s.text(0.0, 1.003, l, fontsize=BASE + 5, fontweight="bold", color=INK, va="top", ha="left")


def title(ax, text, pad=6):
    """Panel title in the sub-figure, left-aligned next to the panel letter."""
    ax.figure.suptitle(text, x=0.045, y=0.995, ha="left", va="top", fontsize=BASE + 0.5, fontweight="bold",
                       color=INK)
    ax.set_title(" ", pad=pad)


def top_legend(ax, ncol=2, y=1.0, x=0.03, **kw):
    """Legend above the axes, left-aligned with the panel (sub-figure) edge, outside the layout."""
    from matplotlib.transforms import blended_transform_factory
    tr = blended_transform_factory(ax.figure.transSubfigure, ax.transAxes)
    kw.setdefault("fontsize", BASE - 1.5)
    lg = ax.legend(loc="lower left", bbox_to_anchor=(x, y), bbox_transform=tr, ncol=ncol, **kw)
    lg.set_in_layout(False)
    return lg


def hgrid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def save(fig, stem: str):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{stem}.pdf")
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=300)
    plt.close(fig)
    print(f"[figures] wrote paper/figures/{stem}.pdf and .png")


def text_color_for(rgba):
    r, g, b = rgba[:3]
    return INK if (0.299 * r + 0.587 * g + 0.114 * b) > 0.55 else "white"


GENOME_NAMES = {
    "ecoli": r"$\it{E.\ coli}$ K-12",
    "ecolo157": r"$\it{E.\ coli}$ O157:H7",
    "bsubtilis": r"$\it{B.\ subtilis}$ 168",
    "paeruginosa": r"$\it{P.\ aeruginosa}$ PAO1",
    "mtb": r"$\it{M.\ tuberculosis}$ H37Rv",
    "hpylori": r"$\it{H.\ pylori}$ 26695",
    "mgenitalium": r"$\it{M.\ genitalium}$ G37",
    "synechocystis": r"$\it{Synechocystis}$ sp. PCC 6803",
    "mjannaschii": r"$\it{M.\ jannaschii}$ DSM 2661",
    "scerevisiae": r"$\it{S.\ cerevisiae}$",
    "celegans": r"$\it{C.\ elegans}$",
    "dmelanogaster": r"$\it{D.\ melanogaster}$",
    "ppatens": r"$\it{P.\ patens}$",
    "athaliana": r"$\it{A.\ thaliana}$",
    "mmusculus": r"$\it{M.\ musculus}$",
}
CALIBRATION_GENOME = "ecoli"


def genome_label(g, star=True):
    name = GENOME_NAMES.get(g, g)
    return name + ("*" if star and g == CALIBRATION_GENOME else "")
