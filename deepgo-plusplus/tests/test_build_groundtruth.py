"""build_groundtruth.py — CAFA6 knowledge-class split + t0 boundary."""
from __future__ import annotations


def _rows(path):
    return {tuple(l.split("\t")) for l in path.read_text().splitlines() if l}


def test_knowledge_classes_and_t0_boundary(tmp_path, builders):
    # acc, GO, F|P|C, evidence, YYYYMMDD ; t0 = 2026-02-02 (20260202)
    annots = builders.write(tmp_path / "annots.tsv", "".join([
        # no-knowledge: nothing pre-t0, gains an MF term after t0
        "Pno\tGO:0010001\tF\tEXP\t20260301\n",
        # limited: one pre-t0 aspect (F, exactly ON t0 -> counts as pre),
        #          gains a BP term after t0 (BP is the test target aspect)
        "Plim\tGO:0010002\tF\tEXP\t20260202\n",
        "Plim\tGO:0020001\tP\tEXP\t20260301\n",
        # partial: all three aspects pre-t0, gains another MF term after t0
        "Ppart\tGO:0010001\tF\tEXP\t20250101\n",
        "Ppart\tGO:0020001\tP\tEXP\t20250101\n",
        "Ppart\tGO:0030001\tC\tEXP\t20250101\n",
        "Ppart\tGO:0010002\tF\tEXP\t20260301\n",
        # pre-t0 only, never gains -> NOT a test target
        "Ppre\tGO:0010001\tF\tEXP\t20250101\n",
    ]))
    out = tmp_path / "gt"
    builders.run_script("build_groundtruth.py", "--annots", annots, "--outdir", out)

    test_proteins = set(out.joinpath("test_proteins.txt").read_text().split())
    assert test_proteins == {"Pno", "Plim", "Ppart"}  # Ppre excluded

    # no-knowledge: only Pno, its post-t0 MF term (full-eval = all terms in aspect)
    assert _rows(out / "gt_no.tsv") == {("Pno", "GO:0010001")}

    # limited: Plim's BP target aspect terms only (F is pre, not a target aspect)
    assert _rows(out / "gt_limited.tsv") == {("Plim", "GO:0020001")}

    # partial: Ppart's full MF set across all dates (pre + post)
    assert _rows(out / "gt_partial.tsv") == {
        ("Ppart", "GO:0010001"), ("Ppart", "GO:0010002")}

    # gt_all is the union of the three classes
    assert _rows(out / "gt_all.tsv") == (
        _rows(out / "gt_no.tsv") | _rows(out / "gt_limited.tsv")
        | _rows(out / "gt_partial.tsv"))


def test_custom_t0_reclassifies(tmp_path, builders):
    # With a much later t0, the post-t0 gain disappears -> protein is no longer
    # a test target (guards the t0 plumbing, not just the default).
    annots = builders.write(tmp_path / "a.tsv",
                            "Px\tGO:0010001\tF\tEXP\t20260301\n")
    out = tmp_path / "gt"
    builders.run_script("build_groundtruth.py", "--annots", annots,
                        "--t0", "20270101", "--outdir", out)
    assert out.joinpath("test_proteins.txt").read_text().split() == []
