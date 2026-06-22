"""build_net_component.py — Net-KNN guilt-by-association over STRING PPI."""
from __future__ import annotations


def _rows(path):
    out = set()
    for l in path.read_text().splitlines():
        p, t, s = l.split("\t")
        out.add((p, t, float(s)))
    return out


def _build(tmp_path, builders):
    # accession <tab> taxon <tab> string_id
    builders.write(tmp_path / "net_index.tsv", "".join([
        "Q1\t9606\t9606.Q1\n",
        "N1\t9606\t9606.N1\n",
        "N2\t9606\t9606.N2\n",
        "Q2\t4932\t4932.Q2\n",   # in a species whose STRING file is corrupt
    ]))
    builders.write(tmp_path / "train_terms.tsv", "".join([
        "EntryID\tterm\taspect\n",
        "N1\tGO:0010001\tF\n",
        "N1\tGO:0020001\tP\n",
        "N2\tGO:0030001\tC\n",
    ]))
    builders.write(tmp_path / "queries.txt", "Q1\nQ2\n")
    sdir = tmp_path / "string"
    builders.write_gz(sdir / "9606.protein.links.v12.0.txt.gz", "".join([
        "protein1 protein2 combined_score\n",
        "9606.Q1 9606.N1 800\n",   # high confidence -> kept
        "9606.Q1 9606.N2 300\n",   # below min-conf 400 -> dropped
    ]))
    # a corrupt (non-gzip) STRING file with the expected name -> must be skipped
    builders.write(sdir / "4932.protein.links.v12.0.txt.gz", "not gzip at all\n")
    out = tmp_path / "net.tsv"
    builders.run_script(
        "build_net_component.py", "--index", tmp_path / "net_index.tsv",
        "--train-terms", tmp_path / "train_terms.tsv",
        "--queries", tmp_path / "queries.txt", "--string-dir", sdir,
        "--out", out, "--min-conf", "400", "--topk", "50")
    return out


def test_neighbor_vote_and_minconf(tmp_path, builders):
    out = _build(tmp_path, builders)
    rows = _rows(out)
    # Only N1 survives min-conf; its two labels transfer, normalised to 1.0
    assert rows == {("Q1", "GO:0010001", 1.0), ("Q1", "GO:0020001", 1.0)}
    # N2's label is gated out by the confidence cutoff
    assert all(t != "GO:0030001" for _q, t, _s in rows)


def test_corrupt_string_file_skipped(tmp_path, builders):
    out = _build(tmp_path, builders)  # would raise if the corrupt file crashed it
    proteins = {p for p, _t, _s in _rows(out)}
    assert "Q2" not in proteins   # its species was skipped, gracefully
    assert "Q1" in proteins       # the healthy species still produced output
