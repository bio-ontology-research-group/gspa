"""build_lit_component.py — BM25 text-kNN with the temporal leak guard.

The decisive property: the QUERY is matched on its NAME fields only, never on
its (post-t0) characterization text. The fixture is built so a query overlaps
corpus protein A by *name* and corpus protein B by *characterization only* — a
correct implementation must transfer A's label and never B's.
"""
from __future__ import annotations


# Filler corpus proteins with distinct vocabulary so the discriminating tokens
# ('alpha'/'kinase', 'zzznovel'/'marker') are RARE -> positive BM25 idf (a term
# present in every document gets non-positive idf and never transfers).
_FILLERS_INDEX = "".join(
    f"C{i}\t1\t\tword{i}a word{i}b\tdesc{i}x desc{i}y\n" for i in range(6))
_FILLERS_TERMS = "".join(f"C{i}\tGO:0099{i:03d}\tF\n" for i in range(6))


def test_query_uses_name_only_not_characterization(tmp_path, builders):
    # accession <tab> taxon <tab> string_id <tab> name_text <tab> char_text
    builders.write(tmp_path / "index.tsv", _FILLERS_INDEX + "".join([
        "A\t1\t\talpha kinase\tbinds atp\n",
        "B\t1\t\tbeta phosphatase\tzzznovel marker\n",
        # query: name overlaps A; the 'zzznovel marker' overlap with B is in
        # CHAR text, which must be ignored for the query.
        "Q\t1\t\talpha kinase enzyme\tzzznovel marker\n",
    ]))
    builders.write(tmp_path / "train_terms.tsv",
                   "EntryID\tterm\taspect\n" + _FILLERS_TERMS + "".join([
                       "A\tGO:0010001\tF\n",
                       "B\tGO:0020002\tP\n",
                   ]))
    builders.write(tmp_path / "queries.txt", "Q\n")
    out = tmp_path / "lit.tsv"
    builders.run_script(
        "build_lit_component.py", "--index", tmp_path / "index.tsv",
        "--train-terms", tmp_path / "train_terms.tsv",
        "--queries", tmp_path / "queries.txt", "--out", out, "--topk", "30")

    rows = {tuple(l.split("\t")[:2]) for l in out.read_text().splitlines() if l}
    assert ("Q", "GO:0010001") in rows         # transferred via NAME match to A
    assert ("Q", "GO:0020002") not in rows      # B matched only on CHAR -> blocked


def test_sharding_partitions_queries(tmp_path, builders):
    # Two queries, two shards -> each shard handles exactly one (no overlap).
    builders.write(tmp_path / "index.tsv", _FILLERS_INDEX + "".join([
        "A\t1\t\talpha kinase\tbinds atp\n",
        "Q0\t1\t\talpha kinase\t\n",
        "Q1\t1\t\talpha kinase\t\n",
    ]))
    builders.write(tmp_path / "train_terms.tsv",
                   "EntryID\tterm\taspect\n" + _FILLERS_TERMS + "A\tGO:0010001\tF\n")
    builders.write(tmp_path / "queries.txt", "Q0\nQ1\n")
    seen = set()
    for i in range(2):
        out = tmp_path / f"lit_{i}.tsv"
        builders.run_script(
            "build_lit_component.py", "--index", tmp_path / "index.tsv",
            "--train-terms", tmp_path / "train_terms.tsv",
            "--queries", tmp_path / "queries.txt", "--out", out,
            "--shard", f"{i}/2")
        seen |= {l.split("\t")[0] for l in out.read_text().splitlines() if l}
    assert seen == {"Q0", "Q1"}   # union of shards covers every query
