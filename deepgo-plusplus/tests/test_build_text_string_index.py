"""build_text_string_index.py — one SwissProt pass -> text + STRING + taxon."""
from __future__ import annotations

ENTRY1 = """\
AC   P12345; Q67890;
OX   NCBI_TaxID=9606;
DE   RecName: Full=Test kinase;
GN   Name=tka;
KW   Kinase; Transferase.
CC   -!- FUNCTION: Phosphorylates stuff in cells.
CC       More function detail here.
DR   STRING; 9606.ENSP00000001; -.
//
"""

ENTRY2 = """\
AC   A00002;
OX   NCBI_TaxID=4932;
DE   SubName: Full=Hypothetical protein;
//
"""


def _index(path):
    out = {}
    for line in path.read_text().splitlines():
        ac, tax, sid, name, char = (line.split("\t") + ["", "", "", "", ""])[:5]
        out[ac] = dict(taxon=tax, sid=sid, name=name, char=char)
    return out


def test_parses_acs_taxon_string_and_text_fields(tmp_path, builders):
    dat = builders.write_gz(tmp_path / "sprot.dat.gz", ENTRY1 + ENTRY2)
    out = tmp_path / "index.tsv"
    builders.run_script("build_text_string_index.py", dat, out)
    idx = _index(out)

    # both primary AND secondary accession map to the same entry
    assert "P12345" in idx and "Q67890" in idx
    for ac in ("P12345", "Q67890"):
        e = idx[ac]
        assert e["taxon"] == "9606"
        assert e["sid"] == "9606.ENSP00000001"
        # name_text = identification fields (DE Full + GN Name)
        assert "kinase" in e["name"].lower() and "tka" in e["name"].lower()
        # char_text = CC FUNCTION free text + KW keywords (characterization)
        assert "phosphorylates" in e["char"].lower()
        assert "function detail" in e["char"].lower()
        assert "Kinase" in e["char"]  # keyword

    # an entry with no STRING xref -> empty string_id, still indexed
    e2 = idx["A00002"]
    assert e2["taxon"] == "4932" and e2["sid"] == ""
    assert "hypothetical" in e2["name"].lower()
