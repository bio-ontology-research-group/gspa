"""Tests for prepare_wp4_perturbations.py (WP4 design §4, §6, §10 row 4). Synthetic inputs only."""
import csv
import json
import math
import os
import random
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import prepare_wp4_perturbations as wp

ORTHO_HEADER = "genome_id\tprotein_id\torthogroup_id\tsubfamily_id\tweight\n"
LAYOUT_HEADER = "protein_id\tcontig\tstart\tend\tstrand\n"


def seq(i):
    return "M" + "ACDEFGHIKLMNPQRSTVWY"[i % 20] * 30 + "K" * (i % 7 + 1)


def claim(protein, function_id, function_type="EC"):
    return json.dumps({"evidence_type": "SEQUENCE_SIMILARITY", "function_id": function_id,
                       "function_type": function_type, "protein_id": protein, "raw_score": 0.5,
                       "source": "diamond_homology"}, sort_keys=True) + "\n"


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_ids(path):
    return [line[1:].split()[0] for line in Path(path).read_text().splitlines() if line.startswith(">")]


class Fixture:
    """Two cohort genomes (G1 target, G2), three peers of G1 (P1 congeneric, P2 near-duplicate, P3 congeneric)."""

    def __init__(self, root):
        root = self.root = Path(root)
        self.sequence = root / "sequence"
        self.layout = root / "layout"
        self.claims = root / "target_claims"
        self.peers = root / "peers"
        for d in (self.sequence / "claims", self.layout, self.claims, self.peers):
            d.mkdir(parents=True)
        # ---- cohort genomes
        self.g1 = {f"a{i}": (1000 * i + 100, 1000 * i + 900) for i in range(1, 9)}   # contig C1, 20 kb, closed
        self.g1["a9"] = (4800, 5300)                                                  # crosses a 5 kb cut
        g2 = {f"b{i}": (2000 * i, 2000 * i + 600) for i in range(1, 4)}
        layouts = {"G1": (self.g1, "C1"), "G2": (g2, "D1")}
        fasta, ortho = [], []
        for genome, (proteins, contig) in layouts.items():
            rows = [LAYOUT_HEADER] + [f"{p}\t{contig}\t{s}\t{e}\t+\n" for p, (s, e) in sorted(proteins.items())]
            (self.layout / f"{genome}.layout.tsv").write_text("".join(rows))
            for i, p in enumerate(sorted(proteins)):
                fasta.append(f">{genome}__{p}\n{seq(i)}\n")
                ortho.append(f"{genome}\t{p}\tfam_{i % 3}\t-\t1\n")
        fasta.append(f">G1__u1\n{seq(40)}\n")                       # unplaced protein of G1
        ortho.append("G1\tu1\tfam_9\t-\t1\n")
        (self.layout / "all_proteins.faa").write_text("".join(fasta))
        self.claim_lines = {"G1": [claim("a2", "1.1.1.1"), claim("a2", "GO:0003824", "GO"), claim("a3", "2.2.2.2"),
                                   claim("a5", "1.2.3.-"), claim("a6", "3.3.3.3"), claim("a7", "EC:4.4.4.4"),
                                   claim("a9", "5.5.5.5"), claim("u1", "6.6.6.6")],
                            "G2": [claim("b1", "1.1.1.1")]}
        for genome, lines in self.claim_lines.items():
            (self.claims / f"{genome}.claims.jsonl").write_text("".join(lines))
        metadata = {
            "G1": {"assemblyStatus": "COMPLETE", "clade": "ncbi_genus_1", "strainCluster": "ncbi_species_11",
                   "contigs": {"C1": {"length": 20000, "closed": True}}, "unplacedProteins": ["u1"],
                   "plausibleReactions": ["rxn1"], "familySearchComplete": True},
            "G2": {"assemblyStatus": "COMPLETE", "clade": "ncbi_genus_2", "strainCluster": "ncbi_species_21",
                   "contigs": {"D1": {"length": 9000, "closed": True}}, "unplacedProteins": [],
                   "plausibleReactions": ["rxn2"], "familySearchComplete": True}}
        # ---- peers
        peer_specs = {"P1": ("ncbi_genus_1", "ncbi_species_12"), "P2": ("ncbi_genus_1", "ncbi_species_11"),
                      "P3": ("ncbi_genus_1", "ncbi_species_13")}
        for k, (accession, (clade, strain)) in enumerate(sorted(peer_specs.items())):
            d = self.peers / accession
            d.mkdir()
            contigs = {"NZ1": {"length": 30000, "closed": True}, "NZ2": {"length": 6000, "closed": False}}
            placed = {}
            for i in range(24):                                     # dense and sparse stretches on NZ1
                start = 200 + i * 600 if i < 16 else 12000 + (i - 16) * 2200
                placed[f"{accession}__p{i:02d}"] = ("NZ1", start, start + 450)
            for i in range(5):
                placed[f"{accession}__q{i}"] = ("NZ2", 300 + i * 1100, 300 + i * 1100 + 500)
            unplaced = [f"{accession}__z0"]
            (d / "layout.tsv").write_text(LAYOUT_HEADER + "".join(
                f"{p}\t{c}\t{s}\t{e}\t-\n" for p, (c, s, e) in sorted(placed.items())))
            (d / "unplaced.txt").write_text("".join(p + "\n" for p in unplaced))
            (d / "contigs.json").write_text(json.dumps(contigs))
            proteins = sorted(placed) + unplaced
            (d / "protein.faa").write_text("".join(f">{p}\n{seq(i + k)}\n" for i, p in enumerate(proteins)))
            for i, p in enumerate(proteins):
                ortho.append(f"{accession}\t{p}\tfam_{(i * 7 + k) % 5}\t{'sub_' + str(i % 2)}\t{1 if i % 4 else 0.5}\n")
            (self.sequence / "claims" / f"{accession}.claims.jsonl").write_text(
                "".join(claim(p, "1.1.1.1") for p in proteins[:3]))
            metadata[accession] = {"assemblyStatus": "COMPLETE", "clade": clade, "strainCluster": strain,
                                   "contigs": contigs, "unplacedProteins": unplaced, "familySearchComplete": True}
        (self.sequence / "cohort_orthogroups.tsv").write_text(ORTHO_HEADER + "".join(ortho))
        (self.sequence / "cohort_metadata.json").write_text(json.dumps({"genomes": metadata}))
        self.metadata = metadata
        self.peers_tsv = root / "peers.tsv"
        self.peers_tsv.write_text(
            "target_taxid\ttarget_accession\tpeer_accession\tpeer_species_taxid\tpeer_organism\trole\n"
            "1\tG1\tP1\t12\tGenus one\tclose_peer\n1\tG1\tP2\t11\tGenus one\tnear_duplicate\n"
            "1\tG1\tP3\t13\tGenus one\tclose_peer\n")
        self.labels = root / "labels" / "labels.tsv"
        self.labels.parent.mkdir()
        self.labels.write_text(
            "genome_id\tprotein_id\tfunction_id\tevidence\ttier\tswissprot_accession\tmatch_basis\n"
            "G1\ta2\tEC:1.1.1.1\tIDA\texperimental\tP1\taccession_and_sequence\n"
            "G1\ta3\tGO:0005524\tIDA\texperimental\tP2\taccession_and_sequence\n"
            "G1\ta6\tEC:3.3.3.3\tECO:0000305\tcurated_inferential\tP3\taccession_and_sequence\n"
            "G2\tb1\tEC:1.1.1.1\tIDA\texperimental\tP4\taccession_and_sequence\n")
        self.groups = root / "pathway-groups.tsv"
        edges = [(f"pw{p}", f"grp{g}", f"rxn{r}") for p, g, rs in
                 [(1, 1, (1, 2, 3)), (1, 2, (3, 4)), (2, 3, (5, 6, 7, 8)), (2, 4, (1, 9)), (3, 5, (2, 6, 10))]
                 for r in rs]
        self.groups.write_text("pathway_id\tgroup_id\treaction_id\tsource\n" +
                               "".join(f"{p}\t{g}\t{r}\tgapsmith-curation\n" for p, g, r in edges))

    def ortho_rows(self, path=None):
        return read_tsv(path or self.sequence / "cohort_orthogroups.tsv")


class TmpCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR"))
        self.fx = Fixture(self.tmp.name)
        self.out = Path(self.tmp.name) / "out"

    def tearDown(self):
        self.tmp.cleanup()


class ShuffleTest(TmpCase):
    def run_shuffle(self, out=None, **kwargs):
        return wp.run_shuffle(self.fx.sequence, self.fx.layout, self.fx.peers, self.fx.groups, out or self.out,
                              half_width=kwargs.pop("half_width", 1500), **kwargs)

    def test_s1_preserves_family_multisets_and_strata_and_never_touches_targets(self):
        receipt = self.run_shuffle()
        original = (self.fx.sequence / "cohort_orthogroups.tsv").read_text().splitlines(keepends=True)
        before = self.fx.ortho_rows()
        moved_total = 0
        for seed in wp.S1_SEEDS:
            path = self.out / "S1" / f"seed{seed}" / "cohort_orthogroups.tsv"
            lines = path.read_text().splitlines(keepends=True)
            after = read_tsv(path)
            self.assertEqual(len(lines), len(original))
            # cohort (target) rows byte-identical, same row order and keys everywhere
            cohort = {"G1", "G2"}
            for old, new, a, b in zip(original[1:], lines[1:], before, after):
                self.assertEqual((a["genome_id"], a["protein_id"]), (b["genome_id"], b["protein_id"]))
                if a["genome_id"] in cohort:
                    self.assertEqual(old, new)
            for genome in ("P1", "P2", "P3"):
                layout = wp.read_layout(self.fx.peers / genome / "layout.tsv")
                strata = wp.window_strata(layout, self.fx.metadata[genome]["contigs"], 1500)
                payload = lambda r: (r["orthogroup_id"], r["subfamily_id"], r["weight"])
                old = [r for r in before if r["genome_id"] == genome]
                new = [r for r in after if r["genome_id"] == genome]
                self.assertEqual(Counter(map(payload, old)), Counter(map(payload, new)))
                for stratum in set(strata.values()):
                    members = {p for p, s in strata.items() if s == stratum}
                    self.assertEqual(Counter(payload(r) for r in old if r["protein_id"] in members),
                                     Counter(payload(r) for r in new if r["protein_id"] in members))
                # unplaced proteins keep their family
                self.assertEqual([payload(r) for r in old if r["protein_id"].endswith("z0")],
                                 [payload(r) for r in new if r["protein_id"].endswith("z0")])
            moved_total += sum(receipt["counts"]["S1"][f"seed{seed}"]["proteins_moved"].values())
        self.assertGreater(moved_total, 0)
        self.assertEqual(receipt["parameters"]["peers"], ["P1", "P2", "P3"])
        self.assertEqual([v["variant"] for v in read_tsv(self.out / "variants.tsv")],
                         [f"S1:{s}" for s in wp.S1_SEEDS] + [f"S2:{s}" for s in wp.S2_SEEDS])
        # claims and layouts are not rewritten by S1
        self.assertFalse(any(p.suffix == ".jsonl" or p.name == "layout.tsv" for p in self.out.rglob("*")))

    def test_strata_mirror_engine_buckets(self):
        layout = wp.read_layout(self.fx.peers / "P1" / "layout.tsv")
        strata = wp.window_strata(layout, self.fx.metadata["P1"]["contigs"], 1500)
        self.assertEqual(strata["P1__p05"][2:], ("FULL", "closed"))      # closed contig -> FULL
        self.assertEqual(strata["P1__q0"][2:], ("TRUNCATED", "open"))    # open contig, window crosses the end
        self.assertEqual(strata["P1__q2"][2:], ("FULL", "open"))
        self.assertEqual(strata["P1__p05"][1], 5)                          # dense stretch: 600 bp spacing
        # window bounds are inclusive on both sides (GenomeLayout.windowAround)
        rows = [{"protein_id": n, "contig": "X", "start": m, "end": m, "strand": "+"}
                for n, m in (("left", 1000), ("centre", 2500), ("right", 4000), ("out", 4001))]
        edge = wp.window_strata(rows, {"X": {"length": 10000, "closed": True}}, 1500)
        self.assertEqual(edge["centre"][1], 3)
        self.assertEqual(strata["P1__p20"][1], 1)                          # sparse stretch

    def test_seeded_and_deterministic(self):
        self.run_shuffle()
        self.run_shuffle(out=self.out.with_name("again"))
        for rel in ("S1/seed1/cohort_orthogroups.tsv", "S2/seed2/pathway-groups.tsv"):
            self.assertEqual((self.out / rel).read_bytes(), (self.out.with_name("again") / rel).read_bytes())
        self.assertNotEqual((self.out / "S1/seed1/cohort_orthogroups.tsv").read_bytes(),
                            (self.out / "S1/seed2/cohort_orthogroups.tsv").read_bytes())

    def test_s2_is_degree_preserving(self):
        receipt = self.run_shuffle()
        original = read_tsv(self.fx.groups)
        group_size = Counter((r["pathway_id"], r["group_id"]) for r in original)
        reaction_count = Counter(r["reaction_id"] for r in original)
        changed = 0
        for seed in wp.S2_SEEDS:
            rows = read_tsv(self.out / "S2" / f"seed{seed}" / "pathway-groups.tsv")
            self.assertEqual(Counter((r["pathway_id"], r["group_id"]) for r in rows), group_size)
            self.assertEqual(Counter(r["reaction_id"] for r in rows), reaction_count)
            keys = [(r["pathway_id"], r["group_id"], r["reaction_id"]) for r in rows]
            self.assertEqual(len(keys), len(set(keys)))
            changed += receipt["counts"]["S2"][f"seed{seed}"]["edges_changed"]
        self.assertGreater(changed, 0)

    def test_refuses_label_inputs(self):
        labels = self.fx.root / "labels.tsv"
        labels.write_text(self.fx.groups.read_text())
        with self.assertRaisesRegex(wp.PerturbationError, "label-file"):
            wp.run_shuffle(self.fx.sequence, self.fx.layout, self.fx.peers, labels, self.out)
        self.assertFalse(self.out.exists())


class FragmentTest(TmpCase):
    def run_fragment(self, genomes=("G1", "G2"), out=None, medians=(5000,)):
        return wp.run_fragment(self.fx.sequence, self.fx.layout, self.fx.claims, list(genomes), out or self.out,
                               medians=medians)

    def test_cut_contigs_and_boundary_rule(self):
        pieces, used = wp.cut_contigs({"C1": {"length": 12000}}, iter([5000, 5000, 5000, 5000]))
        self.assertEqual(pieces, {"C1": [(1, 5000), (5001, 10000), (10001, 12000)]})
        self.assertEqual(used, [5000, 5000, 5000])
        rows = [{"protein_id": "end_at_cut", "contig": "C1", "start": 4000, "end": 5000, "strand": "+"},
                {"protein_id": "crossing", "contig": "C1", "start": 4800, "end": 5300, "strand": "+"},
                {"protein_id": "after_cut", "contig": "C1", "start": 5001, "end": 5500, "strand": "-"}]
        kept, lost = wp.fragment_layout(rows, pieces)
        self.assertEqual(lost, ["crossing"])
        self.assertEqual([(r["protein_id"], r["contig"], r["start"], r["end"]) for r in kept],
                         [("end_at_cut", "C1__frag_1_5000", 4000, 5000),
                          ("after_cut", "C1__frag_5001_10000", 1, 500)])

    def test_fixed_draws_are_recorded_and_shared_by_targets(self):
        self.run_fragment()
        draws = json.loads((self.out / "fragment_draws.json").read_text())
        self.assertEqual(sorted(draws), ["median=5000;seed=1", "median=5000;seed=2"])
        for seed in (1, 2):
            record = draws[f"median=5000;seed={seed}"]
            rng = random.Random(f"wp4-fragment|median=5000|sigma=1.0|min=500|seed={seed}")
            expected = [max(500, int(round(rng.lognormvariate(math.log(5000), 1.0))))
                        for _ in record["draws"]]
            self.assertEqual(record["draws"], expected)
            self.assertTrue(all(d >= 500 for d in record["draws"]))
            self.assertEqual(set(record["draws_used"]), {"G1", "G2"})
            self.assertEqual(record["rng_seed"], f"wp4-fragment|median=5000|sigma=1.0|min=500|seed={seed}")
        self.assertEqual(wp.FRAGMENT_MEDIANS, (100000, 25000, 5000))
        self.assertEqual((wp.FRAGMENT_SIGMA, wp.FRAGMENT_MIN_BP, wp.FRAGMENT_SEEDS), (1.0, 500, (1, 2)))
        # a second run reproduces every output byte
        self.run_fragment(out=self.out.with_name("again"))
        first = json.loads((self.out / "receipt.json").read_text())["outputs"]
        second = json.loads((self.out.with_name("again") / "receipt.json").read_text())["outputs"]
        self.assertEqual(first, second)

    def test_boundary_cds_dropped_everywhere(self):
        receipt = self.run_fragment(genomes=("G1",))
        ortho = self.fx.ortho_rows(self.out / "cohort_orthogroups.tsv")
        metadata = json.loads((self.out / "cohort_metadata.json").read_text())["genomes"]
        source_rows = {r["protein_id"]: r for r in ortho if r["genome_id"] == "G1"}
        for copy_id, counts in receipt["counts"]["genomes"]["G1"].items():
            lost = set(counts["removed_proteins"])
            base = self.out / "genomes" / copy_id
            proteins = set(read_ids(base / "protein.faa"))
            layout = wp.read_layout(base / "layout.tsv")
            claims = [json.loads(l)["protein_id"] for l in (self.out / "claims" / f"{copy_id}.claims.jsonl")
                      .read_text().splitlines()]
            family = {r["protein_id"]: r for r in ortho if r["genome_id"] == copy_id}
            self.assertEqual(proteins, set(self.fx.g1) - lost | {"u1"})
            self.assertEqual({r["protein_id"] for r in layout}, proteins - {"u1"})
            self.assertFalse(lost & set(claims))
            self.assertEqual(Counter(claims), Counter(p for p in
                             (json.loads(l)["protein_id"] for l in self.fx.claim_lines["G1"]) if p not in lost))
            self.assertEqual(set(family), proteins)
            for p, row in family.items():
                self.assertEqual(row["orthogroup_id"], source_rows[p]["orthogroup_id"])
            contigs = json.loads((base / "contigs.json").read_text())
            for row in layout:
                self.assertLessEqual(row["end"], contigs[row["contig"]]["length"])
            meta = metadata[copy_id]
            self.assertEqual(meta["assemblyStatus"], "DRAFT")
            self.assertTrue(all(c["closed"] is False for c in meta["contigs"].values()))
            self.assertEqual((meta["clade"], meta["strainCluster"]), ("ncbi_genus_1", "ncbi_species_11"))
            self.assertNotIn("plausibleReactions", meta)
            self.assertEqual(meta["unplacedProteins"], ["u1"])
        self.assertTrue(any(c["boundary_cds_removed"] for c in receipt["counts"]["genomes"]["G1"].values()))
        # the original target rows stay as they were
        self.assertEqual(metadata["G1"], self.fx.metadata["G1"])
        variants = {r["variant"]: r for r in read_tsv(self.out / "variants.tsv")}
        self.assertEqual(sorted(variants), ["fragment:5k:1", "fragment:5k:2"])
        self.assertEqual(variants["fragment:5k:1"]["sealed_genomes"], "G1,G1-frag5k-s1=G1")
        listed = read_tsv(self.out / "gapsmith_inputs.tsv")
        self.assertEqual({r["genome_id"] for r in listed}, set(receipt["counts"]["genomes"]["G1"]))
        for r in listed:
            self.assertEqual(r["sha256"], wp.sha256_file(self.out / r["proteome"]))

    def test_refuses_label_paths_and_product_names(self):
        with self.assertRaisesRegex(wp.PerturbationError, "label-file"):
            wp.run_fragment(self.fx.sequence, self.fx.layout, self.fx.root / "labels" / "claims", ["G1"], self.out)
        path = self.fx.claims / "G1.claims.jsonl"
        obj = json.loads(self.fx.claim_lines["G1"][0])
        obj["metadata"] = {"product": "alcohol dehydrogenase"}
        path.write_text(json.dumps(obj) + "\n")
        with self.assertRaisesRegex(wp.PerturbationError, "product-name"):
            self.run_fragment(genomes=("G1",), out=self.out.with_name("p"))


class RemovalTest(TmpCase):
    def test_removal_selection_copies_and_label_exception(self):
        receipt = wp.run_removal(self.fx.sequence, self.fx.layout, self.fx.claims, self.fx.labels, ["G1"], self.out)
        selection = read_tsv(self.out / "removal_selection.tsv")
        self.assertEqual(list(selection[0]), wp.REMOVAL_FIELDS)
        labelled = [r["removed_protein"] for r in selection if r["class"] == "labelled"]
        randoms = [r["removed_protein"] for r in selection if r["class"] == "random"]
        self.assertEqual(labelled, ["a2"])
        by_protein = {r["removed_protein"]: r for r in selection}
        self.assertEqual((by_protein["a2"]["function_id"], by_protein["a2"]["function_source"]),
                         ("EC:1.1.1.1", "label"))
        claimed = {"a7": "EC:4.4.4.4", "a9": "EC:5.5.5.5", "u1": "EC:6.6.6.6"}
        self.assertEqual((by_protein[randoms[0]]["function_id"], by_protein[randoms[0]]["function_source"]),
                         (claimed[randoms[0]], "claim"))
        # pool: complete-EC claimed proteins without any label row: a7, a9, u1 (a3 GO-labelled, a5 incomplete EC,
        # a6 curated label)
        self.assertEqual(receipt["counts"]["genomes"]["G1"]["unlabelled_ec_claimed_pool"], 3)
        self.assertEqual(randoms, random.Random("42\tG1").sample(["a7", "a9", "u1"], 1))
        ortho = self.fx.ortho_rows(self.out / "cohort_orthogroups.tsv")
        metadata = json.loads((self.out / "cohort_metadata.json").read_text())["genomes"]
        all_g1 = set(self.fx.g1) | {"u1"}
        for row in selection:
            copy_id, removed = row["copy_genome"], row["removed_protein"]
            self.assertEqual(row["target_genome"], "G1")
            self.assertNotIn(removed, copy_id)
            base = self.out / "genomes" / copy_id
            self.assertEqual(set(read_ids(base / "protein.faa")), all_g1 - {removed})
            self.assertNotIn(removed, {r["protein_id"] for r in wp.read_layout(base / "layout.tsv")})
            claims = [json.loads(l)["protein_id"] for l in (self.out / "claims" / f"{copy_id}.claims.jsonl")
                      .read_text().splitlines()]
            self.assertNotIn(removed, claims)
            self.assertEqual(len(claims), sum(1 for l in self.fx.claim_lines["G1"]
                                              if json.loads(l)["protein_id"] != removed))
            self.assertEqual({r["protein_id"] for r in ortho if r["genome_id"] == copy_id}, all_g1 - {removed})
            meta = metadata[copy_id]
            self.assertEqual((meta["clade"], meta["strainCluster"]), ("ncbi_genus_1", "ncbi_species_11"))
            self.assertNotIn("plausibleReactions", meta)
            self.assertEqual(meta["unplacedProteins"], [] if removed == "u1" else ["u1"])
        # cohort and peer rows keep their bytes as the prefix of the joint file
        source = (self.fx.sequence / "cohort_orthogroups.tsv").read_bytes()
        self.assertTrue((self.out / "cohort_orthogroups.tsv").read_bytes().startswith(source))
        # the labels file is a declared label input, separated from the label-blind inputs
        labels = str(self.fx.labels.resolve())
        self.assertEqual(receipt["label_inputs"], {labels: wp.sha256_file(self.fx.labels)})
        self.assertNotIn(labels, receipt["inputs"])
        self.assertIn("documented exception", receipt["label_exception"])
        # label-derived content appears only in the removals manifest (and the receipt that hashes it)
        for path in self.out.rglob("*"):
            if path.is_file() and path.name not in ("receipt.json", "removal_selection.tsv"):
                self.assertNotIn("EC:1.1.1.1", path.read_text())
        variants = read_tsv(self.out / "variants.tsv")
        copies = sorted(r["copy_genome"] for r in selection)
        self.assertEqual(variants, [{"variant": "removal", "target": "G1", "path": ".",
                                     "sealed_genomes": ",".join(["G1"] + [f"{c}=G1" for c in copies])}])
        self.assertEqual(receipt["variants"], variants)
        self.assertEqual(receipt["mode"], "removal")

    def test_pool_too_small_is_refused(self):
        labels = self.fx.root / "more" / "labels.tsv"
        labels.parent.mkdir()
        labels.write_text(self.fx.labels.read_text() + "".join(
            f"G1\t{p}\tEC:9.9.9.{i}\tIDA\texperimental\tQ{i}\taccession_and_sequence\n"
            for i, p in enumerate(["a1", "a4", "a8"])))
        with self.assertRaisesRegex(wp.PerturbationError, "unlabelled EC-claimed"):
            wp.run_removal(self.fx.sequence, self.fx.layout, self.fx.claims, labels, ["G1"], self.out)


class DuplicateTest(TmpCase):
    def test_variants_share_cluster_and_clade_as_designed(self):
        receipt = wp.run_duplicate(self.fx.sequence, self.fx.layout, self.fx.claims, self.fx.peers,
                                   self.fx.peers_tsv, ["G1"], self.out)
        summary = receipt["counts"]["targets"]["G1"]
        self.assertEqual(summary["tiers"], {"P1": "C", "P2": "ND", "P3": "C"})
        self.assertEqual((summary["top_c_peer"], summary["near_duplicate"]), ("P1", "P2"))

        def meta(variant):
            return json.loads((self.out / "G1" / variant / "cohort_metadata.json").read_text())["genomes"]

        m = meta("i_top_c_copy_same_cluster")["P1-dup1"]
        self.assertEqual((m["clade"], m["strainCluster"]), ("ncbi_genus_1", "ncbi_species_12"))
        m = meta("ii_top_c_copy_new_species")["P1-dup2"]
        self.assertEqual(m["clade"], "ncbi_genus_1")
        self.assertNotIn(m["strainCluster"], {g["strainCluster"] for g in self.fx.metadata.values()})
        m = meta("iii_nd_own_cluster")
        self.assertNotEqual(m["P2"]["strainCluster"], "ncbi_species_11")
        self.assertEqual(m["P2"]["clade"], "ncbi_genus_1")
        self.assertEqual(m["P1"], self.fx.metadata["P1"])
        m = meta("iv_target_copy_congeneric")
        self.assertEqual(m["G1-dup4"]["clade"], "ncbi_genus_1")
        self.assertNotEqual(m["G1-dup4"]["strainCluster"], "ncbi_species_11")
        self.assertEqual(m["G1"], self.fx.metadata["G1"])
        self.assertNotIn("P2", meta("v_without_nd"))
        self.assertEqual((self.out / "G1" / "v_without_nd" / "population.txt").read_text().split(),
                         ["G1", "G2", "P1", "P3"])
        # exact copies: same protein sequences, layout, claims and family payloads under the new genome id
        copy_dir = self.out / "genomes" / "G1" / "i_top_c_copy_same_cluster" / "P1-dup1"
        self.assertEqual(sorted((copy_dir / "protein.faa").read_text().splitlines()),
                         sorted((self.fx.peers / "P1" / "protein.faa").read_text().splitlines()))
        self.assertEqual(wp.read_layout(copy_dir / "layout.tsv"), wp.read_layout(self.fx.peers / "P1" / "layout.tsv"))
        claims = self.out / "claims" / "G1" / "i_top_c_copy_same_cluster" / "P1-dup1.claims.jsonl"
        self.assertEqual(claims.read_bytes(), (self.fx.sequence / "claims" / "P1.claims.jsonl").read_bytes())
        rows = self.fx.ortho_rows(self.out / "G1" / "i_top_c_copy_same_cluster" / "cohort_orthogroups.tsv")
        strip = lambda r: (r["protein_id"], r["orthogroup_id"], r["subfamily_id"], r["weight"])
        self.assertEqual([strip(r) for r in rows if r["genome_id"] == "P1-dup1"],
                         [strip(r) for r in rows if r["genome_id"] == "P1"])
        variants = {r["variant"]: r for r in read_tsv(self.out / "variants.tsv")}
        self.assertEqual(sorted(variants), ["dup:i", "dup:ii", "dup:iii", "dup:iv", "dup:v"])
        self.assertEqual(variants["dup:iv"]["sealed_genomes"], "G1,G1-dup4=G1")
        self.assertEqual(variants["dup:i"]["sealed_genomes"], "G1")
        self.assertEqual(variants["dup:v"]["path"], "G1/v_without_nd")
        target_copy = self.out / "genomes" / "G1" / "iv_target_copy_congeneric" / "G1-dup4"
        self.assertEqual(set(read_ids(target_copy / "protein.faa")), set(self.fx.g1) | {"u1"})

    def test_refuses_label_paths(self):
        bad = self.fx.root / "role_truth.tsv"
        bad.write_text(self.fx.peers_tsv.read_text())
        with self.assertRaisesRegex(wp.PerturbationError, "label-file"):
            wp.run_duplicate(self.fx.sequence, self.fx.layout, self.fx.claims, self.fx.peers, bad, ["G1"], self.out)


class ReceiptTest(TmpCase):
    def test_receipt_hashes_and_verify(self):
        receipt = wp.run_fragment(self.fx.sequence, self.fx.layout, self.fx.claims, ["G1"], self.out,
                                  medians=(5000,))
        self.assertEqual(receipt["schema"], wp.SCHEMA)
        self.assertEqual(receipt["script"]["sha256"], wp.sha256_file(wp.__file__))
        for path, digest in receipt["inputs"].items():
            self.assertEqual(wp.sha256_file(path), digest)
        self.assertIn(str((self.fx.sequence / "cohort_orthogroups.tsv").resolve()), receipt["inputs"])
        self.assertIn(str((self.fx.claims / "G1.claims.jsonl").resolve()), receipt["inputs"])
        for rel, digest in receipt["outputs"].items():
            self.assertEqual(wp.sha256_file(self.out / rel), digest)
        self.assertEqual(receipt["label_inputs"], {})
        self.assertTrue(wp.verify(self.out)["passed"])
        layout = next((self.out / "genomes").glob("*/layout.tsv"))
        layout.write_text(layout.read_text() + "x\tC1\t1\t2\t+\n")
        result = wp.verify(self.out)
        self.assertFalse(result["passed"])
        self.assertIn("outputs differ from the receipt", result["failures"])

    def test_verify_detects_changed_input(self):
        wp.run_shuffle(self.fx.sequence, self.fx.layout, self.fx.peers, self.fx.groups, self.out, s1_seeds=(1,),
                       s2_seeds=(1,))
        self.fx.groups.write_text(self.fx.groups.read_text() + "pw9\tgrp9\trxn99\tx\n")
        result = wp.verify(self.out)
        self.assertFalse(result["passed"])
        self.assertIn(str(self.fx.groups.resolve()), result["failures"])

    def test_outputs_have_bare_headers(self):
        (self.fx.layout / "all_proteins.faa").write_text(
            (self.fx.layout / "all_proteins.faa").read_text().replace(">G1__a1\n", ">G1__a1 alcohol dehydrogenase\n"))
        wp.run_fragment(self.fx.sequence, self.fx.layout, self.fx.claims, ["G1"], self.out, medians=(5000,))
        for path in self.out.rglob("protein.faa"):
            for line in path.read_text().splitlines():
                if line.startswith(">"):
                    self.assertRegex(line, r"^>\S+$")
            self.assertNotIn("dehydrogenase", path.read_text())


if __name__ == "__main__":
    unittest.main()
