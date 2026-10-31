"""Unit tests for the AllTheBacteria pilot helpers (plain unittest, Python >= 3.8).

  /home/leechuck/tmp/gspa-v3-testenv/bin/python -m unittest discover -s benchmark/allthebacteria/shaheen/tests
"""
import csv
import gzip
import importlib.util
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)


def readf(path):
    with open(path) as f:
        return f.read()


def readjson(path):
    with open(path) as f:
        return json.load(f)


def readjsonl(path):
    with open(path) as f:
        return [json.loads(x) for x in f]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(PKG, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


receipt = load("atb_receipt", "receipt.py")
summ = load("atb_summarise", "summarise_receipts.py")
select = load("atb_select", "00_select_pilot.py")
scatter = load("atb_scatter", "atb_scatter.py")

TIMEV = """\tCommand being timed: "pyrodigal -i x.fna -a x.faa"
\tUser time (seconds): 25.50
\tSystem time (seconds): 1.25
\tPercent of CPU this job got: 98%
\tElapsed (wall clock) time (h:mm:ss or m:ss): {wall}
\tMaximum resident set size (kbytes): 123456
\tExit status: {exit}
"""


class TmpCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="atbp-test-")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def path(self, *p):
        return os.path.join(self.tmp, *p)

    def write(self, rel, text):
        p = self.path(rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as f:
            f.write(text)
        return p


class TestTimev(TmpCase):
    def test_parse_minutes_and_hours(self):
        p = self.write("a.time", TIMEV.format(wall="1:02.50", exit=0))
        r = receipt.parse_timev(p)
        self.assertAlmostEqual(r["wall_s"], 62.5)
        self.assertAlmostEqual(r["user_s"], 25.5)
        self.assertAlmostEqual(r["sys_s"], 1.25)
        self.assertEqual(r["maxrss_kb"], 123456)
        self.assertEqual(r["exit"], 0)
        p = self.write("b.time", "Command exited with non-zero status 3\n" + TIMEV.format(wall="2:00:01", exit=3))
        r = receipt.parse_timev(p)
        self.assertAlmostEqual(r["wall_s"], 7201.0)
        self.assertEqual(r["exit"], 3)

    def test_missing_or_truncated(self):
        self.assertIsNone(receipt.parse_timev(self.path("nope.time"))["wall_s"])
        p = self.write("t.time", "\tCommand being timed: \"x\"\n")
        self.assertIsNone(receipt.parse_timev(p)["wall_s"])

    def test_collect_and_aggregate(self):
        root = self.path("times")
        self.write("times/genecall/pyrodigal/s000.G1.time", TIMEV.format(wall="0:10.00", exit=0))
        self.write("times/genecall/pyrodigal/s000.G2.time", TIMEV.format(wall="0:30.00", exit=1))
        self.write("times/genecall/pyrodigal/s000.G3.time", "")          # killed: incomplete
        self.write("times/linclust/s000.time", TIMEV.format(wall="0:05.00", exit=0))  # no phase dir -> main
        recs = receipt.collect_times(root)
        self.assertEqual(len(recs), 4)
        agg = {(g["phase"], g["tool"]): g for g in receipt.aggregate(recs)}
        g = agg[("genecall", "pyrodigal")]
        self.assertEqual((g["n"], g["n_failed"], g["n_incomplete"]), (3, 1, 1))
        self.assertAlmostEqual(g["cpu_s"], 2 * 26.75)
        self.assertAlmostEqual(g["wall_s_median"], 20.0)
        self.assertEqual(g["wall_s_max"], 30.0)
        self.assertIn(("main", "linclust"), agg)

    def test_phases(self):
        p = self.path("phases.tsv")
        for name, ev in (("a", "start"), ("a", "end"), ("b", "start")):
            receipt.main(["phase", "--file", p, "--name", name, "--event", ev])
        ph = receipt.read_phases(p)
        self.assertEqual([x["name"] for x in ph], ["a", "b"])
        self.assertIsNotNone(ph[0]["end"])
        self.assertIsNone(ph[1]["end"])


def make_receipt(path, step, job_id, wall, phases=(), counts=None, kind="job", nodes=1, cores=192, cpu_s=0.0):
    t0 = 1000.0
    rec = {"schema": receipt.SCHEMA, "step": step, "unit": f"job-{job_id}", "kind": kind, "job_id": str(job_id),
           "nodes": nodes, "cores_per_node": cores, "job_start": t0, "job_end": t0 + wall, "wall_s": wall,
           "phases": [{"name": n, "start": t0 + a, "end": t0 + b} for n, a, b in phases],
           "counts": counts or {}, "resources": [{"phase": phases[0][0] if phases else "main", "tool": "x", "cpu_s": cpu_s}]}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(rec, f)


class TestReceiptWrite(TmpCase):
    def test_write_roundtrip(self):
        self.write("times/genecall/pyrodigal/s000.G1.time", TIMEV.format(wall="0:10.00", exit=0))
        phases = self.path("phases.tsv")
        receipt.main(["phase", "--file", phases, "--name", "genecall", "--event", "start"])
        receipt.main(["phase", "--file", phases, "--name", "genecall", "--event", "end"])
        vers = self.write("versions.tsv", "pyrodigal\tpyrodigal v3.7.1\n")
        inp = self.write("in.txt", "abc")
        out = self.path("receipts", "S1", "1.json")
        tsv = self.path("receipts", "S1", "1.resources.tsv.gz")
        receipt.main(["write", "--step", "S1", "--unit", "job-1", "--out", out, "--phases", phases, "--times-root",
                      self.path("times"), "--job-start", "100", "--job-end", "160", "--count", "n_genomes=2",
                      "--param", "shards=s000:s001", "--input-sha", inp, "--versions-file", vers, "--resources-tsv", tsv])
        r = readjson(out)
        self.assertEqual(r["schema"], receipt.SCHEMA)
        self.assertEqual(r["wall_s"], 60.0)
        self.assertEqual(r["counts"]["n_genomes"], 2)
        self.assertEqual(r["params"]["shards"], "s000:s001")
        self.assertEqual(r["versions"]["pyrodigal"], "pyrodigal v3.7.1")
        self.assertEqual(r["inputs"][0]["sha256"], "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
        self.assertEqual(r["resources"][0]["n"], 1)
        with gzip.open(tsv, "rt") as f:
            rows = f.read().splitlines()
        self.assertEqual(len(rows), 2)
        self.assertTrue(rows[1].startswith("genecall\tpyrodigal\ts000.G1\t10.0"))


class TestSummarise(TmpCase):
    def build(self):
        R = self.path("receipts")
        make_receipt(f"{R}/S0/1.json", "S0", 1, 60, counts={"n_genomes": 1000}, cores=1)
        # S1: 2 jobs of 1 h each, 500 genomes each, 2,000 proteins per genome
        make_receipt(f"{R}/S1/2.json", "S1", 2, 3600, phases=[("extract", 0, 600), ("genecall", 600, 3000)],
                     counts={"n_genomes": 500, "n_proteins": 1_000_000}, cpu_s=192 * 2400)
        make_receipt(f"{R}/S1/3.json", "S1", 3, 3600, counts={"n_genomes": 500, "n_proteins": 1_000_000})
        make_receipt(f"{R}/S2/4.json", "S2", 4, 1800, counts={"n_genomes": 1000, "n_proteins": 2_000_000,
                                                             "n_representatives": 200_000})
        make_receipt(f"{R}/S3/5.json", "S3", 5, 7200, counts={"n_representatives": 200_000})
        for i in range(4):
            make_receipt(f"{R}/S3/chunks/c{i}.json", "S3", 5, 100, counts={"n_representatives": 50_000}, kind="unit", nodes=0)
        make_receipt(f"{R}/S4/6.json", "S4", 6, 360, counts={"n_genomes": 1000})
        make_receipt(f"{R}/S5/7.json", "S5", 7, 3600, counts={"n_genomes": 1000})
        make_receipt(f"{R}/S6/8.json", "S6", 8, 3600, phases=[("stage0", 0, 600), ("stage12", 600, 3000)],
                     counts={"n_genomes": 1000})
        return R

    def test_costs_and_projection(self):
        jobs, units = summ.load_receipts(self.build())
        self.assertEqual(len(units), 4)
        res = summ.summarise(jobs, units, n_target=2_000_000, budget=4.5e6)
        st = {r["step"]: r for r in res["steps"]}
        self.assertEqual(res["pilot_genomes"], 1000)
        self.assertAlmostEqual(res["proteins_per_genome"], 2000)
        self.assertAlmostEqual(res["rep_fraction_measured"], 0.1)
        self.assertAlmostEqual(st["S1"]["charged_core_h"], 2 * 192)
        self.assertAlmostEqual(st["S1"]["charged_core_h_per_genome"], 0.384)
        self.assertAlmostEqual(st["S1"]["used_cpu_h"], 128)
        self.assertAlmostEqual(st["S1"]["projected_core_h"], 0.384 * 2_000_000)
        # S3: 2 node-hours on 200,000 representatives
        self.assertAlmostEqual(st["S3"]["charged_core_h_per_rep"], 384 / 200_000)
        self.assertAlmostEqual(st["S3"]["projected_core_h"], 384 / 200_000 * 2_000_000 * 2000 * 0.1)
        self.assertEqual(st["S0"]["projected_core_h"], 0.0)
        total = sum(r["projected_core_h"] for r in res["steps"])
        self.assertAlmostEqual(res["projected_core_h"], total)
        self.assertEqual(res["verdict"], "FAIL" if total > 4.5e6 else "PASS")
        ph = {(p["step"], p["phase"]): p for p in res["phases"]}
        self.assertAlmostEqual(ph[("S1", "extract")]["charged_core_h"], 192 * 600 / 3600)
        self.assertAlmostEqual(ph[("S1", "overhead")]["charged_core_h"], 2 * 192 - 192 * 3000 / 3600)
        self.assertAlmostEqual(ph[("S6", "stage12")]["charged_core_h"], 192 * 2400 / 3600)

    def test_rep_fraction_overhead_budget(self):
        jobs, units = summ.load_receipts(self.build())
        base = summ.summarise(jobs, units, n_target=2_000_000, budget=1e12)
        low = summ.summarise(jobs, units, n_target=2_000_000, rep_fraction=0.05, budget=1e12)
        s3b = [r for r in base["steps"] if r["step"] == "S3"][0]["projected_core_h"]
        s3l = [r for r in low["steps"] if r["step"] == "S3"][0]["projected_core_h"]
        self.assertAlmostEqual(s3l, s3b / 2)
        over = summ.summarise(jobs, units, n_target=2_000_000, overhead_factor=1.3, budget=1e12)
        self.assertAlmostEqual(over["projected_core_h"], 1.3 * base["projected_core_h"])
        self.assertEqual(base["verdict"], "PASS")
        self.assertEqual(summ.summarise(jobs, units, n_target=2_000_000, budget=1.0)["verdict"], "FAIL")

    def test_incomplete(self):
        R = self.build()
        os.remove(os.path.join(R, "S6", "8.json"))
        jobs, units = summ.load_receipts(R)
        res = summ.summarise(jobs, units)
        self.assertEqual(res["verdict"], "INCOMPLETE")
        self.assertEqual(res["missing_steps"], ["S6"])

    def test_sacct_override_and_unreceipted(self):
        R = self.build()
        sacct = self.write("sacct.tsv", "2|atbp-S1-g0|1|7200|COMPLETED\n"
                                        "9|atbp-S3-w1|1|86400|TIMEOUT\n"
                                        "9.batch|batch|1|86400|CANCELLED\n"
                                        "10|other-job|1|100|COMPLETED\n")
        jobs, units = summ.load_receipts(R)
        res = summ.summarise(jobs, units, summ.load_sacct(sacct))
        st = {r["step"]: r for r in res["steps"]}
        self.assertAlmostEqual(st["S1"]["charged_core_h"], 192 * 2 + 192)       # job 2 now 2 h
        self.assertAlmostEqual(st["S3"]["charged_core_h"], 2 * 192 + 24 * 192)  # + killed worker 9
        self.assertEqual(res["unreceipted_jobs"], ["9"])
        ph = {(p["step"], p["phase"]): p for p in res["phases"]}
        self.assertAlmostEqual(ph[("S3", "unreceipted")]["charged_core_h"], 24 * 192)

    def test_cli_outputs(self):
        R = self.build()
        rc = summ.main([R, "--out-tsv", self.path("s.tsv"), "--out-json", self.path("s.json"), "--out-md", self.path("s.md")])
        self.assertEqual(rc, 0)
        with open(self.path("s.tsv")) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual([r["step"] for r in rows], ["S0", "S1", "S2", "S3", "S4", "S5", "S6"])
        self.assertIn("G0 projection", readf(self.path("s.md")))
        self.assertEqual(readjson(self.path("s.json"))["schema"], "gspa.atb-pilot-g0/1")
        self.assertEqual(summ.main([self.path("empty")]), 2)


class TestSelect(TmpCase):
    def build(self):
        cat = ["batch_id\ttar_xz\ttar_xz_url\ttar_xz_md5\ttar_xz_size_MB\tn_genomes\tdominant_species\tdominant_phylum\tmanifest_path"]
        spec = [("w2_0001", 40, "Escherichia coli", "p__Pseudomonadota", "NCBITaxon:562"),
                ("w2_0002", 30, "Klebsiella pneumoniae", "p__Pseudomonadota", "NCBITaxon:999999"),
                ("w2_0003", 50, "Staphylococcus aureus", "p__Bacillota", ""),
                ("w2_0004", 5, "Tiny bug", "p__Bacillota", ""),
                ("w2_0005", 20, "Weird_A thing", "p__Weirdota", "")]
        os.makedirs(self.path("manifests"))
        for bid, n, sp, ph, tid in spec:
            cat.append(f"{bid}\t{bid}.tar.xz\tu\tmd5{bid}\t1\t{n}\t{sp}\t{ph}\t/ibex/x/{bid}.tsv")
            rows = ["sample\tspecies\tdomain\tphylum\tncbi_taxid\tessential_profile\tsource_type\tsource_uri\ttar_xz\tfilename_in_tar_xz"]
            for i in range(n):
                rows.append(f"{bid}_G{i:03d}\t{sp}\td__Bacteria\t{ph}\t{tid}\tprokaryote\tatb\tu\t{bid}.tar.xz\t{bid}/{bid}_G{i:03d}.fa")
            self.write(f"manifests/{bid}.tsv", "\n".join(rows) + "\n")
        self.write("catalog.tsv", "\n".join(cat) + "\n")
        self.write("labels.tsv", "Taxon_ID\tLabel\nNCBITaxon_562\tEscherichia coli\nNCBITaxon_570\tKlebsiella\n"
                                 "NCBITaxon_1239\tFirmicutes\nNCBITaxon_2\tBacteria\n")

    def run_select(self, out, seed=7, shards=3, per=10):
        select.main(["--catalog", self.path("catalog.tsv"), "--manifests", self.path("manifests"), "--labels",
                     self.path("labels.tsv"), "--seed", str(seed), "--shards", str(shards), "--genomes-per-shard",
                     str(per), "--out-dir", self.path(out)])
        with open(self.path(out, "pilot_genomes.tsv")) as f:
            genomes = list(csv.DictReader(f, delimiter="\t"))
        with open(self.path(out, "shards.tsv")) as f:
            shards_rows = list(csv.DictReader(f, delimiter="\t"))
        return genomes, shards_rows, readjson(self.path(out, "selection.json"))

    def test_selection_reproducible_and_stratified(self):
        self.build()
        g1, s1, sel = self.run_select("a")
        g2, s2, _ = self.run_select("b")
        self.assertEqual(g1, g2)
        self.assertEqual(len(g1), 30)
        self.assertEqual(len(s1), 3)
        self.assertEqual(sel["stratum_shards"], {"p__Pseudomonadota": 2, "p__Bacillota": 1, "p__Weirdota": 0})
        self.assertNotIn("w2_0004", [s["batch_id"] for s in s1])          # too small for a shard
        self.assertEqual(sorted(s["shard"] for s in s1), ["s000", "s001", "s002"])
        for s in s1:
            self.assertEqual(sum(1 for g in g1 if g["shard"] == s["shard"]), 10)
        self.assertEqual(len({g["genome_id"] for g in g1}), 30)
        g3, _, _ = self.run_select("c", seed=8)
        self.assertNotEqual([g["genome_id"] for g in g1], [g["genome_id"] for g in g3])

    def test_taxon_gram_profile(self):
        self.build()
        g, _, _ = self.run_select("a")
        by_batch = {}
        for r in g:
            by_batch.setdefault(r["batch_id"], r)
        self.assertEqual((by_batch["w2_0001"]["taxon"], by_batch["w2_0001"]["taxon_rank"]), ("NCBITaxon_562", "manifest_taxid"))
        self.assertEqual((by_batch["w2_0002"]["taxon"], by_batch["w2_0002"]["taxon_rank"]), ("NCBITaxon_570", "genus"))
        self.assertEqual((by_batch["w2_0003"]["taxon"], by_batch["w2_0003"]["taxon_rank"]), ("NCBITaxon_1239", "phylum"))
        self.assertEqual(by_batch["w2_0003"]["gram"], "pos")
        self.assertEqual(by_batch["w2_0001"]["gram"], "neg")
        self.assertTrue(all(r["essential_profile"] == "bacteria" for r in g))
        self.assertEqual(select.assign_taxon({"species": "s__Nothing known", "phylum": "p__Weirdota", "domain": "d__Bacteria"},
                                             {}, set()), ("NCBITaxon_2", "domain"))

    def test_allocation_helpers(self):
        self.assertEqual(select.largest_remainder({"a": 50, "b": 30, "c": 20}, 3), {"a": 1, "b": 1, "c": 1})
        self.assertEqual(select.largest_remainder({"a": 60, "b": 30, "c": 10}, 3), {"a": 2, "b": 1, "c": 0})
        self.assertEqual(sum(select.largest_remainder({"a": 1, "b": 1, "c": 1}, 2).values()), 2)
        import random
        picks = select.weighted_sample(["x", "y", "z"], [1, 1, 0], 3, random.Random(1))
        self.assertEqual(sorted(picks), ["x", "y"])

    def test_not_enough_batches(self):
        self.build()
        with self.assertRaises(SystemExit):
            self.run_select("d", shards=5, per=25)


class TestScatter(TmpCase):
    def build(self, finish=True):
        D = self.path("s000")
        os.makedirs(os.path.join(D, "chunks"))
        self.write("s000/counts.tsv", "GA\t3\nGB\t2\nGC\t0\n")
        with gzip.open(os.path.join(D, "clu_cluster.tsv.gz"), "wt") as f:
            for rep, mem in (("GA|p1", "GA|p1"), ("GA|p1", "GB|q1"), ("GA|p2", "GA|p2"), ("GA|p3", "GA|p3"), ("GB|q2", "GB|q2")):
                f.write(f"{rep}\t{mem}\n")
        self.write("s000/chunks/c0000.faa", ">GA|p1\nM\n>GA|p2\nM\n")
        self.write("s000/chunks/c0001.faa", ">GA|p3\nM\n>GB|q2\nM\n")
        claims = {"c0000": [{"protein_id": "GA|p1", "term": "GO:1", "source": "diamond", "raw_score": 0.9},
                            {"protein_id": "GA|p1", "term": "GO:2", "source": "pfam", "raw_score": 1.0},
                            {"protein_id": "GX|zz", "term": "GO:3", "source": "pfam", "raw_score": 1.0}],
                  "c0001": [{"protein_id": "GB|q2", "term": "GO:4", "source": "interpro", "raw_score": 0.5}]}
        for c, recs in claims.items():
            with gzip.open(os.path.join(D, "chunks", f"{c}.claims.jsonl.gz"), "wt") as f:
                for r in recs:
                    f.write(json.dumps(r) + "\n")
            if finish or c == "c0000":
                self.write(f"s000/chunks/{c}.done", "x\n")
        return D

    def test_scatter(self):
        D = self.build()
        out = self.path("claims")
        rc = scatter.scatter(D, out, os.path.join(D, "coverage.tsv"), handle_limit=1)
        self.assertEqual(rc, 0)
        read = lambda g: readjsonl(os.path.join(out, f"{g}.claims.jsonl"))
        ga, gb, gc = read("GA"), read("GB"), read("GC")
        self.assertEqual(sorted((r["protein_id"], r["term"]) for r in ga), [("p1", "GO:1"), ("p1", "GO:2")])
        self.assertEqual(sorted((r["protein_id"], r["term"]) for r in gb), [("q1", "GO:1"), ("q1", "GO:2"), ("q2", "GO:4")])
        self.assertEqual(gc, [])
        self.assertEqual(ga[0]["source"], "diamond")
        with open(os.path.join(D, "coverage.tsv")) as f:
            cov = {r["genome_id"]: r for r in csv.DictReader(f, delimiter="\t")}
        self.assertEqual((cov["GA"]["n_proteins_clustered"], cov["GA"]["n_proteins_with_claim"], cov["GA"]["n_claims"]), ("3", "1", "2"))
        self.assertEqual((cov["GB"]["n_proteins_clustered"], cov["GB"]["n_proteins_with_claim"], cov["GB"]["n_claims"]), ("2", "2", "3"))
        self.assertEqual(cov["GC"]["n_claims"], "0")

    def test_protein_id_with_special_characters(self):
        D = self.build()
        with gzip.open(os.path.join(D, "clu_cluster.tsv.gz"), "at") as f:
            f.write('GA|p"4\tGA|p"4\n')
        with gzip.open(os.path.join(D, "chunks", "c0001.claims.jsonl.gz"), "at") as f:
            f.write(json.dumps({"protein_id": 'GA|p"4', "term": "GO:9"}) + "\n")
        scatter.scatter(D, self.path("claims"), os.path.join(D, "coverage.tsv"))
        ga = readjsonl(self.path("claims", "GA.claims.jsonl"))
        self.assertIn(('p"4', "GO:9"), [(r["protein_id"], r["term"]) for r in ga])

    def test_unfinished_chunk(self):
        D = self.build(finish=False)
        self.assertEqual(scatter.scatter(D, self.path("claims"), os.path.join(D, "coverage.tsv")), 3)
        self.assertFalse(os.path.exists(os.path.join(D, "coverage.tsv")))


if __name__ == "__main__":
    unittest.main()
