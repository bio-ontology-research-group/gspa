"""Independent C3 check: no reference entry has a panel taxid or a fenced genus (first word or OC lineage)."""
import csv, gzip, json, sys
dat, ref_entries, genomes, ref_report, out = sys.argv[1:6]
rep = json.load(open(ref_report)); genera = set(rep["fenced_genera"]); taxids = set(rep["panel_taxids"])
keep = {r["accession"]: r for r in csv.DictReader(open(ref_entries), delimiter="\t")}
bad_tax, bad_gen, n, acc, oc, os_ = [], [], 0, None, [], ""
with gzip.open(dat, "rt") as fh:
    for line in fh:
        if line.startswith("AC   ") and acc is None: acc = line[5:].split(";")[0].strip()
        elif line.startswith("OS   ") and not os_: os_ = line[5:]
        elif line.startswith("OC   "): oc.append(line[5:].strip())
        elif line.startswith("//"):
            if acc in keep:
                n += 1
                names = {p.strip().rstrip(".").lower() for t in " ".join(oc).split(";") for p in t.split("/") if p.strip()}
                names |= {x.split()[0] for x in names}
                first = os_.replace("Candidatus ", "").replace("[", "").replace("]", "").split()[0].lower()
                if keep[acc]["taxid"] in taxids: bad_tax.append(acc)
                if first in genera or names & genera: bad_gen.append(acc)
            acc, oc, os_ = None, [], ""
res = {"reference_entries_checked": n, "entries_in_table": len(keep), "panel_taxid_entries": len(bad_tax),
       "fenced_genus_entries": len(bad_gen), "examples": (bad_tax + bad_gen)[:20], "pass": not bad_tax and not bad_gen and n == len(keep)}
json.dump(res, open(out, "w"), indent=1); print(res)
