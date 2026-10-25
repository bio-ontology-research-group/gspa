#!/usr/bin/env python3
from record_cases import ROOT, docker, show

doc = ROOT / 'docs/tutorials/resources.md'
destination = doc
if doc.exists(): doc = doc.with_stem(doc.stem + '-refreshed')
show(doc,'init','Download resources and prepare the teaching datasets')
show(doc,'note','Complete [installation](install.md) first. Run commands from the checkout root. This compact CPU tutorial uses about 40 MB of biological input, a small DIAMOND index, and bundled reference tables. It does not download the entire optional resource catalog: GTDB alone exceeds 100 GB compressed. The [download reference](../reference/resources.md) explains production databases, model assets and tools required by individual downloads.')
show(doc,'exec',docker('gspa resources --json > /work/catalog.json\npython3 -c \'import json; d=json.load(open("/work/catalog.json")); print("Catalog entries:", len(d["resources"])); print("Groups:", ", ".join(d["groups"]))\''))
show(doc,'note','Download one model individually, or all bundled models and tables together. The DeepGO-PlusPlus JSONs are integrator coefficients; they are not complete neural checkpoints. The catalog explicitly identifies resources that need licensed or user-supplied assets.')
show(doc,'exec',docker('''gspa download deepgo-plusplus-light --dir /work/resources > /work/download-one.log 2>&1 || { cat /work/download-one.log; exit 1; }
gspa download --group bundled --dir /work/resources > /work/download-bundled.log 2>&1 || { cat /work/download-bundled.log; exit 1; }
python3 -c 'import pathlib; p=pathlib.Path("/work/resources"); import json; ids=json.load(open("/work/catalog.json"))["groups"]["bundled"]; assert all((p/i/".gspa-resource.json").is_file() for i in ids); print("Verified bundled resources:", len(ids)); print("Reusable config:", (p/"gspa-paths.yaml").is_file())' '''))
show(doc,'note','Preview the full optional set before committing storage. The command below really runs the planner, saves its complete output and checks that dry-run created no target directory. For a production machine with enough storage and the required tool executables, remove `--dry-run` from `gspa download --all --dir /work/all-resources`. Manual/licensed entries are reported as unresolved, and failures return a nonzero status.')
show(doc,'exec',docker('''gspa download --all --dry-run --dir /work/all-resources > /work/download-plan.txt
test ! -e /work/all-resources
printf 'Full resource plan saved; dry-run created no download directory.\\n' '''))
show(doc,'note','Fetch versioned nucleotide records from [ENA](https://www.ebi.ac.uk/ena/browser/home) and a small reviewed-protein table from [UniProt](https://rest.uniprot.org/documentation). `data.py` records every URL and SHA256; `data-lock.json` rejects upstream changes instead of silently changing the example. UniProt is a live query: retain `.tutorial-work/data/` for exact replay. A fresh fetch after UniProt changes requires a deliberate lock refresh and regenerated outputs. Reference FASTA headers retain the actual GO/EC identifiers returned by UniProt.')
show(doc,'exec',docker('python3 /tutorial/data.py'))
show(doc,'note','Build a DIAMOND database from those annotation-bearing headers. This reference deliberately includes close/self matches to keep the example small and interpretable. It is not a held-out evaluation. The helper also writes minimal JSON configuration files (valid YAML) that select DIAMOND and turn off optional tools; these files do not override the cellular preset assumptions.')
show(doc,'exec',docker('''diamond makedb --in /work/data/reference.faa --db /work/data/reference --threads 2 > /work/makedb.log 2>&1 || { cat /work/makedb.log; exit 1; }
test -s /work/data/reference.dmnd
python3 -m json.tool /work/data/cpu.json'''))
show(doc,'note','You can now follow [bacteria](bacteria.md), [archaea](archaea.md), [eukaryotes](eukaryote.md), [one incomplete bin](mag.md), [mixed metagenomes](metagenome.md), [phage](virus.md), or [a prokaryote of unspecified domain](prokaryote.md). The optional [CPU/GPU neural tutorial](neural.md) downloads a real ESM-2 model and compares inference on both devices.\n\nVerify this document with `uvx showboat verify docs/tutorials/resources.md`.')

if doc != destination: doc.replace(destination)
