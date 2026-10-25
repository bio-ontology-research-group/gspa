# Bacterial genome: E. coli K-12

*2026-10-24T14:24:43Z by Showboat 0.6.1*
<!-- showboat-id: f43e1cb2-09aa-4b98-8c9a-943ef887b3bc -->

A complete E. coli K-12 MG1655 chromosome, ENA [U00096.3](https://www.ebi.ac.uk/ena/browser/view/U00096.3). Prodigal calls CDS, DIAMOND transfers GO/EC evidence, and the bacterial preset applies cellular genome priors and taxon enforcement.

First complete [installation](install.md) and [data and resource setup](resources.md). Run every command from the checkout root; each command starts a new CPU container. Outputs persist in `.tutorial-work/bacteria/`.

The compact teaching reference contains reviewed proteins from these same organisms, including self/close matches. These runs prove the workflow executes; they do not estimate accuracy. We use one DIAMOND hit per query and disable optional domain/neural predictors. No ontology was downloaded, so ontology-based quality metrics and ancestor relaxation are not demonstrated here; integration uses bundled constraint tables.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:cpu bash -eu <<'GSPA'
gspa integrate --preset bacteria --print-config > /work/bacteria-settings.json
python3 - <<'PYCONFIG'
import json
settings = json.load(open('/work/bacteria-settings.json'))
for key in ['--taxon', '--essential-profile', '--enable-priors', '--enforce-stage2', '--lite']:
    print(key + ': ' + str(settings[key]))
PYCONFIG
GSPA
```

```output
--taxon: bacteria
--essential-profile: bacteria
--enable-priors: essentiality,coherence,consistency,anatomical,gap_filling,genomic_context
--enforce-stage2: True
--lite: True
```

Call genes where needed and search the annotation-bearing DIAMOND database. Logs retain the tool commands and diagnostics. Conversion keeps unique direct DIAMOND GO/EC claims; operon-propagated suggestions in the annotation TSV are excluded, because operon coordinates are supplied separately to integration. This avoids counting the same context twice. Its teaching score is percent identity / 100 clipped to [0.01, 0.99], not a calibrated likelihood. For production, use held-out calibration and a broader reference, and add independent predictors.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:cpu bash -eu <<'GSPA'
mkdir -p /work/bacteria/gspa_out
gspa annotate --preset bacteria --genome /work/data/bacteria.fna \
  -c /work/data/cpu.json --db /work/data/reference.dmnd \
  -o /work/bacteria/annot -t 2 > /work/bacteria/annotate.log 2>&1 || { tail -60 /work/bacteria/annotate.log; exit 1; }
python3 /tutorial/results.py claims bacteria
GSPA
```

```output
23,189 unique DIAMOND GO/EC claims from 3,876 proteins
Teaching score: DIAMOND identity fraction, clipped to [0.01, 0.99]; not a fitted probability.
```

Integrate the claims and retain provenance. Without a GO ontology, cellular consistency enforcement can remove incompatible terms but cannot relax them to an ancestor. Mixed-community and phage runs deliberately use evidence alone. The results below are counts of hypotheses, not validated functions.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:cpu bash -eu <<'GSPA'
gspa integrate --preset bacteria --operons /work/bacteria/annot/operons_for_integrate.tsv \
  --claims /work/bacteria/claims.jsonl --out /work/bacteria/gspa_out/integrated.tsv \
  --provenance /work/bacteria/gspa_out/provenance.json \
  --enforcement-actions-out /work/bacteria/gspa_out/actions.tsv \
  > /work/bacteria/integrate.log 2>&1 || { tail -60 /work/bacteria/integrate.log; exit 1; }
python3 /tutorial/results.py summary bacteria
GSPA
```

```output
Proteins: 4,319; with GO/EC search evidence: 3,876
Integrated hypotheses: GO=21,226; EC=1,939
Hypotheses at posterior >= 0.5: 23,161
These counts describe this teaching reference; they are not accuracy or completeness estimates.
```

```bash {image}
![Actual annotation counts for bacteria](.tutorial-work/bacteria/summary.svg)
```

![Actual annotation counts for bacteria](72b12fd3-2026-10-24.svg)

Build the interactive GSPA browser. The adapter supplies its expected metadata layout from the actual GSPA CDS coordinates (or protein IDs for yeast); it does not run Prokka or invent gene positions. GO/EC IDs remain visible because this compact example has no term-label ontology. The annotations table works locally; the genome tab may need internet access for igv.js.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:cpu bash -eu <<'GSPA'
python3 /tutorial/results.py browser bacteria
gspa visualize -w /work/bacteria -g bacteria \
  --go-obo /work/data/no-labels.obo --ec2go /work/data/no-ec-labels.txt \
  --fasta /work/data/bacteria.fna \
  > /work/bacteria/visualize.log 2>&1 || { tail -60 /work/bacteria/visualize.log; exit 1; }
test -s /work/bacteria/gspa_out/bacteria_browser.html
printf 'Interactive browser generated successfully.\n'
GSPA
```

```output
Browser metadata prepared from real CDS coordinates / supplied protein identifiers.
Interactive browser generated successfully.
```

Open `.tutorial-work/bacteria/gspa_out/bacteria_browser.html` in your browser. Inspect `integrated.tsv` for final hypotheses, `provenance.json` for supporting evidence, and `actions.tsv` when enforcement acts. Keep the logs and `data/downloads.json` alongside results. For an actual project, replace the input and reference database; retain the biological assumptions appropriate to your unit of analysis.

Re-run this document with `uvx showboat verify docs/tutorials/bacteria.md` after installation and resource setup.
