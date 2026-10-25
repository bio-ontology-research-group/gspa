# One incomplete archaeal bin

*2026-10-24T14:25:51Z by Showboat 0.6.1*
<!-- showboat-id: d059ec3b-7b26-401c-a249-b9ec6a932ee6 -->

A constructed teaching bin retains three of every five 20 kb chunks from the Pyrococcus reference genome. This approximately 60% retained sequence is not a real environmental MAG or a measured completeness score. The mag preset uses meta gene calling and avoids completeness/coherence/gap-filling priors. We explicitly provide the known archaeal domain.

First complete [installation](install.md) and [data and resource setup](resources.md). Run every command from the checkout root; each command starts a new CPU container. Outputs persist in `.tutorial-work/mag/`.

The compact teaching reference contains reviewed proteins from these same organisms, including self/close matches. These runs prove the workflow executes; they do not estimate accuracy. We use one DIAMOND hit per query and disable optional domain/neural predictors. No ontology was downloaded, so ontology-based quality metrics and ancestor relaxation are not demonstrated here; integration uses bundled constraint tables.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:cpu bash -eu <<'GSPA'
gspa integrate --preset mag --print-config > /work/mag-settings.json
python3 - <<'PYCONFIG'
import json
settings = json.load(open('/work/mag-settings.json'))
for key in ['--taxon', '--essential-profile', '--enable-priors', '--enforce-stage2', '--lite']:
    print(key + ': ' + str(settings[key]))
PYCONFIG
GSPA
```

```output
--taxon: None
--essential-profile: None
--enable-priors: consistency,anatomical,genomic_context
--enforce-stage2: False
--lite: True
```

Call genes where needed and search the annotation-bearing DIAMOND database. Logs retain the tool commands and diagnostics. Conversion keeps unique direct DIAMOND GO/EC claims; operon-propagated suggestions in the annotation TSV are excluded, because operon coordinates are supplied separately to integration. This avoids counting the same context twice. Its teaching score is percent identity / 100 clipped to [0.01, 0.99], not a calibrated likelihood. For production, use held-out calibration and a broader reference, and add independent predictors.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:cpu bash -eu <<'GSPA'
mkdir -p /work/mag/gspa_out
gspa annotate --preset mag -k archaea --genome /work/data/mag.fna \
  -c /work/data/cpu.json --db /work/data/reference.dmnd \
  -o /work/mag/annot -t 2 > /work/mag/annotate.log 2>&1 || { tail -60 /work/mag/annotate.log; exit 1; }
python3 /tutorial/results.py claims mag
GSPA
```

```output
2,287 unique DIAMOND GO/EC claims from 384 proteins
Teaching score: DIAMOND identity fraction, clipped to [0.01, 0.99]; not a fitted probability.
```

Integrate the claims and retain provenance. Without a GO ontology, cellular consistency enforcement can remove incompatible terms but cannot relax them to an ancestor. Mixed-community and phage runs deliberately use evidence alone. The results below are counts of hypotheses, not validated functions.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:cpu bash -eu <<'GSPA'
gspa integrate --preset mag --taxon archaea --operons /work/mag/annot/operons_for_integrate.tsv \
  --claims /work/mag/claims.jsonl --out /work/mag/gspa_out/integrated.tsv \
  --provenance /work/mag/gspa_out/provenance.json \
  --enforcement-actions-out /work/mag/gspa_out/actions.tsv \
  > /work/mag/integrate.log 2>&1 || { tail -60 /work/mag/integrate.log; exit 1; }
python3 /tutorial/results.py summary mag
GSPA
```

```output
Proteins: 1,307; with GO/EC search evidence: 384
Integrated hypotheses: GO=1,973; EC=314
Hypotheses at posterior >= 0.5: 1,263
These counts describe this teaching reference; they are not accuracy or completeness estimates.
```

```bash {image}
![Actual annotation counts for mag](.tutorial-work/mag/summary.svg)
```

![Actual annotation counts for mag](63bc6344-2026-10-24.svg)

Build the interactive GSPA browser. The adapter supplies its expected metadata layout from the actual GSPA CDS coordinates (or protein IDs for yeast); it does not run Prokka or invent gene positions. GO/EC IDs remain visible because this compact example has no term-label ontology. The annotations table works locally; the genome tab may need internet access for igv.js.

```bash
docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \
  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \
  gspa-tutorial:cpu bash -eu <<'GSPA'
python3 /tutorial/results.py browser mag
gspa visualize -w /work/mag -g mag \
  --go-obo /work/data/no-labels.obo --ec2go /work/data/no-ec-labels.txt \
  --fasta /work/data/mag.fna \
  > /work/mag/visualize.log 2>&1 || { tail -60 /work/mag/visualize.log; exit 1; }
test -s /work/mag/gspa_out/mag_browser.html
printf 'Interactive browser generated successfully.\n'
GSPA
```

```output
Browser metadata prepared from real CDS coordinates / supplied protein identifiers.
Interactive browser generated successfully.
```

Open `.tutorial-work/mag/gspa_out/mag_browser.html` in your browser. Inspect `integrated.tsv` for final hypotheses, `provenance.json` for supporting evidence, and `actions.tsv` when enforcement acts. Keep the logs and `data/downloads.json` alongside results. For an actual project, replace the input and reference database; retain the biological assumptions appropriate to your unit of analysis.

Re-run this document with `uvx showboat verify docs/tutorials/mag.md` after installation and resource setup.
