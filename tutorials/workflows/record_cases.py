#!/usr/bin/env python3
"""Record each real workflow through Showboat (never hand-edit captured Markdown)."""
import argparse
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
DESCRIPTIONS = {
'bacteria': ('Bacterial genome: E. coli K-12', 'A complete E. coli K-12 MG1655 chromosome, ENA [U00096.3](https://www.ebi.ac.uk/ena/browser/view/U00096.3). Prodigal calls CDS, DIAMOND transfers GO/EC evidence, and the bacterial preset applies cellular genome priors and taxon enforcement.'),
'archaea': ('Archaeal genome: Pyrococcus furiosus', 'The complete Pyrococcus furiosus DSM 3638 genome, ENA [AE009950.1](https://www.ebi.ac.uk/ena/browser/view/AE009950.1). The archaeal preset selects archaeal essentials and constraints; it does not treat this organism as a bacterium.'),
'eukaryote': ('Eukaryotic proteome: budding yeast', 'The reviewed Saccharomyces cerevisiae S288c proteome from [UniProt](https://www.uniprot.org/proteomes/UP000002311). Supply translated proteins or trusted splice-aware gene models for eukaryotes. This run uses the supplied proteins, so it makes no new CDS predictions and has no genomic coordinates.'),
'mag': ('One incomplete archaeal bin', 'A constructed teaching bin retains three of every five 20 kb chunks from the Pyrococcus reference genome. This approximately 60% retained sequence is not a real environmental MAG or a measured completeness score. The mag preset uses meta gene calling and avoids completeness/coherence/gap-filling priors. We explicitly provide the known archaeal domain.'),
'metagenome': ('Mixed community assembly', 'A constructed unbinned community of E. coli, Pyrococcus and lambda sequences split into separate 20 kb contigs. It demonstrates community handling, not a real metagenomic sample. The metagenome preset uses meta gene calling and disables pooled single-organism priors, taxon enforcement and operon inference. Analyze actual bins individually with the MAG workflow.'),
'virus': ('Virus/phage: bacteriophage lambda', 'The complete lambda genome, ENA [J02459.1](https://www.ebi.ac.uk/ena/browser/view/J02459.1). There is currently no virus workflow preset: this tutorial records explicit viral settings. It uses Prodigal meta mode, no cellular completeness assumptions, no operon/pathway inference and evidence-only integration. This gene-calling example is for a bacterial dsDNA phage with genetic code 11; it is not a universal RNA-virus or alternative-code workflow. Supply externally called proteins for those cases.'),
'prokaryote': ('Prokaryote with domain unspecified', 'The same E. coli teaching genome is analyzed without providing its domain, to show the prokaryote preset. It does not silently choose bacterial essentials or assert bacterial taxonomy. In real work, prefer bacteria or archaea when the domain is known.')}


def show(doc, action, text=None):
    cmd = ['uvx', 'showboat', action, str(doc)]
    if action == 'exec': cmd += ['bash']
    if text is not None: cmd += [text]
    subprocess.run(cmd, cwd=ROOT, check=True)


def docker(code):
    return ('docker run --rm -i --user "$(id -u):$(id -g)" --cpus 2 --memory 6g \\\n'
            '  -v "$PWD/.tutorial-work:/work" -v "$PWD/tutorials/workflows:/tutorial:ro" \\\n'
            "  gspa-tutorial:cpu bash -eu <<'GSPA'\n" + code.rstrip() + '\nGSPA')


def record(case):
    title, intro = DESCRIPTIONS[case]
    doc = ROOT / f'docs/tutorials/{case}.md'
    destination = doc
    if doc.exists(): doc = doc.with_stem(doc.stem + '-refreshed')
    show(doc,'init',title)
    show(doc,'note',intro + '\n\nFirst complete [installation](install.md) and [data and resource setup](resources.md). Run every command from the checkout root; each command starts a new CPU container. Outputs persist in `.tutorial-work/' + case + '/`.')
    show(doc,'note','The compact teaching reference contains reviewed proteins from these same organisms, including self/close matches. These runs prove the workflow executes; they do not estimate accuracy. We use one DIAMOND hit per query and disable optional domain/neural predictors. No ontology was downloaded, so ontology-based quality metrics and ancestor relaxation are not demonstrated here; integration uses bundled constraint tables.')
    if case != 'virus':
        show(doc, 'exec', docker(f"""gspa integrate --preset {case} --print-config > /work/{case}-settings.json
python3 - <<'PYCONFIG'
import json
settings = json.load(open('/work/{case}-settings.json'))
for key in ['--taxon', '--essential-profile', '--enable-priors', '--enforce-stage2', '--lite']:
    print(key + ': ' + str(settings[key]))
PYCONFIG"""))
    else:
        show(doc, 'exec', docker('python3 -m json.tool /work/data/virus.json'))
    inp = '--proteins /work/data/eukaryote.faa' if case == 'eukaryote' else f'--genome /work/data/{"bacteria" if case == "prokaryote" else case}.fna'
    preset = '-k virus --mag' if case == 'virus' else '--preset ' + case + (' -k archaea' if case == 'mag' else '')
    config = 'virus' if case == 'virus' else 'cpu'
    show(doc,'note','Call genes where needed and search the annotation-bearing DIAMOND database. Logs retain the tool commands and diagnostics. Conversion keeps unique direct DIAMOND GO/EC claims; operon-propagated suggestions in the annotation TSV are excluded, because operon coordinates are supplied separately to integration. This avoids counting the same context twice. Its teaching score is percent identity / 100 clipped to [0.01, 0.99], not a calibrated likelihood. For production, use held-out calibration and a broader reference, and add independent predictors.')
    code = f'''mkdir -p /work/{case}/gspa_out
gspa annotate {preset} {inp} \\
  -c /work/data/{config}.json --db /work/data/reference.dmnd \\
  -o /work/{case}/annot -t 2 > /work/{case}/annotate.log 2>&1 || {{ tail -60 /work/{case}/annotate.log; exit 1; }}
python3 /tutorial/results.py claims {case}'''
    show(doc,'exec',docker(code))
    opts = '--lite --enable-priors "" --no-enforce-stage2' if case == 'virus' else '--preset ' + case + (' --taxon archaea' if case == 'mag' else '')
    if case in ('bacteria','archaea','mag','prokaryote'):
        opts += f' --operons /work/{case}/annot/operons_for_integrate.tsv'
    show(doc,'note','Integrate the claims and retain provenance. Without a GO ontology, cellular consistency enforcement can remove incompatible terms but cannot relax them to an ancestor. Mixed-community and phage runs deliberately use evidence alone. The results below are counts of hypotheses, not validated functions.')
    code = f'''gspa integrate {opts} \\
  --claims /work/{case}/claims.jsonl --out /work/{case}/gspa_out/integrated.tsv \\
  --provenance /work/{case}/gspa_out/provenance.json \\
  --enforcement-actions-out /work/{case}/gspa_out/actions.tsv \\
  > /work/{case}/integrate.log 2>&1 || {{ tail -60 /work/{case}/integrate.log; exit 1; }}
python3 /tutorial/results.py summary {case}'''
    show(doc,'exec',docker(code))
    show(doc,'image',f'![Actual annotation counts for {case}](.tutorial-work/{case}/summary.svg)')
    show(doc,'note','Build the interactive GSPA browser. The adapter supplies its expected metadata layout from the actual GSPA CDS coordinates (or protein IDs for yeast); it does not run Prokka or invent gene positions. GO/EC IDs remain visible because this compact example has no term-label ontology. The annotations table works locally; the genome tab may need internet access for igv.js.')
    code = f'''python3 /tutorial/results.py browser {case}
gspa visualize -w /work/{case} -g {case} \\
  --go-obo /work/data/no-labels.obo --ec2go /work/data/no-ec-labels.txt \\
  {('--fasta /work/data/' + ('bacteria' if case == 'prokaryote' else case) + '.fna') if case != 'eukaryote' else ''} \\
  > /work/{case}/visualize.log 2>&1 || {{ tail -60 /work/{case}/visualize.log; exit 1; }}
test -s /work/{case}/gspa_out/{case}_browser.html
printf 'Interactive browser generated successfully.\\n' '''
    show(doc,'exec',docker(code))
    show(doc,'note',f'Open `.tutorial-work/{case}/gspa_out/{case}_browser.html` in your browser. Inspect `integrated.tsv` for final hypotheses, `provenance.json` for supporting evidence, and `actions.tsv` when enforcement acts. Keep the logs and `data/downloads.json` alongside results. For an actual project, replace the input and reference database; retain the biological assumptions appropriate to your unit of analysis.\n\nRe-run this document with `uvx showboat verify docs/tutorials/{case}.md` after installation and resource setup.')

    if doc != destination: doc.replace(destination)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cases', nargs='*', choices=list(DESCRIPTIONS))
    args = ap.parse_args()
    for case in args.cases or DESCRIPTIONS: record(case)
