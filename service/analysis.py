"""Full, opt-out analysis profile shared by all portal submissions."""
from __future__ import annotations
import csv
import hashlib
from functools import lru_cache
import subprocess
import json
import os
from pathlib import Path
import re
import shutil
import time
import urllib.request
import zipfile
from urllib.parse import unquote
from service import jobs as j


def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


@lru_cache(maxsize=1)
def resource_provenance():
    root = j.RESOURCES
    paths = [root/p for p in ('diamond/swissprot.dmnd','go/go.owl','ec2go/ec2go','kegg_pathways.tsv','kegg_modules.tsv','reasoner-cache/has_part_pairs.tsv')]
    paths += list((root/'fba_data').glob('*'))
    hashes = {}
    for path in paths:
        if path.is_file():
            with path.open('rb') as f: hashes[str(path.relative_to(root))] = hashlib.file_digest(f,'sha256').hexdigest()
    return {'sha256': hashes,
        'versions': {tool:subprocess.check_output([tool,'--version'],stderr=subprocess.STDOUT,text=True).strip() for tool in ('diamond','amrfinder','antismash')},
        'amrfinder_database': (root/'tool-data/amrfinder/latest').resolve().name,
        'antismash_database_directories': sorted(str(p.relative_to(root/'tool-data/antismash')) for p in (root/'tool-data/antismash').glob('*/*') if p.is_dir()),
        'fba_numerics': {'numpy':'2.2.6','scipy':'1.15.3','libsbml':'5.21.1'},
        'psortb_image': 'brinkmanlab/psortb_commandline:1.0.2',
        'pathway_sources': ['https://rest.kegg.jp/link/pathway/ec','https://rest.kegg.jp/link/enzyme/module'],
        'metabolic_reference': 'GSPA PSB revision reference_data: ModelSEED v6RE universal matrix, tight Gram-specific bounds and medium.pkl'}


def run_full(folder):
    info = j.state(folder, status='running', started=time.time(), stage='Preparing full analysis')
    kind, preset = info['kind'], info['preset']
    enabled = {key: info.get('analyses', {}).get(key, True) for key in j.DEFAULT_ANALYSES}
    supplied = kind == 'proteins'
    prok = preset in ('bacteria', 'archaea', 'prokaryote', 'mag', 'metagenome')
    out, browser = folder / 'annot', folder / 'gspa_out'
    browser.mkdir(exist_ok=True)
    input_path = folder / ('sample.faa' if supplied else 'sample.fna')
    resources = j.RESOURCES
    records = {key: {'status': 'pending' if on else 'skipped', 'detail': '' if on else 'Disabled at submission'} for key, on in enabled.items()}
    def record(key, status, detail, **extra):
        records[key] = dict(status=status, detail=detail, **extra)
        (browser / 'analyses.json').write_text(json.dumps(records, indent=2))
    def applicable(key, condition, why):
        if not enabled[key]: return False
        if not condition: record(key, 'not_applicable', why)
        return condition
    def done(key, detail, **extra): record(key, 'complete', detail, **extra)
    (browser / 'analyses.json').write_text(json.dumps(records, indent=2))
    want_operons = applicable('operons', prok and not supplied, 'Requires prokaryotic CDS coordinates.')
    config = {'predictors': {'domains': {'enabled': False},
              'similarity': {'enabled': enabled['homology'] or enabled['names'],
                 'database': str(resources / 'diamond/swissprot.dmnd'), 'assignProductNames': enabled['names']},
              'operons': {'enabled': want_operons},
              'neural': {'deepGoPlusPlusLight': {'minScore': 0.5}}},
              'integration': {'enabled': False}, 'quality': {'enabled': False}}
    (folder / 'config.json').write_text(json.dumps(config, indent=2))
    ann = [j.BIN, 'annotate', '--proteins' if supplied else '--genome', str(input_path),
           '-o', str(out), '-t', j.THREADS, '-c', str(folder / 'config.json'),
           '--base-predictor', 'light', '--neural-sidecar', '/app/benchmark/neural/run_neural_predictors.py',
           '--deepgo-plusplus-light-assets', str(j.ASSETS), '--min-score', '0.5']
    if kind == 'gff3': ann += ['--gff3', str(folder / 'sample.gff3')]
    ann += ['-k', 'virus', '--mag'] if preset == 'virus' else ['--preset', preset]
    j.command(folder, 'Predicting GO, enzyme functions and homology names', ann)
    j.verify_prediction(folder / 'run.log')
    if enabled['homology'] or enabled['names']:
        log = (folder / 'run.log').read_text()
        if 'diamond failed:' in log or not re.search(r'diamond: \d+ annotations assigned', log):
            raise RuntimeError('Swiss-Prot homology analysis did not complete; inspect the run log.')
    proteins = table(next(out.glob('*_proteins.tsv')))
    annotation_path = next(out.glob('*_annotations.tsv'))
    rows = table(annotation_path)
    if not enabled['homology']:
        # DIAMOND may still run solely to assign names. Its functional claims
        # must also disappear from the browser's raw/tool-only annotation layer.
        with annotation_path.open() as fh:
            fields = next(csv.reader(fh, delimiter='\t'))
        rows = [r for r in rows if r['source'] != 'diamond']
        with annotation_path.open('w') as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
            writer.writeheader(); writer.writerows(rows)
    faa = next(out.glob('*.faa'))
    gff = next(out.glob('*.gff3'))
    if enabled['homology']: done('homology', 'DIAMOND against reviewed UniProtKB/Swiss-Prot; GO and EC transfer.', ec_annotations=sum(r['type']=='EC' for r in rows))
    if enabled['names']: done('names', 'Input products retained; otherwise a putative Swiss-Prot product requires ≥40% identity and ≥80% query and subject coverage.', assigned=sum('homology' in p['product_source'] for p in proteins))
    if want_operons:
        for name in ('operons.tsv',):
            if (out / name).exists(): shutil.copyfile(out / name, browser / name)
        done('operons', 'Candidate operons: consecutive CDS on the same strand with an intergenic gap ≤300 bp. Functional labels are inferred, not established operon names.')

    if enabled['amr']:
        dest = folder / 'amrfinder_out'; dest.mkdir(exist_ok=True)
        amrfile = dest / 'sample.amr.tsv'
        j.command(folder, 'Screening AMR, stress and virulence proteins with AMRFinderPlus',
            ['amrfinder', '-p', str(faa), '--plus', '-d', str(resources / 'tool-data/amrfinder/latest'),
             '--threads', j.THREADS, '-o', str(amrfile)])
        text = amrfile.read_text().replace('Protein identifier\t', 'Protein id\t', 1)
        amrfile.write_text(text)
        hits = table(amrfile)
        done('amr', 'AMRFinderPlus protein search; acquired/intrinsic determinants and Plus hits. Organism-specific resistance mutations were not tested. A hit alone does not establish phenotypic resistance.', hits=len(hits), amr_hits=sum(r.get('Type')=='AMR' for r in hits), stress_hits=sum(r.get('Type')=='STRESS' for r in hits), virulence_hits=sum(r.get('Type')=='VIRULENCE' for r in hits))

    if applicable('bgc', prok and not supplied, 'BGC detection requires genomic CDS coordinates and a supported bacterial/archaeal context.'):
        dest = folder / 'antismash_out'
        # GFF importers can merge alternative CDS sharing one Parent gene. Supply
        # one GenBank CDS per protein, with the existing translated sequence.
        from Bio import SeqIO
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        seqs = {pid: seq for pid, seq in j.fasta(faa)}
        segments = {}
        for line in gff.read_text().split('##FASTA')[0].splitlines():
            fields = line.split('\t')
            if len(fields) != 9 or fields[2] != 'CDS': continue
            attrs = {k:unquote(v) for k,v in (a.split('=',1) for a in fields[8].split(';') if '=' in a)}
            pid = attrs.get('protein_id',attrs.get('ID',attrs.get('locus_tag')))
            if pid in seqs: segments.setdefault(pid,[]).append((int(fields[3])-1,int(fields[4])))
        records_by_id = {rec.id: rec for rec in SeqIO.parse(str(input_path),'fasta')}
        for rec in records_by_id.values(): rec.annotations['molecule_type']='DNA'
        for protein in proteins:
            pid=protein['protein_id']; spans=sorted(set(segments[pid]))
            strand=1 if protein['strand']=='+' else -1
            qualifiers={'locus_tag':[pid],'protein_id':[pid],'translation':[seqs[pid].rstrip('*')],
                        'product':[protein['product'] or 'uncharacterized protein'],
                        'original_locus_tag':[protein['locus_tag']]}
            if protein['gene']: qualifiers['gene']=[protein['gene']]
            if any(a[1]>b[0] for a,b in zip(spans,spans[1:])):
                # A programmed frameshift can reuse genomic bases. Preserve the
                # exact translated sequence and its outer genomic span explicitly.
                location=FeatureLocation(spans[0][0],max(e for _,e in spans),strand=strand)
                qualifiers['note']=['Overlapping CDS segments; supplied GSPA translation retained. See original GFF3 for segment coordinates.']
            else:
                parts=[FeatureLocation(start,end,strand=strand) for start,end in (spans if strand==1 else reversed(spans))]
                location=parts[0] if len(parts)==1 else CompoundLocation(parts)
            records_by_id[protein['contig']].features.append(SeqFeature(location,type='CDS',qualifiers=qualifiers))
        genbank = folder/'sample.gbk'
        SeqIO.write(list(records_by_id.values()),str(genbank),'genbank')
        j.command(folder, 'Detecting biosynthetic gene clusters with antiSMASH',
            ['antismash', str(genbank), '--genefinding-tool', 'none', '--taxon', 'bacteria', '--minimal', '--cpus', j.THREADS,
             '--databases', str(resources / 'tool-data/antismash'), '--output-dir', str(dest)])
        # The native region file is retained. Normalize coordinates to 1-based inclusive for the browser.
        native = json.loads((dest / 'sample.json').read_text())
        from antismash.common.secmet.locations import location_from_string
        regions = []
        for rec in native['records']:
            entry = {'seq_id': rec['id'], 'regions': []}
            for feature in rec['features']:
                if feature['type'] != 'region': continue
                location = location_from_string(feature['location'])
                entry['regions'].append({'idx': len(entry['regions'])+1, 'start': int(location.start)+1,
                    'end': int(location.end), 'products': feature['qualifiers'].get('product', [])})
            regions.append(entry)
        (dest / 'regions.js').write_text('var recordData = '+json.dumps(regions)+';')
        done('bgc', 'antiSMASH rule-based cluster detection using supplied CDS; minimal mode excludes optional similarity comparisons.', regions=sum(len(r['regions']) for r in regions))

    aspects, hashes = j.reference_metadata()
    psort_claims = []
    if applicable('localization', preset in ('bacteria', 'archaea'), 'PSORTb requires a bacterial Gram type or archaeal context. GO cellular-component predictions remain available.'):
        mode = 'a' if preset == 'archaea' else ('n' if info.get('gram', 'neg') == 'neg' else 'p')
        j.state(folder, stage='Predicting subcellular localization with PSORTb')
        url = os.environ.get('GSPA_PSORTB_URL', 'http://psortb:8000') + '/predict?gram=' + mode
        request = urllib.request.Request(url, data=faa.read_bytes(), headers={'Content-Type':'text/plain'})
        with urllib.request.urlopen(request, timeout=j.TIMEOUT) as response: raw = response.read()
        dest = folder / 'psortb_out'; dest.mkdir(exist_ok=True)
        (dest / 'sample.tsv').write_bytes(raw)
        loc_go = {'Cytoplasmic':'GO:0005737', 'CytoplasmicMembrane':'GO:0005886', 'Periplasmic':'GO:0042597',
                  'OuterMembrane':'GO:0019867', 'Extracellular':'GO:0005576', 'Cellwall':'GO:0005618'}
        loc_rows = table(dest / 'sample.tsv')
        if not loc_rows or 'Localization' not in loc_rows[0]:
            raise RuntimeError('PSORTb returned an unrecognized or empty output; inspect its TSV.')
        for row in loc_rows:
            for location in row['Localization'].replace(' ', '').split('/'):
                term = loc_go.get(location)
                if term:
                    score = min(1.0, float(row['Score']) / 10)
                    psort_claims.append(dict(protein_id=row['SeqID'].split()[0],type='GO',value=term,source='psortb',score=score,aspect='CC'))
        rows.extend(psort_claims)
        done('localization', 'PSORTb 3; scores divided by 10 for integration. Tied best locations retained; unknown predictions remain in the raw results.', assigned=len({r['protein_id'] for r in psort_claims}), gram=mode)

    seen = set()
    with (folder / 'claims.jsonl').open('w') as fh:
        for row in rows:
            if row['type'] not in ('GO', 'EC') or row['source'] == 'operon': continue
            key = row['protein_id'], row['type'], row['value'], row['source']
            if key in seen: continue
            seen.add(key); score = float(row['score'])
            fh.write(json.dumps(dict(protein_id=key[0], function_type=key[1], function_id=key[2], source=key[3],
                raw_score=score, calibrated_prob=min(1-1e-6,max(1e-6,score)),
                go_aspect=aspects.get(key[2], row.get('aspect', ''))))+'\n')
    integrate = [j.BIN, 'integrate', '--claims', str(folder/'claims.jsonl'), '--out', str(browser/'integrated.tsv'),
        '--provenance', str(browser/'provenance.json'), '--enforcement-actions-out', str(browser/'actions.tsv'),
        '--go-owl', str(resources/'go/go.owl'), '--ec2go', str(resources/'ec2go/ec2go')]
    integrate += ['--lite', '--enable-priors', '', '--no-enforce-stage2'] if preset == 'virus' else ['--preset', preset]
    pathway_args = ['--pathways', str(resources/'kegg_pathways.tsv'), '--modules', str(resources/'kegg_modules.tsv')]
    # KEGG module definitions with step logic: evidence-gated module completeness, module-based refinement,
    # and the per-module report (integrated.tsv.modules.json) the results page shows.
    if (resources/'kegg_module_defs.json').exists(): pathway_args += ['--module-defs', str(resources/'kegg_module_defs.json')]
    if enabled['pathways']: integrate += pathway_args
    if want_operons and (out/'operons_for_integrate.tsv').exists(): integrate += ['--operons', str(out/'operons_for_integrate.tsv')]
    want_fba = applicable('metabolism', preset == 'bacteria', 'The installed biomass templates cover bacteria. A validated organism-specific model is required for archaea, eukaryotes, viruses and mixed/unknown domains.')
    if want_fba:
        integrate += ['--enable-fba-refiner', '--fba-data-dir', str(resources/'fba_data'),
            '--fba-sidecar', '/app/scripts/stoichiometric_locus_refiner.py', '--fba-python-bin', '/opt/fba/bin/python', '--fba-gram', info.get('gram','neg'), '--fba-medium', 'minimal']
    if seen:
        j.command(folder, 'Integrating evidence and reconstructing the metabolic model' if want_fba else 'Integrating functional evidence', integrate)
    else:
        (browser/'integrated.tsv').write_text('protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\tn_supporting\tpriors_fired\n')
        (browser/'provenance.json').write_text('[]')
    if want_fba:
        fba = browser/'metabolism/fba_output.json'
        if fba.exists():
            result = json.loads(fba.read_text())
            if not (browser/'metabolism/metabolic_model.xml').exists(): raise RuntimeError('FBA did not export its model.')
            j.command(folder, 'Testing metabolic nutrient dependencies', ['/opt/fba/bin/python', '/app/service/metabolic_conditions.py', '--model', str(browser/'metabolism/metabolic_model.xml'), '--out', str(browser/'metabolism/conditions.json')])
            done('metabolism', 'Stoichiometric reconstruction and flux-certified repair from EC evidence. Biomass flux is conditional on the selected template and supplemented medium, not an observed growth rate.', **result.get('model', {}))
        else: record('metabolism','insufficient_evidence','No enzyme evidence was available to reconstruct a model.')
    integrated = table(browser/'integrated.tsv')
    if seen and not integrated: raise RuntimeError('Integration discarded all predictions; inspect the run log.')

    if enabled['quality'] or enabled['pathways']:
        gaf = folder/'integrated.gaf'
        ec_map = {}
        for line in (resources/'ec2go/ec2go').read_text().splitlines():
            m = re.match(r'(EC:[\d.-]+).*;\s*(GO:\d+)', line)
            if m: ec_map[m[1]] = m[2]
        with gaf.open('w') as fh:
            fh.write('!gaf-version: 2.2\n')
            emitted=set()
            for row in integrated:
                if float(row['posterior_prob']) < 0.5: continue
                term = row['function_id']; term = ec_map.get(term, term)
                if not term.startswith('GO:') or (row['protein_id'], term) in emitted: continue
                emitted.add((row['protein_id'],term))
                aspect = {'BP':'P','MF':'F','CC':'C'}.get(aspects.get(term), 'F')
                fh.write('\t'.join(['GSPA',row['protein_id'],row['protein_id'],'',term,'GSPA:prediction','IEA','',aspect,'','','protein','taxon:1',time.strftime('%Y%m%d'),'GSPA','',''])+'\n')
        evaluate = [j.BIN, 'evaluate', '-i', str(gff), '-a', str(gaf), '--go-owl', str(resources/'go/go.owl'),
            '-k', preset, '-o', str(browser/'quality_gspa.json')]
        if supplied: evaluate += ['--proteins', str(faa)]
        if enabled['pathways']: evaluate += pathway_args + ['--ec2go', str(resources/'ec2go/ec2go')]
        if not enabled['quality']: evaluate += ['--lite']
        else: evaluate += ['--reasoner-cache', str(resources/'reasoner-cache')]
        if preset in ('mag', 'metagenome', 'virus'): evaluate += ['--mag']
        j.command(folder, 'Assessing annotation quality and pathway coverage', evaluate)
        if enabled['quality']: done('quality', 'GSPA annotation completeness, ontology/taxon consistency and functional coherence using integrated scores ≥0.5. These assess annotations, not assembly contamination or experimental accuracy.')
        if enabled['pathways']: done('pathways', 'KEGG module completeness is assessed per step from the module definitions (alternative enzymes, complexes, steps without a GO mapping) and counted only for modules the evidence activates; map coverage counts distinct EC-mapped functions as a breadth measure.' if (resources/'kegg_module_defs.json').exists() else 'KEGG maps and modules projected through EC→GO. Coverage measures distinct mapped functions; alternative KEGG module branches are not a Boolean completeness test.')

    metadata = folder/'prokka_out'; metadata.mkdir(exist_ok=True)
    with (metadata/'sample.tsv').open('w') as fh, (metadata/'sample.gff').open('w') as gf:
        writer=csv.writer(fh,delimiter='\t'); writer.writerow(['locus_tag','ftype','length_bp','gene','EC_number','COG','product','product_source','original_locus_tag','homology_accession','identity','query_coverage','subject_coverage'])
        gf.write('##gff-version 3\n')
        for p in proteins:
            pid=p['protein_id']
            writer.writerow([pid,'CDS',int(p['length_aa'])*3,p['gene'],'','',p['product'],p['product_source'],p['locus_tag'],p['homology_accession'],p['identity'],p['query_coverage'],p['subject_coverage']])
            if not supplied: gf.write('\t'.join([p['contig'],'GSPA','CDS',p['start'],p['end'],'.',p['strand'],'0','ID='+pid+';locus_tag='+pid])+'\n')
    # Manifest records the exact enabled profile and all data receipts; model inputs and raw outputs are archived too.
    receipts = {str(p.relative_to(resources)): json.loads(p.read_text()) for p in resources.glob('*/.gspa-resource.json')}
    (folder/'manifest.json').write_text(json.dumps({'profile':'full-v1','kind':kind,'preset':preset,'input_name':info.get('filename','sample'),'gram':info.get('gram','neg'),
        'analyses':records,'options':enabled,'asset_sha256':hashes,'reference_receipts':receipts,'resources':resource_provenance(),
        'model_sha256':hashlib.sha256(Path('/app/deepgo-plusplus/models/deepgo_plusplus_light_fast.json').read_bytes()).hexdigest(),
        'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.glob('sample.*')},
        'created':info['created'], 'translation_note':'Existing CDS models are reused for GFF3. Standard code 1/11 translation with genomic-coordinate Sec/Pyl/TERM exceptions; pseudogene CDS excluded.'},indent=2))
    viz=[j.BIN,'visualize','-w',str(folder),'-g','sample','--go-obo',str(j.ASSETS/'go.obo'),'--ec2go',str(resources/'ec2go/ec2go')]
    if not supplied: viz += ['--fasta',str(input_path)]
    j.command(folder,'Building the interactive biological report',viz)
    if not (browser/'sample_browser.html').is_file(): raise RuntimeError('The results browser was not generated.')
    with zipfile.ZipFile(folder/'results.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for sub in ('annot','gspa_out','prokka_out','amrfinder_out','antismash_out','psortb_out'):
            for path in (folder/sub).rglob('*'):
                if path.is_file(): archive.write(path,path.relative_to(folder))
        for name in ('claims.jsonl','config.json','manifest.json','run.log','integrated.gaf','sample.gff3','sample.fna','sample.faa','sample.gbk'):
            if (folder/name).is_file(): archive.write(folder/name,name)
    j.state(folder,status='complete',stage='Results ready',finished=time.time(),proteins=len(proteins),
        annotated_proteins=len({r['protein_id'] for r in integrated}),annotations=len(integrated),
        sources=sorted({r['source'] for r in rows if r['source']!='operon'}),analysis_results=records)
