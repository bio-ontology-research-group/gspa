"""Joint multi-genome cohort job for the web service: the CLI batch route, run by the worker.

A cohort job is one uploaded ZIP bundle: a manifest ``cohort.tsv`` with the columns of the Nextflow cohort
module (``gspa-nf/modules/cohort.nf``, ``cohortColumns()``), the per-genome claims and layout files it names,
``cohort_metadata.json`` and ``cohort_orthogroups.tsv``. The worker writes the same typed workflow JSON the
Nextflow module composes and runs ``gspa integrate --workflow-config … --batch-manifest … --batch-threads …``:
the two-pass joint cohort inference is the engine's (``IntegrateCommand.runBatchManifest``); nothing of it
is implemented here. The service does not compute orthogroups: like the CLI batch and Nextflow routes it
takes them as supplied cohort evidence.

The kind is offered only when the deployment sets ``GSPA_COHORT_RESOURCES`` (a JSON file with the cohort
resources, same keys as ``cohort.nf --cohort_resources``). Single-genome jobs never touch this module.
"""
from __future__ import annotations
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time
import zipfile

from service import jobs as j

# The manifest columns of the Nextflow cohort module (cohortColumns); `out_path` is owned by the service.
COLUMNS = ['genome_id', 'claims_path', 'out_path', 'taxon', 'essential_profile', 'operons_path',
           'gaps_path', 'max_iter', 'enforce_stage2', 'enforce_mode', 'enable_priors',
           'genome_layout_path', 'orthogroups_path', 'cluster_consensus_path', 'fba_medium', 'fba_gram']
REQUIRED_COLUMNS = ('genome_id', 'claims_path', 'genome_layout_path')
# The resource keys the Nextflow cohort module accepts (cohortResources); metadata and orthogroups come
# from the bundle, the rest from the deployment.
RESOURCE_KEYS = ['goOwlFile', 'ec2goFile', 'pathwaysFile', 'thetaFile', 'learnedModelFile',
                 'externalContextFile', 'reactionGraphFile', 'reactionEcAliasesFile', 'diffusionMetsFile',
                 'cohortMetadataFile', 'cohortOrthogroupsFile', 'cohortRoleEvidenceFile',
                 'cohortInteractionEdgesFile', 'cohortReactionDirectionsFile', 'moduleDefsFile', 'moduleFiles',
                 'cohortPathwaysFile',
                 'taxonomyFile', 'taxonConstraintsFile', 'complexTermsFile', 'essentialFunctionsFile', 'overlayFile', 'reasonerCacheDir']
BUNDLE_KEYS = ('cohortMetadataFile', 'cohortOrthogroupsFile')
REQUIRED_RESOURCES = ('reactionGraphFile', 'reactionEcAliasesFile')
# Defaults of the Nextflow cohort module (cohort.nf `settings`).
SETTINGS = {'lite': True, 'essentialProfile': 'bacteria',
            'enabledPriors': 'essentiality,coherence,gap_filling,genomic_context',
            'cohortConfig': {'holdout': 'CLADE', 'contextMode': 'UNIVERSAL'}}
MANIFEST, METADATA, ORTHOGROUPS = 'cohort.tsv', 'cohort_metadata.json', 'cohort_orthogroups.tsv'
ID_RE = re.compile(r'[A-Za-z0-9][A-Za-z0-9._-]{0,63}')
MAX_MEMBERS = 2000
MAX_GENOMES = 64
EMPTY = ('', '-')


class BundleError(ValueError):
    pass


def offered(environ=None) -> bool:
    environ = os.environ if environ is None else environ
    return bool(environ.get('GSPA_COHORT_RESOURCES'))


def resources(environ=None) -> dict:
    """Deployment-wide cohort resources: validated against the Nextflow module's allowed keys."""
    environ = os.environ if environ is None else environ
    path = Path(environ['GSPA_COHORT_RESOURCES'])
    spec = json.loads(path.read_text())
    unknown = set(spec) - set(RESOURCE_KEYS)
    if unknown:
        raise ValueError(f'unknown cohort resources: {sorted(unknown)}')
    if set(spec) & set(BUNDLE_KEYS):
        raise ValueError('cohort metadata and orthogroups come from the uploaded bundle, not from the deployment')
    missing = [k for k in REQUIRED_RESOURCES if k not in spec]
    if missing:
        raise ValueError(f'cohort resources require {missing}')
    out = {}
    for key, value in spec.items():
        values = value if isinstance(value, list) else [value]
        if key != 'moduleFiles' and len(values) != 1:
            raise ValueError(f'{key} requires one file')
        resolved = [str((path.parent / v).resolve()) if not Path(v).is_absolute() else str(Path(v)) for v in values]
        for r in resolved:
            if not Path(r).exists():
                raise ValueError(f'cohort resource {key} does not exist: {r}')
        out[key] = resolved if isinstance(value, list) else resolved[0]
    return out


def _safe_member(name: str) -> str:
    if not name or name.endswith('/'):
        raise BundleError(f'bundle contains a directory entry or empty name: {name!r}')
    p = Path(name)
    if p.is_absolute() or any(part in ('..', '') for part in p.parts) or name.startswith(('/', '\\')) or '\\' in name:
        raise BundleError(f'unsafe path in bundle: {name!r}')
    return name


def extract(zip_path: Path, dest: Path) -> list:
    """Extract the bundle into ``dest`` (regular files only, relative paths only)."""
    with zipfile.ZipFile(zip_path) as archive:
        members = [m for m in archive.infolist() if not m.is_dir()]
        if len(members) > MAX_MEMBERS:
            raise BundleError(f'bundle has more than {MAX_MEMBERS} files')
        names = []
        for m in members:
            name = _safe_member(m.filename)
            if m.file_size > j.MAX_BYTES * 4:
                raise BundleError(f'{name} is larger than the allowed size')
            names.append(name)
        if len(set(names)) != len(names):
            raise BundleError('duplicate file names in bundle')
        dest.mkdir(parents=True, exist_ok=True)
        for m, name in zip(members, names):
            target = dest / name
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(m) as src, target.open('wb') as out:
                shutil.copyfileobj(src, out)
    return names


def _inside(bundle: Path, value: str, what: str) -> Path:
    if value in EMPTY:
        raise BundleError(f'{what} is required')
    if re.search(r'[\t\r\n\x00]', value):
        raise BundleError(f'invalid {what}')
    p = Path(value)
    if p.is_absolute() or any(part == '..' for part in p.parts):
        raise BundleError(f'{what} must be a relative path inside the bundle: {value!r}')
    target = bundle / p
    if not target.is_file():
        raise BundleError(f'{what} not found in bundle: {value!r}')
    return target


def read_bundle(bundle: Path) -> dict:
    """Validated manifest rows (sorted by genome id, paths resolved inside the bundle), metadata, orthogroups."""
    manifest = bundle / MANIFEST
    for name in (MANIFEST, METADATA, ORTHOGROUPS):
        if not (bundle / name).is_file():
            raise BundleError(f'bundle must contain {name}')
    with manifest.open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        header = reader.fieldnames or []
        unknown = set(header) - set(COLUMNS)
        if unknown:
            raise BundleError(f'unknown cohort columns: {sorted(unknown)}')
        for col in REQUIRED_COLUMNS:
            if col not in header:
                raise BundleError(f'cohort.tsv must have the column {col}')
        raw_rows = list(reader)
    if not raw_rows:
        raise BundleError('cohort input is empty')
    if len(raw_rows) > MAX_GENOMES:
        raise BundleError(f'at most {MAX_GENOMES} genomes per cohort job')
    ids = [(r.get('genome_id') or '').strip() for r in raw_rows]
    if any(not ID_RE.fullmatch(i) for i in ids) or len(set(ids)) != len(ids):
        raise BundleError('cohort needs distinct, safe genome IDs')
    rows = []
    for raw in sorted(raw_rows, key=lambda r: r['genome_id'].strip()):
        if None in raw or any(v is None for v in raw.values()):
            raise BundleError(f'{raw.get("genome_id")}: malformed manifest row')
        row = {'genome_id': raw['genome_id'].strip()}
        for key in COLUMNS:
            if key in ('genome_id', 'out_path'):
                continue
            value = (raw.get(key) or '').strip()
            if key in ('claims_path', 'genome_layout_path'):
                row[key] = _inside(bundle, value, f'{row["genome_id"]}: {key}')
            elif value not in EMPTY:
                if re.search(r'[\t\r\n\x00]', value):
                    raise BundleError(f'{row["genome_id"]}: invalid {key}')
                row[key] = _inside(bundle, value, f'{row["genome_id"]}: {key}') if key.endswith('_path') else value
        rows.append(row)
    try:
        metadata = json.loads((bundle / METADATA).read_text())
        genomes = set(metadata['genomes'])
    except (ValueError, KeyError, TypeError):
        raise BundleError(f'{METADATA} is not a cohort metadata document')
    if genomes != set(ids):
        raise BundleError('cohort_metadata.json must describe exactly the genomes of cohort.tsv')
    with (bundle / ORTHOGROUPS).open() as fh:
        first = fh.readline().rstrip('\n').split('\t')
    if len(first) < 3 or first[0].strip() != 'genome_id' or first[1].strip() != 'protein_id' or first[2].strip().lower() not in ('orthogroup_id', 'orthogroup'):
        raise BundleError(f'{ORTHOGROUPS} must start with the header genome_id, protein_id, orthogroup_id')
    return {'rows': rows, 'genome_ids': [r['genome_id'] for r in rows], 'metadata': bundle / METADATA, 'orthogroups': bundle / ORTHOGROUPS}


def count_proteins(rows) -> int:
    n = 0
    for row in rows:
        with row['genome_layout_path'].open() as fh:
            n += max(0, sum(1 for _ in fh) - 1)
    return n


def accept(folder: Path) -> dict:
    """Submission: extract and validate the uploaded bundle; raises BundleError (→ 400) on any defect."""
    try:
        extract(folder / 'cohort.zip', folder / 'bundle')
    except zipfile.BadZipFile:
        raise BundleError('Upload a ZIP bundle with cohort.tsv, the claims and layout files, cohort_metadata.json and cohort_orthogroups.tsv.')
    bundle = read_bundle(folder / 'bundle')
    return {'genomes': bundle['genome_ids'], 'proteins': count_proteins(bundle['rows'])}


def workflow_config(*, resources: dict, metadata: Path, orthogroups: Path, out_dir: Path,
                    feature_export: bool = False) -> dict:
    """The typed workflow JSON of one cohort job: the Nextflow cohort module's settings and keys."""
    config = dict(SETTINGS)
    config['cohortConfig'] = dict(SETTINGS['cohortConfig'])
    for key in sorted(resources):
        if key == 'overlayFile':
            continue      # applied by the CLI through --workflow-overlay, never merged here
        config[key] = list(resources[key]) if isinstance(resources[key], list) else resources[key]
    config['cohortMetadataFile'] = str(metadata)
    config['cohortOrthogroupsFile'] = str(orthogroups)
    config['cohortInference'] = True
    config['provenanceFile'] = str(out_dir / 'provenance.json')
    config['suggestionsOut'] = str(out_dir / 'suggestions.tsv')
    if feature_export:
        config['featureExportDir'] = str(out_dir / 'features')
    return config


def batch_manifest(rows, out_dir: Path) -> str:
    """The CLI batch manifest: the module's columns, per-genome output under ``out_dir/<genome>/integrated.tsv``."""
    lines = ['\t'.join(COLUMNS)]
    for row in rows:
        values = []
        for key in COLUMNS:
            if key == 'out_path':
                values.append(str(out_dir / row['genome_id'] / 'integrated.tsv'))
            else:
                values.append(str(row.get(key, '-')))
        lines.append('\t'.join(values))
    return '\n'.join(lines) + '\n'


def integrate_command(binary, config_file: Path, manifest: Path, threads, overlay: Path | None = None) -> list[str]:
    command = [str(binary), 'integrate', '--workflow-config', str(config_file), '--batch-manifest', str(manifest),
               '--batch-threads', str(threads)]
    if overlay is not None:
        command += ['--workflow-overlay', str(overlay)]
    return command


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


# The batch route names a genome's companion files after its table (``integrated.tsv.provenance.json`` …); the
# renderer reads the single-job layout (``gspa_out/provenance.json``). Only the names differ.
BATCH_TO_RUN_NAMES = {'integrated.tsv.provenance.json': 'provenance.json', 'integrated.tsv.enforcement.tsv': 'actions.tsv'}


def render_run_dir(run: Path, target: Path) -> Path:
    """A run directory for the renderer: every batch artifact linked under the name the single-job layout uses."""
    if target.is_symlink():
        target.unlink()
    target.mkdir(parents=True, exist_ok=True)
    for old in target.iterdir():
        if old.is_symlink():
            old.unlink()
    for src in sorted(run.iterdir()):
        if not src.is_file():
            continue
        link = target / BATCH_TO_RUN_NAMES.get(src.name, src.name)
        link.symlink_to(os.path.relpath(src, target))
    return target


def render_genome(folder: Path, genome: str, layout: Path, go_obo: Path, ec2go: str, python: str | None) -> Path:
    """Results browser of one genome through the packaged `gspa visualize` (same renderer as single jobs)."""
    run = folder / 'gspa_out' / 'cohort' / genome
    work = folder / 'render' / genome
    adapter = work / 'prokka_out'
    adapter.mkdir(parents=True, exist_ok=True)
    render_run_dir(run, work / 'gspa_out')
    with layout.open(newline='') as fh:
        coords = {r['protein_id']: r for r in csv.DictReader(fh, delimiter='\t')}
    with (adapter / f'{genome}.tsv').open('w') as tsv, (adapter / f'{genome}.gff').open('w') as gff:
        writer = csv.writer(tsv, delimiter='\t')
        writer.writerow(['locus_tag', 'ftype', 'length_bp', 'product'])
        gff.write('##gff-version 3\n')
        for pid, r in coords.items():
            writer.writerow([pid, 'CDS', int(r['end']) - int(r['start']) + 1, ''])
            gff.write('\t'.join([r['contig'], 'GSPA', 'CDS', r['start'], r['end'], '.', r['strand'], '0', f'ID={pid};locus_tag={pid}']) + '\n')
    out = run / 'browser.html'
    argv = [j.BIN, 'visualize', '--workdir', str(work), '--run-dir', 'gspa_out', '--genome-id', genome,
            '--out', str(out), '--go-obo', str(go_obo), '--ec2go', str(ec2go)]
    if python:
        argv += ['--python', python]
    j.command(folder, f'Building the results browser of {genome}', argv)
    if not out.is_file():
        raise RuntimeError(f'The results browser of {genome} was not generated.')
    return out


def run_cohort(folder: Path):
    info = j.state(folder, status='running', started=time.time(), stage='Preparing the joint cohort run')
    res = resources()
    bundle = read_bundle(folder / 'bundle')
    out_dir = folder / 'gspa_out' / 'cohort'
    out_dir.mkdir(parents=True, exist_ok=True)
    feature_export = str(os.environ.get('GSPA_FEATURE_EXPORT', '')).strip().lower() in ('1', 'true', 'yes', 'on')
    config = workflow_config(resources=res, metadata=bundle['metadata'], orthogroups=bundle['orthogroups'],
                             out_dir=out_dir, feature_export=feature_export)
    config_file = folder / 'cohort.workflow.json'
    config_file.write_text(json.dumps(config, indent=2))
    manifest = folder / 'batch.tsv'
    manifest.write_text(batch_manifest(bundle['rows'], out_dir))
    overlay = Path(res['overlayFile']) if 'overlayFile' in res else None
    j.command(folder, f'Integrating {len(bundle["rows"])} genomes jointly (cohort context, two passes)',
              integrate_command(j.BIN, config_file, manifest, j.THREADS, overlay))
    genomes = {}
    for row in bundle['rows']:
        table = out_dir / row['genome_id'] / 'integrated.tsv'
        if not table.is_file():
            raise RuntimeError(f'No prediction table was written for {row["genome_id"]}; inspect the run log.')
        with table.open() as fh:
            genomes[row['genome_id']] = {'annotations': max(0, sum(1 for _ in fh) - 1)}
    go_obo = Path(os.environ.get('GSPA_GO_OBO', str(j.ASSETS / 'go.obo')))
    ec2go = res.get('ec2goFile') or str(j.RESOURCES / 'ec2go/ec2go')
    for row in bundle['rows']:
        render_genome(folder, row['genome_id'], row['genome_layout_path'], go_obo, ec2go, os.environ.get('GSPA_VIZ_PYTHON'))
    receipts = {k: ([_sha(Path(p)) for p in v] if isinstance(v, list) else (_sha(Path(v)) if Path(v).is_file() else None)) for k, v in res.items()}
    bundle_files = {str(p.relative_to(folder / 'bundle')): _sha(p) for p in sorted((folder / 'bundle').rglob('*')) if p.is_file()}
    (folder / 'manifest.json').write_text(json.dumps({
        'profile': 'cohort-v1', 'kind': 'cohort', 'genomes': bundle['genome_ids'], 'settings': {k: config[k] for k in SETTINGS},
        'cohort_resources': res, 'cohort_resources_sha256': receipts, 'bundle_sha256': bundle_files,
        'engine': j.BIN, 'batch_threads': j.THREADS, 'feature_export': feature_export, 'created': info['created'],
        'note': 'The joint cohort inference is the engine\'s batch route (gspa integrate --batch-manifest); orthogroups and cohort metadata are supplied evidence.'}, indent=2))
    with zipfile.ZipFile(folder / 'results.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in (folder / 'gspa_out').rglob('*'):
            if path.is_file():
                archive.write(path, path.relative_to(folder))
        for name in ('bundle/' + MANIFEST, 'bundle/' + METADATA, 'cohort.workflow.json', 'batch.tsv', 'manifest.json', 'run.log'):
            if (folder / name).is_file():
                archive.write(folder / name, name)
    j.state(folder, status='complete', stage='Results ready', finished=time.time(), genomes=bundle['genome_ids'],
            genome_results=genomes, annotations=sum(g['annotations'] for g in genomes.values()))
