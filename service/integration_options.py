"""Canonical integration-workflow configuration for the web route.

The portal does not choose integration stages or scoring itself: it writes the same typed
``IntegrationWorkflowConfig`` JSON that ``gspa integrate --workflow-config`` and
``gspa annotate --integration-config`` read (keys are the Java property names, listed in
``gspa-core/src/main/resources/workflow/integration-workflow-config.keys.json``; unknown keys
are rejected by the CLI). Stage order, priors and enforcement are decided by the shared
workflow, never here.

Optional stages (learned model, feature export, metabolic experiments) are switched on per
deployment through environment variables read by :func:`deployment_options`; with none set, the
configuration written for a job is identical to the one written before these options existed.
"""
from __future__ import annotations
import json
import os
from pathlib import Path

# Image layout. ``GSPA_APP_ROOT`` lets the same worker run outside the container (tests, route
# comparisons) without changing the paths the container uses.
APP_ROOT = os.environ.get('GSPA_APP_ROOT', '/app')
FBA_SIDECAR = f'{APP_ROOT}/scripts/stoichiometric_locus_refiner.py'
FBA_PYTHON = os.environ.get('GSPA_FBA_PYTHON', '/opt/fba/bin/python')
FBA_MEDIUM = 'minimal'
EXPERIMENT_SIDECAR = f'{APP_ROOT}/scripts/metabolic_experiments.py'
EXPERIMENT_MODES = ('manifest', 'evidence_guided')
EXPERIMENT_STRATEGIES = ('fixed', 'random', 'evidence')   # qwen (language-model choice) is not offered on the web route


def _truthy(value) -> bool:
    return str(value).strip().lower() in ('1', 'true', 'yes', 'on') if value is not None else False


def deployment_options(environ=None) -> dict:
    """Optional stages a deployment switches on for every web job.

    ``GSPA_LEARNED_MODEL``: path of a frozen learned model (``learnedModelFile``); the shared workflow
    validates its contract and abstains explicitly when the inputs are unsupported.
    ``GSPA_FEATURE_EXPORT``: truthy to export the learned-stage feature table next to the browser.
    ``GSPA_EXPERIMENT_MANIFEST``: path of a frozen metabolic-experiment manifest; optional
    ``GSPA_EXPERIMENT_MODE`` (manifest | evidence_guided), ``GSPA_EXPERIMENT_STRATEGY``
    (fixed | random | evidence), ``GSPA_EXPERIMENT_SIDECAR`` and ``GSPA_EXPERIMENT_PYTHON``.
    Unset variables add nothing. Invalid values raise so a misconfigured deployment fails at
    start-up rather than silently running jobs without the requested stage.
    """
    environ = os.environ if environ is None else environ
    options: dict = {}
    if environ.get('GSPA_LEARNED_MODEL'):
        options['learned_model'] = Path(environ['GSPA_LEARNED_MODEL'])
    if _truthy(environ.get('GSPA_FEATURE_EXPORT')):
        options['feature_export'] = True
    if environ.get('GSPA_EXPERIMENT_MANIFEST'):
        mode = environ.get('GSPA_EXPERIMENT_MODE', 'manifest')
        strategy = environ.get('GSPA_EXPERIMENT_STRATEGY', 'evidence')
        if mode not in EXPERIMENT_MODES or strategy not in EXPERIMENT_STRATEGIES:
            raise ValueError(f'GSPA_EXPERIMENT_MODE must be one of {EXPERIMENT_MODES} and '
                             f'GSPA_EXPERIMENT_STRATEGY one of {EXPERIMENT_STRATEGIES}')
        options['experiments'] = {'manifest': Path(environ['GSPA_EXPERIMENT_MANIFEST']), 'mode': mode, 'strategy': strategy,
                                  'sidecar': environ.get('GSPA_EXPERIMENT_SIDECAR', EXPERIMENT_SIDECAR),
                                  'python': environ.get('GSPA_EXPERIMENT_PYTHON', FBA_PYTHON)}
    return options


def workflow_config(*, preset: str, resources: Path, claims: Path, browser: Path,
                    pathways: bool, operons: Path | None = None, fba_gram: str | None = None,
                    learned_model: Path | None = None, feature_export: bool = False,
                    experiments: dict | None = None) -> dict:
    """Settings of one web job.

    ``virus`` has no workflow preset: it integrates sequence evidence without priors or
    enforcement. Every other preset supplies its priors/taxon/enforcement defaults through
    ``--preset`` on the command line (precedence: preset < this file < explicit flags).
    The three optional stages are off unless the deployment enables them (``deployment_options``).
    """
    config: dict = {
        'claimsFile': str(claims),
        'outFile': str(browser / 'integrated.tsv'),
        'provenanceFile': str(browser / 'provenance.json'),
        'enforcementActionsOut': str(browser / 'actions.tsv'),
        'goOwlFile': str(resources / 'go/go.owl'),
        'ec2goFile': str(resources / 'ec2go/ec2go'),
    }
    if preset == 'virus':
        config.update(lite=True, enabledPriors='', enforceStage2=False)
    if pathways:
        config['pathwaysFile'] = str(resources / 'kegg_pathways.tsv')
        config['moduleFiles'] = [str(resources / 'kegg_modules.tsv')]
        # KEGG module definitions with step logic: evidence-gated completeness, module refinement
        # and the per-module report the results page shows.
        if (resources / 'kegg_module_defs.json').exists():
            config['moduleDefsFile'] = str(resources / 'kegg_module_defs.json')
    if operons is not None:
        config['operonsFile'] = str(operons)
    if fba_gram is not None:
        config.update(enableFbaRefiner=True, fbaDataDir=str(resources / 'fba_data'),
                      fbaSidecar=FBA_SIDECAR, fbaPythonBin=FBA_PYTHON, fbaGram=fba_gram, fbaMedium=FBA_MEDIUM)
    if learned_model is not None:
        config['learnedModelFile'] = str(learned_model)
    if feature_export:
        config['featureExportDir'] = str(browser / 'features')
    if experiments:
        config.update(experimentManifestFile=str(experiments['manifest']), experimentMode=experiments['mode'],
                      experimentStrategy=experiments['strategy'], experimentSidecar=str(experiments['sidecar']),
                      experimentPythonBin=str(experiments['python']))
    return config


INTEGRATED_TYPES = ('GO', 'EC')


def build_claims(raw: Path, out: Path, *, exclude_sources=(), extra_claims=(), aspects=None) -> set:
    """Claims file for ``gspa integrate`` from the annotator's raw claims.

    ``raw`` is the JSONL written by ``gspa annotate --claims-out`` (the typed EvidenceClaim contract:
    protein id, function, source, evidence type, raw and calibrated scores, metadata). Lines are
    carried verbatim; the only selection is the portal's scope (GO and EC claims, not operon context,
    minus ``exclude_sources``) and the only edit is filling a missing GO aspect from the reference.
    ``extra_claims`` (PSORTb rows: ``type``, ``value``, ``score``, ``aspect``) are appended with a raw
    score and no probability, so the integrator calibrates them with the same table as the rest.
    Returns the set of (protein, type, function, source) keys written.
    """
    aspects = aspects or {}
    seen = set()
    with out.open('w') as fh:
        for line in raw.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            key = (rec['protein_id'], rec['function_type'], rec['function_id'], rec['source'])
            if rec['function_type'] not in INTEGRATED_TYPES or rec['source'] == 'operon' \
                    or rec['source'] in exclude_sources or key in seen:
                continue
            if rec['function_type'] == 'GO' and not rec.get('go_aspect') and aspects.get(rec['function_id']):
                rec['go_aspect'] = aspects[rec['function_id']]
                line = json.dumps(rec)
            seen.add(key)
            fh.write(line + '\n')
        for row in extra_claims:
            key = (row['protein_id'], row['type'], row['value'], row['source'])
            if key in seen:
                continue
            seen.add(key)
            fh.write(json.dumps(dict(protein_id=row['protein_id'], function_type=row['type'],
                                     function_id=row['value'], source=row['source'], raw_score=float(row['score']),
                                     go_aspect=aspects.get(row['value'], row.get('aspect', '')))) + '\n')
    return seen


def integrate_command(binary, preset: str, config_file: Path) -> list[str]:
    """The ``gspa integrate`` argv for a job: all settings travel in ``config_file``."""
    command = [str(binary), 'integrate', '--workflow-config', str(config_file)]
    if preset != 'virus':
        command += ['--preset', preset]
    return command


def write_workflow_config(folder: Path, config: dict) -> Path:
    path = folder / 'integration.workflow.json'
    path.write_text(json.dumps(config, indent=2))
    return path
