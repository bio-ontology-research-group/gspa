"""Canonical integration-workflow configuration for the web route.

The portal does not choose integration stages or scoring itself: it writes the same typed
``IntegrationWorkflowConfig`` JSON that ``gspa integrate --workflow-config`` and
``gspa annotate --integration-config`` read (keys are the Java property names, listed in
``gspa-core/src/main/resources/workflow/integration-workflow-config.keys.json``; unknown keys
are rejected by the CLI). Stage order, priors and enforcement are decided by the shared
workflow, never here.
"""
from __future__ import annotations
import json
from pathlib import Path

FBA_SIDECAR = '/app/scripts/stoichiometric_locus_refiner.py'
FBA_PYTHON = '/opt/fba/bin/python'
FBA_MEDIUM = 'minimal'


def workflow_config(*, preset: str, resources: Path, claims: Path, browser: Path,
                    pathways: bool, operons: Path | None = None, fba_gram: str | None = None) -> dict:
    """Settings of one web job.

    ``virus`` has no workflow preset: it integrates sequence evidence without priors or
    enforcement. Every other preset supplies its priors/taxon/enforcement defaults through
    ``--preset`` on the command line (precedence: preset < this file < explicit flags).
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
