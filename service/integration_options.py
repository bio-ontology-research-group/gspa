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
