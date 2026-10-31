#!/usr/bin/env python3
"""Bind supplied model topology and genotype mappings without reading outcomes."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import metabolic_experiments as ME
import prepare_metabolic_benchmarks as P
from prepare_nutrient_environments import sha


def register(inputs, models):
    records = [json.loads(line) for line in inputs.read_text().splitlines()]
    result = {'schema': 'gspa.nutrient_model_preflight/1', 'inputs_sha256': sha(inputs),
              'purpose': 'Known-genotype phenotype evaluation; model GPRs are not independent function truth.',
              'implementation_sha256': sha(Path(__file__)), 'models': {}}
    for species, path in sorted(models.items()):
        selected = [r for r in records if r['species'] == species]
        genes = sorted({g for r in selected for o in r['organisms'] for g in o['deleted_genes']})
        if not selected:
            raise ValueError('no input cases for model species ' + species)
        model, index = ME.load_sbml_model(path), P.index_sbml(path)
        mappings = []
        for gene in genes:
            locus, method = P.map_gene(gene, index['symbol_index'])
            canonical = model.gene_alias.get(locus)
            mappings.append({'symbol': gene, 'locus': locus, 'mapping': method, 'canonical': canonical,
                             'blocked_reactions': sorted(model.blocked_by_genes({canonical})) if canonical else None,
                             'reactions_with_locus': P.reactions_with_locus(index, locus) if locus else None})
        result['models'][species] = {'path': str(path.resolve()), 'sha256': sha(path),
                                    'metabolites': len(model.met_ids), 'reactions': len(model.rxn_ids),
                                    'gene_products': len(set(model.gene_alias.values())), 'objectives': model.objective_terms,
                                    'maintenance_bound': float(model.lb[model.rxn_index['R_ATPM']]) if 'R_ATPM' in model.rxn_index else None,
                                    'genes': mappings, 'case_count': len(selected)}
    result['unmapped_species'] = sorted({r['species'] for r in records} - set(models))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--ec-model', type=Path, required=True)
    parser.add_argument('--bs-model', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output exists; preserve previous preflight')
    report = register(args.inputs, {'Ec': args.ec_model, 'Bs': args.bs_model})
    args.out.write_text(json.dumps(report, indent=2) + '\n')
