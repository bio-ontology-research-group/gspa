"""Write the geNomad / CheckV adapter-certification specs for a Shaheen run (python 3.7).

Usage: make_genomic_specs.py <cert-root> <predictor> <sif> <db> <singularity-bin-dir> <native-sidecar-tsv> <python>
The genuine outcome is derived from the native sidecar TSV of the direct run: data rows
mean SUCCESS, a header-only table means SUCCESS_ZERO (CheckV on a bacterial chromosome).
"""
import json
import sys

root, predictor, sif, db, sing_dir, native, python = sys.argv[1:8]

if predictor == 'gapseq':
    # Nucleotide pathway route: gapseq find -M nucl on the genome. The adapter's documented
    # contract is that nucleotide hits cannot be attributed to proteins, so a genuine run with
    # hits is a typed FAILED with that diagnostic and zero claims; a header-only table is SUCCESS_ZERO.
    rows = [l for l in open(native).read().splitlines()[1:] if l.strip()]
    base = {'similarity': {'enabled': False}, 'domains': {'enabled': False}, 'operons': {'enabled': False},
            'localization': {'enabled': False}, 'antifam': {'enabled': False},
            'pathway': {'enabled': True, 'tool': 'gapseq', 'target': 'genome', 'executablePath': sif, 'threads': 32, 'taxonomy': 'Bacteria'}}
    spec = {'adapter': 'gapseq', 'mode': 'genome',
            'input': {'genome': root + '/bw25113/GCF_000750555.1.fna', 'gff3': root + '/bw25113/GCF_000750555.1.gff'},
            'config': {'input': {'kingdom': 'bacteria'}, 'predictors': base, 'quality': {'enabled': False}},
            'path_dirs': [db], 'tool': 'gapseq', 'stub_kind': 'gapseq',
            'executable_config_key': ['predictors', 'pathway', 'executablePath'],
            'native_output': native, 'oracle_reader': 'none',
            'genuine_outcome': 'FAILED' if rows else 'SUCCESS_ZERO',
            'genuine_detail_contains': ['nucleotide pathway hits require coordinate-to-protein attribution'] if rows else [],
            'genuine_claims': 'none', 'nondeterministic': True, 'malformed_column': 2, 'dropped_protein': 'none',
            # `annotate --threads` overrides pathway.threads; with 1 thread the genome-wide tblastn ran past
            # the harness's 3600 s limit (job 15745883), so the genuine runs get 32 threads and 2 h.
            'cli_threads': 32, 'timeout_seconds': 7200,
            'cases': {'header_only': {'expected': 'SUCCESS_ZERO', 'claims': 'none'},
                      'truncated_row': {'expected': 'FAILED', 'claims': 'none'},
                      'unavailable': {'expected': 'UNAVAILABLE', 'claims': 'none'}},
            'bound': [sif],
            'note': 'Shaheen: gapseq 2.1.0 (conda env /project/k10018/sw/envs/gapseq) find -M nucl on BW25113; the genome route is a documented typed failure when hits exist'}
    out = root + '/results/gapseq/spec.json'
    with open(out, 'w') as f:
        json.dump(spec, f, indent=1)
    print(out, 'genuine_outcome', spec['genuine_outcome'], len(rows), 'native rows')
    sys.exit(0)

G = 'NZ_CP009273.1'
region_type = 'prophage' if predictor == 'genomad' else 'viral_contig'


def tsv(rows):
    head = 'contig_id\tregion_start\tregion_end\tregion_type\tscore\tattributes\n'
    return head + ''.join('%s\t%d\t%d\t%s\t1.0000\tk=v%d\n' % (c, s, e, region_type, i) for i, (c, s, e) in enumerate(rows))


base = {'similarity': {'enabled': False}, 'domains': {'enabled': False}, 'pathway': {'enabled': False},
        'operons': {'enabled': False}, 'localization': {'enabled': False}, 'antifam': {'enabled': False}}
viral = {'genomicSidecar': root + '/run_genomic_predictors.py', 'pythonExecutable': python}
if predictor == 'genomad':
    viral['genomad'] = {'enabled': True, 'dbPath': db, 'sif': sif, 'minScore': 0.7}
else:
    viral['checkv'] = {'enabled': True, 'dbPath': db, 'sif': sif, 'threads': 32, 'minScore': 0.5}
base['viral'] = viral
rows = [l for l in open(native).read().splitlines()[1:] if l.strip()]
cases = {'header_only': {'expected': 'SUCCESS_ZERO', 'claims': 'none'},
         'unavailable': {'expected': 'UNAVAILABLE', 'claims': 'none'},
         'region_over_cds': {'expected': 'SUCCESS', 'output': tsv([(G, 2700, 2900)]), 'expected_proteins': ['WP_001264707.1', 'WP_000241662.1']},
         'region_outside_contig': {'expected': 'PARTIAL', 'output': tsv([(G, 2700, 2900), (G, 9000000, 9000500)]), 'expected_proteins': ['WP_001264707.1', 'WP_000241662.1']},
         'region_unknown_contig': {'expected': 'PARTIAL', 'output': tsv([('contig_not_in_genome', 2700, 2900)]), 'expected_proteins': []}}
if rows:
    cases['malformed'] = {'expected': 'FAILED', 'claims': 'none'}
    cases['truncated_row'] = {'expected': 'FAILED', 'claims': 'none'}
else:  # derive the injected rows from a synthetic region instead
    bad = tsv([(G, 2700, 2900)]).splitlines()[1].split('\t')
    bad[4] = 'not-a-number'
    cases['malformed'] = {'expected': 'FAILED', 'claims': 'none', 'output': tsv([]) + '\t'.join(bad) + '\n'}
    cases['truncated_row'] = {'expected': 'FAILED', 'claims': 'none', 'output': tsv([]) + '\t'.join(tsv([(G, 2700, 2900)]).splitlines()[1].split('\t')[:-2]) + '\n'}
spec = {'adapter': predictor, 'mode': 'genome',
        'input': {'genome': root + '/bw25113/GCF_000750555.1.fna', 'gff3': root + '/bw25113/GCF_000750555.1.gff'},
        'config': {'input': {'kingdom': 'bacteria'}, 'predictors': base, 'quality': {'enabled': False}},
        'path_dirs': [sing_dir, python.rsplit('/', 1)[0]], 'tool': 'run_genomic_predictors.py', 'stub_kind': 'genomic-sidecar',
        'executable_config_key': ['predictors', 'viral', 'genomicSidecar'],
        'native_output': native, 'oracle_reader': 'sidecar_regions',
        'genome_oracle': {'gff3': root + '/bw25113/GCF_000750555.1.gff'},
        'genuine_outcome': 'SUCCESS' if rows else 'SUCCESS_ZERO', 'malformed_column': 4, 'dropped_protein': 'none',
        'cases': cases, 'bound': [sif, root + '/run_genomic_predictors.py', root + '/db/' + ('genomad_db.sha256' if predictor == 'genomad' else 'checkv_db.sha256')],
        'attribution_known_gaps': {k: 'gspa.io.CdsTranslator groups every CDS sharing a protein_id as segments of one protein, so the %d identical copies of this transposase collapse to the first locus; copies inside the region are not attributed (engine finding, not repaired in this round)' % n for k, n in [('WP_000019403.1', 9), ('WP_085947771.1', 5)]},
        'note': 'Shaheen: %s via singularity through the GSPA genomic sidecar; proteins translated from the GFF CDS' % predictor}
out = root + '/results/' + predictor + '/spec.json'
with open(out, 'w') as f:
    json.dump(spec, f, indent=1)
print(out, 'genuine_outcome', spec['genuine_outcome'], len(rows), 'native rows')
