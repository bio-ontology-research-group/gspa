import gzip
import hashlib
import json
import pickle
import shutil
import tempfile
import tracemalloc
import unittest
from pathlib import Path

import numpy as np

from benchmark.v3 import prepare_price as P
from benchmark.v3 import prepare_price_context as C
from benchmark.v3.prepare_price_layout import digest

ASM = 'GCA_T.1'
COLUMNS = ['EC:1.1.1.1', '1.1.1.2', '1.1.1.n11', '1.1', '2.7.1.1', '3.1.1.1']
FASTA = ['P1.1', 'P2.1', 'P3.1', 'P4.1', 'U.1', 'W1.2', 'B.1', 'M1.1', 'M1.2']
PLACED = ['P1.1', 'P2.1', 'P3.1', 'P4.1', 'W1.2', 'B.1', 'M1.1', 'M1.2']
INDEX = ['P1.1', 'P2.1', 'P3.1', 'P4.1', 'U.1', 'W1.1', 'B', 'X', 'M1']
VALUES = np.array([
    [.9, .4, 3e-8, .5, 0, 0],            # P1.1 target: weak preliminary EC and a hierarchical partial
    [.8, 0, 0, 0, 2e-4, .3],             # P2.1 context
    [1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9],  # P3.1 target, all weak, dense
    [0, 0, 0, 0, 0, 0],                  # P4.1 sparse zero
    [.2, .1, 0, 0, 0, 0],                # U.1 acquired, no coordinate
    [.7, 0, 0, 0, 0, 0],                 # W1.1: only W1.2 acquired, never moved to another version
    [.6, .5, 0, 0, 0, 0],                # B -> B.1 (unique version)
    [.4, 0, 0, 0, 0, 0],                 # X absent from the FASTA
    [.3, 0, 0, 0, 0, 0],                 # M1: two versions acquired, ambiguous
], dtype=np.float32)
MANIFEST_COLS = ['price_entry', 'host_assembly', 'proteome_protein_id', 'ground_truth_ecs', 'metabolic_stratum']


def write_layout(root, asm=ASM, fasta=FASTA, placed=PLACED, unavailable=None):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    seq = root / f'{asm}.sequence_hashes.tsv'
    seq.write_text('protein_id\tlength\tsequence_sha256\n' + ''.join(f'{p}\t1\t0\n' for p in fasta))
    lay = root / f'{asm}.layout.tsv'
    lay.write_text('protein_id\tcontig\tstart\tend\tstrand\n' + ''.join(f'{p}\tc\t1\t9\t+\n' for p in placed))
    report = {'schema': 'gspa.price_layout/1', 'all_proteins_sha256': 'x',
              'assemblies': {asm: {'proteins': len(fasta), 'unavailable_loci': unavailable or {'U.1': 'missing_CDS'},
                                   'sequence_hashes_sha256': digest(seq), 'layout_sha256': digest(lay)}}}
    (root / 'report.json').write_text(json.dumps(report))
    return str(root)


def write_matrix(path, index=INDEX, columns=COLUMNS, values=VALUES):
    np.savez(path, index=np.array(index), columns=np.array(columns), values=values)
    return str(path), digest(path)


def compact_stats(index, columns, values, entries, vocab):
    vset = {P.norm_ec(v) for v in vocab}
    mask = np.array([P.norm_ec(c) in vset for c in columns])
    out = {}
    for entry, mid in entries.items():
        row = values[index.index(mid)].astype(float)
        out[entry] = {'matched_id': mid, 'n_nonzero_vocab': int(((row != 0) & mask).sum()), 'max': float(row[mask].max())}
    return out


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.layout = write_layout(self.tmp / 'layout')
        self.parts = str(self.tmp / 'parts')

    def targets_file(self, predictors=tuple(P.PREDICTORS), entries=None, name='targets.json', stats_edit=None):
        entries = entries or {'T1': 'P1.1', 'T3': 'P3.1'}
        stats = compact_stats(INDEX, COLUMNS, VALUES, entries, COLUMNS)
        if stats_edit:
            stats_edit(stats)
        spec = {'schema': C.SCHEMA_TARGETS, 'layout_report_sha256': digest(Path(self.layout) / 'report.json'),
                'ec_vocab': COLUMNS, 'aliases': {ASM: dict(entries)}, 'conflicts': {},
                'targets': {p: {ASM: stats} for p in predictors}}
        path = self.tmp / name
        path.write_text(json.dumps(spec))
        return str(path)

    def matrix(self, name='m.npz', **kw):
        return write_matrix(self.tmp / name, **kw)

    def part(self, predictor='clean', retention=None, targets=None, parts=None, **kw):
        path, sha = kw.pop('matrix', None) or self.matrix(f'{predictor}.npz')
        return C.export_part(predictor, ASM, path, sha, self.layout, targets or self.targets_file(), parts or self.parts,
                             retention, **kw)

    @staticmethod
    def lines(rec):
        return [json.loads(x) for x in (Path(rec['part_dir']) / 'claims.jsonl').read_text().splitlines()]


class ExportPart(Base):
    def test_targets_keep_all_nonzero_context_all_nonzero_default(self):
        rec = self.part()
        claims = self.lines(rec)
        by = {}
        for c in claims:
            by.setdefault(c['protein_id'], {})[c['function_id']] = c['raw_score']
        self.assertEqual(set(by['P1.1']), {'1.1.1.1', '1.1.1.2', '1.1.1.n11'})  # preliminary kept, '1.1' excluded
        self.assertEqual(len(by['P3.1']), 6 - 1)                                  # all weak, partial '1.1' excluded
        self.assertNotIn('P4.1', by)                                              # sparse zero row has no claims
        self.assertEqual(rec['counts']['dropped_hierarchical_partial_nonzero'], 2)
        self.assertEqual(rec['counts']['rows_all_zero'], 1)
        self.assertTrue(rec['target_candidates']['full_candidate_coverage'])
        self.assertTrue(rec['context']['full_candidate_coverage'])
        self.assertTrue(rec['target_coverage']['verified'])

    def test_claim_lines_equal_the_existing_claim_serialisation(self):
        rec = self.part('graphec')
        text = (Path(rec['part_dir']) / 'claims.jsonl').read_text().splitlines()
        for line in text:
            c = json.loads(line)
            want = P._claim(c['protein_id'], c['function_id'], 'graphec', c['raw_score'], C.ORIGIN)
            self.assertEqual(line, json.dumps(want, sort_keys=True, separators=(',', ':')))
        self.assertEqual({json.loads(x)['evidence_type'] for x in text}, {'STRUCTURE_DEEPLEARNING'})

    def test_join_exact_first_bare_unique_only_and_unresolved_kept_intact(self):
        rec = self.part()
        ids = {c['protein_id'] for c in self.lines(rec)}
        self.assertIn('B.1', ids)          # bare id with one acquired version
        self.assertNotIn('B', ids)
        self.assertIn('W1.1', ids)         # wrong version is NOT moved to W1.2
        self.assertNotIn('W1.2', ids)
        self.assertIn('M1', ids)           # ambiguous bare stays as is
        self.assertIn('X', ids)
        counts = rec['join_map']['status_counts']
        self.assertEqual(counts, {'exact': 5, 'bare_unique_version': 1, 'unresolved_absent': 1,
                                  'unresolved_ambiguous_bare': 1, 'unresolved_version_mismatch': 1})
        table = (Path(rec['part_dir']) / 'join_map.tsv').read_text()
        self.assertIn('W1.1\tW1.1\tunresolved_version_mismatch\tnot_acquired\t1', table)
        self.assertIn('U.1\tU.1\texact\tacquired_unplaced\t1', table)

    def test_unplaced_only_from_report_or_explicit_declaration(self):
        rec = self.part(unplaced_declarations={ASM: ['X']})
        self.assertEqual(rec['unplaced']['acquired_unplaced_in_claims'], ['U.1'])
        self.assertEqual(rec['unplaced']['declared_unplaced_in_claims'], ['X'])
        self.assertEqual(rec['unplaced']['undeclared_unresolved_in_claims'], ['M1', 'W1.1'])

    def test_context_pruning_is_explicit_and_never_full_coverage(self):
        rec = self.part(retention=C.Retention(min_score=0.5, top_k=1))
        ctx, tgt = rec['counts']['context'], rec['counts']['targets']
        self.assertGreater(ctx['pruned_away'], 0)
        self.assertTrue(rec['context']['pruned'])
        self.assertFalse(rec['context']['full_candidate_coverage'])
        self.assertTrue(rec['target_candidates']['full_candidate_coverage'])   # targets are never pruned
        self.assertEqual(tgt['retained'], tgt['nonzero_before'])
        self.assertNotEqual(rec['claims']['sha256'], rec['claims_before_retention_sha256'])
        self.assertEqual(rec['retention']['policy'], 'threshold_topk')

    def test_default_retention_before_equals_after_hash(self):
        rec = self.part()
        self.assertEqual(rec['claims']['sha256'], rec['claims_before_retention_sha256'])
        self.assertEqual(rec['claims']['sha256'], digest(Path(rec['part_dir']) / 'claims.jsonl'))

    def test_invalid_scores_refused_and_nothing_left_behind(self):
        for bad in (float('nan'), float('inf'), -0.1, 1.5):
            values = VALUES.copy()
            values[1, 0] = bad
            with self.assertRaisesRegex(C.ExportError, 'NaN/inf/negative'):
                self.part(matrix=self.matrix('bad.npz', values=values), parts=str(self.tmp / 'bad_parts'))
            self.assertEqual([p.name for p in (self.tmp / 'bad_parts' / 'clean').iterdir()], [])

    def test_duplicate_axes_alias_collision_and_ec_collision(self):
        with self.assertRaisesRegex(C.ExportError, 'duplicate protein ids'):
            C.build_join_map(['a.1', 'a.1'], {'a.1'})
        with self.assertRaisesRegex(C.ExportError, 'alias collision'):
            C.build_join_map(['B', 'B.1'], {'B.1'})
        with self.assertRaisesRegex(C.ExportError, 'duplicate normalised EC'):
            C.column_plan(['EC:1.1.1.1', '1.1.1.1'])
        with self.assertRaisesRegex(C.ExportError, 'invalid EC'):
            C.column_plan(['1.1.1.1', 'banana'])

    def test_completed_part_is_immutable(self):
        self.part()
        with self.assertRaisesRegex(C.ExportError, 'immutable'):
            self.part()

    def test_wrong_hash_refused_before_pickle_is_loaded(self):
        marker = self.tmp / 'ran'

        class Evil:
            def __reduce__(self):
                return (Path(marker).write_text, ('loaded',))

        path = self.tmp / 'm.pkl'
        path.write_bytes(pickle.dumps(Evil()))
        with self.assertRaisesRegex(C.ExportError, 'refusing to load'):
            C.open_matrix(str(path), '0' * 64, trust_pickle=True)
        self.assertFalse(marker.exists())
        art = self.tmp / 'a.pkl.gz'
        with gzip.open(art, 'wb') as fh:
            fh.write(pickle.dumps(Evil()))
        with self.assertRaisesRegex(C.ExportError, 'refusing to unpickle'):
            C.build_targets(str(art), '0' * 64, self.layout)
        self.assertFalse(marker.exists())

    def test_missing_full_matrix_is_an_error_not_a_compact_fallback(self):
        with self.assertRaisesRegex(C.ExportError, 'no compact fallback'):
            C.export_part('clean', ASM, str(self.tmp / 'absent.npz'), '0' * 64, self.layout, self.targets_file(), self.parts)
        self.assertFalse((Path(self.parts) / 'clean' / ASM).exists())

    def test_pickle_dataframe_needs_trust_and_loads(self):
        import pandas as pd
        df = pd.DataFrame(VALUES, index=INDEX, columns=COLUMNS)
        path = self.tmp / 'df.pkl'
        df.to_pickle(path)
        with self.assertRaisesRegex(C.ExportError, 'trust_pickle'):
            C.open_matrix(str(path), digest(path))
        rec = self.part(matrix=(str(path), digest(path)), trust_pickle=True)
        self.assertTrue(rec['target_coverage']['verified'])

    def test_interrupted_stream_leaves_no_part(self):
        class Broken(C.ArraySource):
            def rows(self):
                for i, row in super().rows():
                    if i == 3:
                        raise KeyboardInterrupt
                    yield i, row
        with self.assertRaises(KeyboardInterrupt):
            C.export_part_from_source(Broken(INDEX, COLUMNS, VALUES), 'clean', ASM, self.layout, self.targets_file(), self.parts)
        left = list((Path(self.parts) / 'clean').iterdir())
        self.assertEqual(left, [])
        rec = self.part()                                                   # a clean rerun is possible
        self.assertEqual(rec['status'], 'complete')

    def test_stream_is_memory_bounded(self):
        rows, cols = 4000, 1500
        columns = [f'1.1.{i // 100 + 1}.{i % 100 + 1}' for i in range(cols)]
        fasta = [f'Q{i}.1' for i in range(rows)]
        lay = write_layout(self.tmp / 'big_layout', 'GCA_BIG.1', fasta, fasta, {})

        class Gen:
            def __init__(self, index, columns):
                self.index, self.columns = index, columns

            def rows(self):
                rng = np.random.default_rng(0)
                for i in range(rows):
                    row = np.zeros(cols, dtype=np.float32)
                    row[rng.integers(0, cols, 25)] = rng.random(25, dtype=np.float32) * .9 + .01
                    yield i, row

        spec = {'schema': C.SCHEMA_TARGETS, 'layout_report_sha256': digest(Path(lay) / 'report.json'), 'ec_vocab': None,
                'aliases': {}, 'conflicts': {}, 'targets': {}}
        (self.tmp / 'big.json').write_text(json.dumps(spec))
        tracemalloc.start()
        rec = C.export_part_from_source(Gen(fasta, columns), 'clean', 'GCA_BIG.1', lay, str(self.tmp / 'big.json'), self.parts)
        _cur, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        self.assertGreater(rec['claims']['claims'], 90000)
        self.assertGreater(rec['claims']['bytes'], 15_000_000)
        self.assertLess(peak, 6_000_000)       # the dense matrix would need 24 MB; the file is not held either

    def test_export_part_has_no_truth_parameters(self):
        import inspect
        for fn in (C.export_part, C.export_part_from_source, C._stream):
            names = ' '.join(inspect.signature(fn).parameters)
            for forbidden in ('manifest', 'label', 'truth', 'ground'):
                self.assertNotIn(forbidden, names)

    def test_target_coverage_mismatch_with_compact_stats_is_recorded(self):
        def edit(stats):
            stats['T1']['n_nonzero_vocab'] += 1
        rec = self.part(targets=self.targets_file(stats_edit=edit))
        self.assertFalse(rec['target_coverage']['verified'])
        self.assertEqual(rec['target_coverage']['mismatches'][0]['entry'], 'T1')

    def test_absent_target_unverifies_coverage(self):
        rec = self.part(targets=self.targets_file(entries={'T1': 'P1.1', 'T9': 'P3.1'},
                                                  stats_edit=lambda s: s['T9'].update(matched_id='GONE.1')))
        self.assertEqual(rec['target_coverage']['targets_absent_from_matrix'], ['T9'])
        self.assertFalse(rec['target_coverage']['verified'])


class Merge(Base):
    ROWS = [('T1', ASM, 'P1.1', 'EC:1.1.1.1;1.1.1.-'), ('T3', ASM, 'P3.1', '3.1.1.1')]

    def manifest(self, rows=None, name='manifest.tsv'):
        path = self.tmp / name
        path.write_text('\t'.join(MANIFEST_COLS) + '\n' + ''.join(
            '\t'.join([e, a, p, t, 'S']) + '\n' for e, a, p, t in (rows or self.ROWS)))
        return str(path)

    def build_parts(self, predictors=P.POOL4, **kw):
        t = kw.pop('targets', None) or self.targets_file()
        for p in predictors:
            self.part(p, targets=t, **kw)
        return t

    def merge(self, t, rows=None, out='merged', variants=('clean', 'ensemble_pool4'), manifest=None):
        return C.merge_parts(self.parts, t, manifest or self.manifest(rows), str(self.tmp / out), variants)

    def test_merge_pools_validates_and_writes_labels_and_manifests(self):
        t = self.build_parts()
        rep = self.merge(t)
        self.assertEqual(rep['status'], 'complete')
        pooled = Path(rep['claims']['ensemble_pool4'][ASM]['path'])
        n = sum(rep['parts'][f'{p}/{ASM}']['claims'] for p in P.POOL4)
        self.assertEqual(len(pooled.read_text().splitlines()), n)
        keys = {(c['source'], c['protein_id'], c['function_id']) for c in map(json.loads, pooled.read_text().splitlines())}
        self.assertEqual(len(keys), n)
        self.assertEqual(rep['claims']['clean'][ASM]['kind'], 'part (referenced, not copied)')
        man = json.loads(Path(rep['manifests']['ensemble_pool4']).read_text())
        self.assertEqual(man['genomes'][0]['claims'], str(pooled))
        self.assertIn('keeps every nonzero score', man['description'])
        self.assertIn('Target candidates: all nonzero predictor scores', man['description'])
        self.assertTrue(rep['target_candidates_full'] and rep['context_full_candidate_coverage'])
        self.assertEqual(json.loads((self.tmp / 'merged' / 'unplaced_proteins.json').read_text())[ASM], ['U.1'])
        self.assertIn('M1', rep['undeclared_unresolved_claim_ids'][ASM])

    def test_ec_prefix_normalisation_same_names_in_labels_and_claims(self):
        t = self.build_parts(('clean',))
        rep = C.merge_parts(self.parts, t, self.manifest(), str(self.tmp / 'm1'), ('clean',))
        label = (self.tmp / 'm1' / 'labels' / f'{ASM}.ec.tsv').read_text().splitlines()
        self.assertEqual(label, ['protein_id\tec', 'P1.1\t1.1.1.1', 'P3.1\t3.1.1.1'])   # 'EC:' dropped, partial dropped
        claims = {(c['protein_id'], c['function_id']) for c in map(json.loads, Path(rep['claims']['clean'][ASM]['path']).read_text().splitlines())}
        self.assertLessEqual({tuple(x.split('\t')) for x in label[1:]}, claims)

    def test_truth_mutation_changes_labels_never_claims(self):
        t = self.build_parts()
        a = self.merge(t, out='a')
        flipped = [('T1', ASM, 'P1.1', '9.9.9.9'), ('T3', ASM, 'P3.1', '1.1.1.2;2.7.1.1')]
        b = self.merge(t, rows=flipped, out='b', manifest=self.manifest(flipped, 'm2.tsv'))
        for v in a['claims']:
            self.assertEqual(a['claims'][v][ASM]['sha256'], b['claims'][v][ASM]['sha256'])
        self.assertEqual({k: r['claims_sha256'] for k, r in a['parts'].items()}, {k: r['claims_sha256'] for k, r in b['parts'].items()})
        self.assertNotEqual(a['labels'][ASM]['sha256'], b['labels'][ASM]['sha256'])

    def test_merge_refuses_until_all_selected_parts_complete(self):
        t = self.build_parts(('clean', 'dpz', 'mapred'))
        with self.assertRaisesRegex(C.ExportError, 'not complete.*enzbert'):
            self.merge(t)
        self.assertFalse((self.tmp / 'merged').exists())
        self.part('enzbert', targets=t)
        self.merge(t)

    def test_merge_refuses_tampered_claims_and_mixed_retention_and_unverified_coverage(self):
        t = self.build_parts(('clean',))
        f = Path(self.parts) / 'clean' / ASM / 'claims.jsonl'
        f.write_text(f.read_text().replace('"raw_score":', '"raw_score":0', 1))
        with self.assertRaisesRegex(C.ExportError, 'differs from its recorded'):
            self.merge(t, variants=('clean',))
        shutil.rmtree(self.parts)
        self.part('clean', targets=t)
        self.part('dpz', targets=t, retention=C.Retention(top_k=1))
        with self.assertRaisesRegex(C.ExportError, 'retention_sha256'):
            self.merge(t, variants=('clean', 'dpz'))
        shutil.rmtree(self.parts)
        bad = self.targets_file(name='bad.json', stats_edit=lambda s: s['T3'].update(max=0.5))
        self.part('clean', targets=bad)
        with self.assertRaisesRegex(C.ExportError, 'not verified'):
            C.merge_parts(self.parts, bad, self.manifest(), str(self.tmp / 'm3'), ('clean',))

    def test_pruned_context_is_labelled_in_manifest_and_report(self):
        t = self.build_parts(('clean',), retention=C.Retention(min_score=0.5))
        rep = C.merge_parts(self.parts, t, self.manifest(), str(self.tmp / 'm4'), ('clean',))
        self.assertFalse(rep['context_full_candidate_coverage'])
        man = json.loads(Path(rep['manifests']['clean']).read_text())
        self.assertIn('PRUNED', man['description'])
        self.assertEqual(man['declared_controls']['context_retention']['min_score'], 0.5)

    def test_merge_output_is_immutable(self):
        t = self.build_parts(('clean',))
        C.merge_parts(self.parts, t, self.manifest(), str(self.tmp / 'm5'), ('clean',))
        with self.assertRaisesRegex(C.ExportError, 'immutable'):
            C.merge_parts(self.parts, t, self.manifest(), str(self.tmp / 'm5'), ('clean',))

    def test_scan_rejects_duplicate_ec_and_split_protein_blocks(self):
        rec = self.part('clean')
        f = Path(rec['part_dir']) / 'claims.jsonl'
        lines = f.read_text().splitlines(True)
        dup = self.tmp / 'dup.jsonl'
        dup.write_text(lines[0] + lines[0])
        with self.assertRaisesRegex(C.ExportError, 'duplicate EC'):
            C._scan_claims(dup, 'clean')
        split = self.tmp / 'split.jsonl'
        other = next(x for x in lines if json.loads(x)['protein_id'] != json.loads(lines[0])['protein_id'])
        split.write_text(lines[0] + other + lines[1])
        with self.assertRaisesRegex(C.ExportError, 'two blocks'):
            C._scan_claims(split, 'clean')


class Cli(Base):
    def test_cli_export_part_with_verification_and_missing_matrix(self):
        path, sha = self.matrix()
        ver = self.tmp / 'verification.json'
        ver.write_text(json.dumps({'files': [{'predictor': 'clean', 'assembly': ASM, 'destination': path,
                                              'destination_sha256': sha, 'verified': True}]}))
        args = ['export-part', '--predictor', 'clean', '--assembly', ASM, '--verification', str(ver), '--layout', self.layout,
                '--targets', self.targets_file(), '--out', self.parts]
        self.assertEqual(C.main(args), 0)
        self.assertEqual(C.main(args), 2)                       # immutable
        ver.write_text(json.dumps({'files': [{'predictor': 'clean', 'assembly': ASM, 'destination': path + '.gone',
                                              'destination_sha256': sha, 'verified': True}]}))
        self.assertEqual(C.main(args[:-1] + [str(self.tmp / 'other')]), 2)


if __name__ == '__main__':
    unittest.main()
