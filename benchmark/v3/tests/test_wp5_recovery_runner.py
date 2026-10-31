"""WP5 step 4 (design §5, §8.4; rulings §13.5, §13.6): runner on the production selection loop.

The engine's per-cell LP (``simulate_cell``) is replaced by a fake that books the same solver calls;
selection, budgets, masking re-check and the Qwen loop are the production code.
"""
import copy
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wp5_fixtures as F  # noqa: E402
import prepare_wp5_cases as P5  # noqa: E402
import run_wp5_recovery as R  # noqa: E402
import metabolic_experiments as ME  # noqa: E402
import masked_assignment_experiments as MAE  # noqa: E402

PROBLEM = 'W001-manuscript_regulatory'


class Base(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix='wp5-run-', dir=os.environ.get('TMPDIR')))
        self.addCleanup(shutil.rmtree, self.root, True)
        self.world = F.runner_world(self.root)
        self.public, self.record = self.world['public'], self.world['record']
        self.out = self.root / 'runs'
        patcher = F.engine_patch()
        patcher.start()
        self.addCleanup(patcher.stop)

    def model(self):
        return ME.load_sbml_model(self.public / 'models/W001.xml')

    def run_stages(self, chat=None, stages=('exhaustive', 'deterministic', 'qwen')):
        R.make_protocol(self.public, self.out)
        for stage in stages:
            R.execute(self.public, self.out, stage, chat_factory=(lambda: chat) if chat else None)
        return {p.relative_to(self.out / 'runs').as_posix(): P5.read_json(p)
                for p in sorted((self.out / 'runs').rglob('*.json'))}


class ProductionEngineRoute(Base):
    def test_every_arm_runs_through_run_assignment_experiments_without_the_target(self):
        calls = []
        real = MAE.run_assignment_experiments

        def spy(manifest, original, overrides=None, chat=None):
            calls.append((manifest['assignment_recovery']['mask_gene'], dict(overrides)))
            self.assertNotIn('G_b0001', original.gene_alias)
            self.assertNotIn('b0001', json.dumps(manifest))
            return real(manifest, original, overrides, chat)
        with mock.patch.object(MAE, 'run_assignment_experiments', spy):
            runs = self.run_stages(chat=F.reply_chat(['{"experiment_id": "C2"}'] * 20))
        self.assertEqual(len(calls), len(runs))
        self.assertEqual({c[0] for c in calls}, {R.PLACEHOLDER})
        self.assertEqual({r['engine_entry'] for r in runs.values()}, {R.ENGINE_ENTRY})
        self.assertEqual(len(runs), 1 + 2 * (3 + 10) + 2)

    def test_placeholder_masking_is_a_no_op_on_the_public_fingerprint(self):
        model = self.model()
        fp = 'masked:' + MAE.fingerprint(model)
        self.assertEqual(fp, self.record['public_problem']['model'])
        masked = MAE.mask_model(R.placeholder_model(model), R.PLACEHOLDER)
        self.assertEqual(masked.path, fp)
        with self.assertRaises(P5.WP5Error):
            R.placeholder_model(R.placeholder_model(model))

    def test_drifted_masked_model_or_problem_is_refused(self):
        model = self.model()
        record = copy.deepcopy(self.record)
        record['public_problem']['hypotheses'][1]['evidence']['sequence'] = 0.1
        with self.assertRaisesRegex(P5.WP5Error, 'differs from the frozen protocol'):
            R.run_arm(record, model, 'evidence', 3)
        drifted = copy.deepcopy(model)
        drifted.gpr = dict(drifted.gpr, R_DECOY1=[['G_b0003']])
        with self.assertRaisesRegex(P5.WP5Error, 'different public problem'):
            R.run_arm(self.record, drifted, 'evidence', 3)


class EqualBudgetsAndHash(Base):
    def test_identical_public_problem_hash_and_equal_caps_across_arms(self):
        runs = self.run_stages(chat=F.reply_chat(['{"experiment_id": "C2"}'] * 20))
        protocol = P5.read_json(self.out / 'protocol.json')
        entry = protocol['problems'][PROBLEM]
        c = entry['cost_per_condition']
        self.assertEqual(c, 2 + 4 + 4)
        self.assertEqual(c, MAE.condition_cost(self.record['public_problem'], self.record['public_problem']['conditions'][0]))
        self.assertEqual({r['public_problem_sha256'] for r in runs.values()}, {self.record['public_problem_sha256']})
        for name, run in runs.items():
            self.assertEqual(run['protocol_sha256'], P5.sha_file(self.out / 'protocol.json'))
            if run['arm'] == 'exhaustive':
                self.assertEqual((run['max_experiments'], run['max_solver_calls']), (4, 4 * c))
                self.assertEqual(run['visited_conditions'], ['C0', 'C1', 'C2', 'C3'])
                continue
            k = run['k']
            self.assertEqual((run['max_experiments'], run['max_solver_calls']), (k, k * c), name)
            self.assertEqual((run['config']['max_experiments'], run['config']['max_solver_calls']), (k, k * c))
            self.assertEqual(run['solver_calls_used'] + run['solver_calls_unused'], k * c)
            self.assertEqual(run['solver_calls_used'], c * len(run['visited_conditions']))
            self.assertLessEqual(len(run['visited_conditions']), k)

    def test_minimal_medium_visits_only_the_baseline_and_fixed_panel_follows_the_table(self):
        runs = self.run_stages(stages=('exhaustive', 'deterministic'))
        for k in (3, 6):
            mm = runs[f'{PROBLEM}/minimal_medium-k{k}-s0.json']
            self.assertEqual(mm['visited_conditions'], ['C0'])
            self.assertEqual(mm['solver_calls_unused'], (k - 1) * 10)   # never reallocated
            fp = runs[f'{PROBLEM}/fixed_panel-k{k}-s0.json']
            self.assertEqual(fp['visited_conditions'], ['C0', 'C1', 'C2'])
            ev = runs[f'{PROBLEM}/evidence-k{k}-s0.json']
            self.assertEqual(ev['visited_conditions'][0], 'C3')          # the only evidence-supported condition

    def test_arm_cells_must_equal_the_exhaustive_reference(self):
        self.run_stages(stages=('exhaustive',))
        path = self.out / 'runs' / R.run_name(PROBLEM, 'exhaustive', 0, 0)
        exhaustive = P5.read_json(path)
        exhaustive['cells'][0]['contrast'] = 'knockout_no_growth'
        P5.write_json(path, exhaustive)
        with self.assertRaisesRegex(P5.WP5Error, 'differ from the exhaustive'):
            R.execute(self.public, self.out, 'deterministic')
        mm = P5.read_json(self.out / 'runs' / R.run_name(PROBLEM, 'minimal_medium', 3, 0))
        self.assertEqual(mm['status'], 'cell_mismatch')

    def test_deterministic_stage_needs_the_exhaustive_reference_first(self):
        R.make_protocol(self.public, self.out)
        with self.assertRaisesRegex(P5.WP5Error, 'exhaustive stage first'):
            R.execute(self.public, self.out, 'deterministic')

    def test_shards_partition_the_runs(self):
        R.make_protocol(self.public, self.out)
        protocol = P5.read_json(self.out / 'protocol.json')
        runs = R.plan_runs(protocol, 'deterministic')
        shards = [R.shard_of(runs, f'{i}/4') for i in range(4)]
        self.assertEqual(sorted(x for s in shards for x in s), sorted(runs))
        with self.assertRaises(P5.WP5Error):
            R.shard_of(runs, '4/4')


class QwenArm(Base):
    def test_requests_are_aliased_and_carry_no_gene_id_organism_or_strain(self):
        seen = []
        replies = ['{"experiment_id": "%s"}' % c for c in ('C2', 'C1', 'C0') + ('C2', 'C1', 'C0', 'C3')]
        runs = self.run_stages(chat=F.reply_chat(replies, seen))
        run = runs[f'{PROBLEM}/qwen-k3-s0.json']
        self.assertEqual(run['status'], 'completed')
        self.assertEqual(run['visited_conditions'], ['C2', 'C1', 'C0'])
        self.assertEqual(runs[f'{PROBLEM}/qwen-k6-s0.json']['visited_conditions'], ['C2', 'C1', 'C0', 'C3'])
        self.assertTrue(seen)
        text = json.dumps(seen)
        for leak in ('G_b0003', 'G_b0005', 'b0003', 'G_b0001', 'b0001', 'pgiX', 'coli', 'K-12', 'MG1655', 'BW25113'):
            self.assertNotIn(leak, text)
        self.assertIn('MASKED_TARGET', text)
        self.assertRegex(text, r'\bP[0-9]+\b')
        self.assertRegex(text, r'\bH[0-9]+\b')
        self.assertRegex(text, r'\bC[0-9]+\b')
        transcript = P5.read_jsonl(self.out / run['transcript'])
        self.assertEqual(run['transcript_sha256'], P5.sha_file(self.out / run['transcript']))
        self.assertEqual(len(transcript), run['agent_usage']['requests'])
        self.assertEqual(run['agent_usage']['transport']['prompt_tokens'], 10)

    def test_first_engine_prompt_equals_the_audit_reconstruction(self):
        seen = []
        model = self.model()
        R.run_arm(self.record, model, 'qwen', 3, chat=F.reply_chat(['{"experiment_id": "C1"}'] * 5, seen))
        first = json.dumps(seen[0], sort_keys=True).encode()
        public = self.record['public_problem']
        self.assertEqual(first, R.first_prompt(public, 3))
        self.assertEqual(first, R.first_prompt(public, 3, model.gene_alias))

    def test_agent_failure_or_invalid_reply_stops_the_arm_without_credit(self):
        model = self.model()
        for replies in ([RuntimeError('endpoint down')], ['not json'], ['{"experiment_id": "C99"}'],
                        ['{"experiment_id": "C1"}', TimeoutError('slow')]):
            run = R.run_arm(self.record, model, 'qwen', 3, chat=F.reply_chat(replies))
            self.assertEqual(run['status'], 'agent_failed', replies)
            self.assertEqual(run['selection'][-1]['status'], 'agent_failed')
            self.assertEqual(len(run['visited_conditions']), len(replies) - 1)
            self.assertEqual(run['agent_usage']['requests'], len(replies))

    def test_qwen_transport_settings_are_frozen(self):
        self.assertEqual(R.QWEN['timeout_s'], 600)
        self.assertEqual((R.QWEN['temperature'], R.QWEN['seed'], R.QWEN['enable_thinking']), (0.0, 0, False))
        captured = {}

        class Response:
            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def read(self):
                return json.dumps({'choices': [{'message': {'content': '{}'}}], 'usage': {}}).encode()

        def urlopen(request, timeout):
            captured.update(json.loads(request.data), timeout=timeout)
            return Response()
        with mock.patch('urllib.request.urlopen', urlopen), mock.patch.object(ME, 'load_api_key', lambda: 'k'):
            R.qwen_chat()([{'role': 'user', 'content': 'x'}])
        self.assertEqual((captured['timeout'], captured['temperature'], captured['seed']), (600, 0.0, 0))
        self.assertEqual(captured['chat_template_kwargs'], {'enable_thinking': False})


class PublicInputsOnly(Base):
    def test_no_private_file_is_opened_by_any_stage(self):
        with F.refuse_paths(F.no_private) as opened:
            self.run_stages(chat=F.reply_chat(['{"experiment_id": "C2"}'] * 20))
        self.assertTrue(opened)
        self.assertFalse([p for p in opened if F.no_private(p)])

    def test_private_paths_and_unaudited_or_changed_public_inputs_are_refused(self):
        with self.assertRaisesRegex(P5.WP5Error, 'only reads public'):
            R.load_public(self.world['bundle'] / 'private')
        with self.assertRaisesRegex(P5.WP5Error, 'only reads public'):
            R.main(['protocol', '--public', str(self.public), '--out', str(self.world['bundle'] / 'private/runs')])
        audit = self.public / 'masking_audit.json'
        report = P5.read_json(audit)
        P5.write_json(audit, dict(report, status='fail'))
        with self.assertRaisesRegex(P5.WP5Error, 'masking audit'):
            R.make_protocol(self.public, self.out)
        audit.unlink()
        with self.assertRaisesRegex(P5.WP5Error, 'no masking audit'):
            R.make_protocol(self.public, self.out)
        F.pass_audit(self.world['bundle'])
        model = self.public / 'models/W001.xml'
        model.write_text(model.read_text() + '\n')
        with self.assertRaisesRegex(P5.WP5Error, 'changed after construction'):
            R.make_protocol(self.public, self.out)

    def test_protocol_is_frozen_and_bound_to_the_implementation(self):
        protocol = R.make_protocol(self.public, self.out)
        self.assertEqual(protocol['engine_entry'], R.ENGINE_ENTRY)
        self.assertEqual(protocol['budgets'], [3, 6])
        self.assertEqual(protocol['problems'][PROBLEM]['panels'],
                         {'minimal_medium': ['C0'], 'fixed_panel': ['C0', 'C1', 'C2']})
        self.assertIn('production selection loop', (self.out / 'README.md').read_text())
        self.assertEqual(R.make_protocol(self.public, self.out), protocol)
        with mock.patch.object(R, 'implementation_hashes', lambda: {'x': 'changed'}):
            with self.assertRaisesRegex(P5.WP5Error, 'implementation changed'):
                R.execute(self.public, self.out, 'exhaustive')


class AuditHelperParity(Base):
    """``planned_selection`` (audit-only, LP-free) equals the engine's visits for fixed, random, evidence."""

    def test_planned_selection_matches_the_engine(self):
        model = self.model()
        public = self.record['public_problem']
        panels = R.panels(self.record)
        for k in R.BUDGETS:
            for arm, seeds in (('minimal_medium', [0]), ('fixed_panel', [0]), ('evidence', [0]),
                               ('random', list(R.RANDOM_SEEDS))):
                for seed in seeds:
                    run = R.run_arm(self.record, model, arm, k, seed)
                    self.assertEqual(run['visited_conditions'],
                                     R.planned_selection(public, arm, k, seed, panels.get(arm)), (arm, k, seed))
        with self.assertRaises(P5.WP5Error):
            R.planned_selection(public, 'qwen', 3)


if __name__ == '__main__':
    unittest.main()
