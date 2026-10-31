import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import check_source_arm_context_models as C


class SourceModelProjectionTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home() / 'tmp')))
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.model = self.root / 'model.json'
        self.model.write_text(json.dumps(dict(contract={'features': ['fn.in_complex_term', 'seq.max_raw_any']},
                                             display_threshold={'value': .5})))
        self.sha = C.S.sha(self.model)
        self.old, self.new, self.records, self.outcomes, self.before, self.after = [], [], [], [], {}, {}

    def candidate(self, protein, kind='EC', term='EC:1.1.1.1', role='CATALYTIC', probability=.7,
                  unsupported=False, complex_term=False, prior=.4, action='rescored', projected=True):
        key = '|'.join((protein, kind, term))
        ck = key + '|' + role
        row = dict(genome_id='g', protein_id=protein, function_type=kind, function_id=term,
                   function_key=key, candidate_key=ck, role=role, sequence_sources=['clean'],
                   sequence_contracts=[], cohort_context_mode='UNIVERSAL',
                   features={'fn.in_complex_term': float(complex_term), 'seq.max_raw_any': .3})
        reason = 'fixture_unsupported' if unsupported else None
        score = dict(sha256=self.sha, score_status='unsupported' if unsupported else 'scored',
                     unsupported_reason=reason, probability=None if unsupported else probability,
                     margin=None if unsupported else 1., above_display_threshold=None if unsupported else probability >= .5)
        digest = hashlib.sha256(b''.join(struct.pack('>d', row['features'][n])
                                         for n in ('fn.in_complex_term', 'seq.max_raw_any'))).hexdigest()
        record = dict(candidate_key=ck, function_key=key, role=role, sequence_sources=['clean'],
                      sequence_contracts=[], projected=projected, feature_digest=digest,
                      locus_probability=None, locus_null_probability=None, **{k:v for k,v in score.items() if k != 'sha256'})
        self.old.append(row)
        self.new.append(dict(copy.deepcopy(row), model=score))
        self.records.append(record)
        existing = next((o for o in self.outcomes if o['function_key'] == key), None)
        if existing is None:
            existing = dict(function_key=key, action=action, reason=None, prior_probability=prior,
                            projected_role=role if projected else None,
                            learned_probability=probability if projected else None, role_probabilities={})
            self.outcomes.append(existing)
        if not unsupported:
            existing['role_probabilities'][role] = probability
        if prior is not None:
            self.before[key] = prior
        if action in ('rescored', 'added'):
            self.after[key] = probability if projected else self.after.get(key, prior)
        elif action in ('unchanged_unsupported', 'scored_only') and prior is not None:
            self.after[key] = prior
        return record, existing

    def save(self):
        def rows(path, values):
            path.write_text(''.join(json.dumps(r) + '\n' for r in values))
        rows(self.root / 'old.jsonl', self.old)
        rows(self.root / 'new.jsonl', self.new)
        for name, values in (('before.tsv', self.before), ('after.tsv', self.after)):
            with (self.root / name).open('w') as stream:
                writer = csv.writer(stream, delimiter='\t')
                writer.writerow(['protein_id', 'type', 'function_id', 'posterior_prob'])
                for key, value in values.items():
                    writer.writerow([*key.split('|'), value])
        counts = dict(candidates=len(self.records), unsupported=sum(r['score_status']=='unsupported' for r in self.records))
        counts.update({a:sum(o['action']==a for o in self.outcomes) for a in
                       ('rescored', 'added', 'unchanged_unsupported', 'retracted', 'removed_no_projection', 'scored_only')})
        report = dict(schema='gspa.learned-stage-report/1', model={'sha256': self.sha},
                      candidates=self.records, outcomes=self.outcomes, counts=counts,
                      settings=dict(retract_below=None, add_unretained_claims=False, add_context_candidates=True))
        (self.root / 'after.tsv.learned.json').write_text(json.dumps(report))

    def audit(self):
        self.save()
        # This fixture tests the projection auditor, not a fitted model. Existing
        # trainer/JVM tests cover numerical implementations on real frozen models.
        reasons = {tuple(r['features'].values()): self.records[i]['unsupported_reason'] for i,r in enumerate(self.new)}
        with patch.object(C.F.T, 'unsupported_input_reason', side_effect=lambda m,v,*a: reasons[tuple(v)]), \
             patch.object(C.F.T, 'predict_proba', return_value=.7):
            return C.audit_genome(self.root/'before.tsv', self.root/'after.tsv', self.root/'old.jsonl',
                                  self.root/'new.jsonl', self.model)

    def test_rescores_and_new_admissions_are_both_accounted_for(self):
        self.candidate('retained')
        self.candidate('new', prior=None, action='added')
        result = self.audit()
        self.assertEqual(result['added'], 1)
        self.assertEqual(result['actions'], {'rescored':1, 'added':1})

    def test_display_threshold_tie_preserves_jvm_decision(self):
        self.model.write_text(json.dumps(dict(contract={'features': ['fn.in_complex_term', 'seq.max_raw_any']},
                                             display_threshold={'value': .7})))
        self.sha = C.S.sha(self.model)
        record, _ = self.candidate('p', probability=.7)
        record['above_display_threshold'] = False
        self.new[0]['model']['above_display_threshold'] = False
        result = self.audit()
        self.assertEqual(result['display_ties'], 1)
        self.assertEqual(result['display_tie_candidates'], [record['candidate_key']])

    def test_display_threshold_tie_still_requires_jvm_flags_to_agree(self):
        self.model.write_text(json.dumps(dict(contract={'features': ['fn.in_complex_term', 'seq.max_raw_any']},
                                             display_threshold={'value': .7})))
        self.sha = C.S.sha(self.model)
        record, _ = self.candidate('p', probability=.7)
        record['above_display_threshold'] = False
        self.new[0]['model']['above_display_threshold'] = True
        with self.assertRaisesRegex(ValueError, 'display decision differs'):
            self.audit()

    def test_non_tie_display_mismatch_is_rejected(self):
        record, _ = self.candidate('p', probability=.7)
        record['above_display_threshold'] = False
        self.new[0]['model']['above_display_threshold'] = False
        with self.assertRaisesRegex(ValueError, 'display decision differs'):
            self.audit()

    def test_accessory_ec_does_not_replace_catalytic_projection(self):
        self.candidate('p')
        self.candidate('p', role='ACCESSORY', projected=False)
        self.assertEqual(self.audit()['functions'], 1)
        self.outcomes[0]['projected_role'] = 'ACCESSORY'
        with self.assertRaisesRegex(ValueError, 'projected role differs'):
            self.audit()

    def test_ec_with_only_accessory_role_removes_legacy_catalytic_assignment(self):
        self.candidate('p', role='ACCESSORY', action='removed_no_projection', projected=False)
        self.assertEqual(self.audit()['removed'], 1)
        self.after.update(self.before)
        with self.assertRaisesRegex(ValueError, 'not removed'):
            self.audit()

    def test_go_complex_projects_subunit_role(self):
        self.candidate('p', kind='GO', term='GO:0000001', role='COMPLEX_SUBUNIT', complex_term=True)
        self.assertEqual(self.audit()['actions'], {'rescored':1})

    def test_unsupported_candidate_preserves_rule_fallback(self):
        self.candidate('p', unsupported=True, action='unchanged_unsupported', projected=False)
        self.assertEqual(self.audit()['reasons'], {'fixture_unsupported':1})
        self.after[next(iter(self.after))] = .8
        with self.assertRaisesRegex(ValueError, 'fallback prediction changed'):
            self.audit()

    def test_taxon_forbidden_go_keeps_preceding_output(self):
        _, outcome = self.candidate('p', kind='GO', term='GO:0000001', role='UNSPECIFIED', action='scored_only')
        outcome['reason'] = C.TAXON_REASON
        self.assertEqual(self.audit()['actions'], {'scored_only':1})

    def test_unreported_output_change_rejected(self):
        self.candidate('p')
        self.after['stray|EC|EC:2.2.2.2'] = .7
        with self.assertRaisesRegex(ValueError, 'unreported output change'):
            self.audit()

    def test_input_mutation_rejected_even_with_equal_score(self):
        self.candidate('p')
        self.new[0]['features']['seq.max_raw_any'] = .8
        with self.assertRaisesRegex(ValueError, 'changed input features'):
            self.audit()

    def test_order_independent_candidate_comparison(self):
        self.candidate('p')
        self.candidate('q')
        self.new.reverse()
        self.assertEqual(self.audit()['candidates'], 2)

    def test_missing_report_candidate_rejected(self):
        self.candidate('p')
        self.records.clear()
        self.save()
        with self.assertRaisesRegex(ValueError, 'coverage differs'):
            C.audit_genome(self.root/'before.tsv', self.root/'after.tsv', self.root/'old.jsonl',
                           self.root/'new.jsonl', self.model)


class SourceProductionPlansTest(unittest.TestCase):
    def test_all_forty_models_keep_their_own_sources_and_context(self):
        import run_source_arm_context_comparison as R
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as folder:
            root=Path(folder);protocol=dict(genomes=['g'],plans={})
            for source in R.ARMS:
                protocol['plans'][source]={}
                for context in R.CONTEXTS:
                    manifest=root/(source+'-'+context+'.tsv')
                    manifest.write_text('genome_id\tclaims_path\tout_path\ng\t'+source+'.claims\tbefore.tsv\n')
                    command=['java','-jar','engine.jar','integrate']
                    if context!='SEQONLY':command+=['--cohort','--cohort-context-mode',context,'--features-out','old.features']
                    protocol['plans'][source][context]=dict(manifest=str(manifest),command=command,features='original')
            before=copy.deepcopy(protocol)
            with patch.object(C.F.V,'describe_run',return_value={'inputs':{}}) as describe:
                plans,_=C.make_plans(protocol,root/'training',root/'production',[])
            self.assertEqual(protocol,before)
            self.assertEqual(describe.call_count,40)
            models=set()
            for source,contexts in plans.items():
                for context,kinds in contexts.items():
                    for kind,plan in kinds.items():
                        models.add(plan['model'])
                        self.assertEqual(Path(plan['model']),root/'training'/source/context/'models'/f'model_{kind}.json')
                        self.assertEqual('--cohort' in plan['command'],context!='SEQONLY')
                        self.assertEqual('--features-out' in plan['command'],context!='SEQONLY')
                        with Path(plan['manifest']).open() as stream:
                            row=next(csv.DictReader(stream,delimiter='\t'))
                        self.assertEqual(row['claims_path'],source+'.claims')
                        self.assertEqual(Path(row['out_path']).parent,Path(plan['inference']))
            self.assertEqual(len(models),40)


class SourceProductionResumeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home() / 'tmp')))
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.training = self.root / 'training'
        self.out = self.root / 'replay'
        self.training.mkdir()
        self.protocol = dict(genomes=['g'], plans={'source': {}})
        models = []
        for context in ('SEQONLY', 'ORGANISM', 'UNIVERSAL'):
            manifest = self.root / (context + '.tsv')
            manifest.write_text('genome_id\tclaims_path\tout_path\ng\tclaims.jsonl\tbefore.tsv\n')
            self.protocol['plans']['source'][context] = dict(
                manifest=str(manifest), command=['java', '-jar', 'engine.jar', 'integrate'],
                features=str(self.root / 'baseline-features'))
            for kind in C.F.KINDS:
                model = self.training / 'source' / context / 'models' / ('model_' + kind + '.json')
                model.parent.mkdir(parents=True, exist_ok=True)
                model.write_text('{}')
                models.append(model)
        for name in ('protocol.json', 'registration.json', 'evaluation.json'):
            (self.training / name).write_text('{}' if name != 'registration.json' else '{"inputs": {}}')
        (self.training / 'export.json').write_text(json.dumps(dict(status='verified', outputs={})))
        (self.training / 'training.json').write_text(json.dumps(dict(status='verified', outputs=C.F.hashes(models))))

    def invoke(self, resume=False, run_batch=None, arms=None, partial_receipt=None):
        import run_source_arm_context_comparison as R
        run_batch = run_batch or MagicMock()
        with patch.object(R, 'verify', return_value=self.protocol), \
             patch.object(R.C, 'sources', return_value=[]), \
             patch.object(C.F, 'KINDS', ('logistic',)), \
             patch.object(C.F.V, 'describe_run', return_value={'inputs': {}}), \
             patch.object(C.F.V, 'run_batch', run_batch), \
             patch.object(R, 'sequence_stage_audit') as sequence_audit, \
             patch.object(C, 'audit_genome', return_value={'candidates': 0}) as audit:
            result = C.run(self.training, self.out, resume=resume, arms=arms,
                           partial_receipt=partial_receipt)
        return result, run_batch, sequence_audit, audit

    def test_resume_reuses_audits_and_only_reruns_incomplete_arm(self):
        self.invoke()
        reused = self.out / 'source/ORGANISM/logistic'
        audited = self.out / 'source/SEQONLY/logistic'
        rerun = self.out / 'source/UNIVERSAL/logistic'
        (audited / 'audit.json').unlink()
        inference = audited / 'inference'
        for name in ('g.integrated.tsv', 'g.integrated.tsv.learned.json', 'batch_stats.tsv'):
            (inference / name).write_text(name)
        (rerun / 'audit.json').unlink()
        result, run_batch, sequence_audit, audit = self.invoke(resume=True)
        self.assertEqual(run_batch.call_count, 1)
        self.assertEqual(Path(run_batch.call_args.args[1]), rerun / 'manifest.tsv')
        self.assertEqual(audit.call_count, 2)
        self.assertEqual(sequence_audit.call_count, 1)
        self.assertEqual(result['resume'], dict(reused=['source/ORGANISM/logistic'],
            audited=['source/SEQONLY/logistic'], run=['source/UNIVERSAL/logistic'], auditor_change=None))
        saved = json.loads((self.out / 'verification.json').read_text())
        self.assertEqual(saved['resume'], result['resume'])
        self.assertTrue(reused.joinpath('audit.json').is_file())

    def test_partial_resume_processes_only_selected_arms_and_writes_partial_receipt(self):
        self.invoke()
        selected = ['source/ORGANISM/logistic', 'source/SEQONLY/logistic']
        audited = self.out / 'source/SEQONLY/logistic'
        (audited / 'audit.json').unlink()
        for name in ('g.integrated.tsv', 'g.integrated.tsv.learned.json', 'batch_stats.tsv'):
            (audited / 'inference' / name).write_text(name)
        untouched = self.out / 'source/UNIVERSAL/logistic/audit.json'
        before = untouched.read_bytes()
        (self.out / 'verification.json').unlink()
        receipt = self.root / 'partial.json'
        with patch.object(Path, 'rglob', side_effect=AssertionError('partial replay hashed output tree')):
            result, run_batch, sequence_audit, audit = self.invoke(resume=True, arms=selected,
                                                                   partial_receipt=receipt)
        run_batch.assert_not_called()
        self.assertEqual(audit.call_count, 1)
        sequence_audit.assert_called_once()
        self.assertEqual(untouched.read_bytes(), before)
        self.assertFalse((self.out / 'verification.json').exists())
        self.assertEqual(result['schema'], 'gspa.source-context-production-partial/1')
        self.assertEqual(result['status'], 'verified')
        self.assertEqual(result['arms'], selected)
        self.assertEqual(result['dispositions'], dict(reused=['source/ORGANISM/logistic'],
            audited=['source/SEQONLY/logistic'], run=[]))
        self.assertEqual(set(result['summaries']), {'source'})
        self.assertEqual(set(result['summaries']['source']), {'ORGANISM', 'SEQONLY'})
        self.assertEqual(result['summaries']['source']['SEQONLY']['logistic'], {'g': {'candidates': 0}})
        self.assertIsNone(result['auditor_change'])
        registration = json.loads((self.out / 'registration.json').read_text())
        self.assertEqual(result['inputs'], registration['inputs'])
        self.assertEqual(json.loads(receipt.read_text()), result)

    def test_partial_resume_rejects_unknown_duplicate_and_nonresume_arms_before_running(self):
        self.invoke()
        receipt = self.root / 'partial.json'
        for selected, message in [(['source/UNKNOWN/logistic'], 'unknown replay arm'),
                                  (['source/SEQONLY/logistic', 'source/SEQONLY/logistic'], 'duplicate replay arm')]:
            run_batch = MagicMock()
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                self.invoke(resume=True, arms=selected, run_batch=run_batch, partial_receipt=receipt)
            run_batch.assert_not_called()
        run_batch = MagicMock()
        with self.assertRaisesRegex(ValueError, '--arm requires --resume'):
            self.invoke(arms=['source/SEQONLY/logistic'], run_batch=run_batch, partial_receipt=receipt)
        run_batch.assert_not_called()
        self.assertFalse(receipt.exists())

    def test_partial_resume_requires_a_receipt(self):
        with self.assertRaisesRegex(ValueError, '--arm requires --partial-receipt'):
            self.invoke(resume=True, arms=['source/SEQONLY/logistic'])
        self.assertFalse(self.out.exists())

    def test_final_resume_reuses_audits_written_by_partial_runs(self):
        self.invoke()
        for arm in ('SEQONLY', 'ORGANISM'):
            (self.out / 'source' / arm / 'logistic/audit.json').unlink()
        (self.out / 'verification.json').unlink()
        selected = ['source/SEQONLY/logistic', 'source/ORGANISM/logistic']
        partial, run_batch, _, _ = self.invoke(resume=True, arms=selected,
                                               partial_receipt=self.root / 'partial.json')
        self.assertEqual(partial['dispositions']['run'],
                         ['source/ORGANISM/logistic', 'source/SEQONLY/logistic'])
        self.assertEqual(run_batch.call_count, 2)
        result, run_batch, _, audit = self.invoke(resume=True)
        run_batch.assert_not_called()
        audit.assert_not_called()
        self.assertEqual(result['resume'], dict(reused=[
            'source/ORGANISM/logistic', 'source/SEQONLY/logistic', 'source/UNIVERSAL/logistic'],
            audited=[], run=[], auditor_change=None))
        self.assertTrue((self.out / 'verification.json').is_file())

    def test_resume_tolerates_stale_auditor_hash_and_records_change(self):
        self.invoke()
        auditor = str(Path(C.__file__).resolve())
        registered = '0' * 64
        registration = json.loads((self.out / 'registration.json').read_text())
        registration['inputs'][auditor] = registered
        C.F.save(self.out / 'registration.json', registration)
        result, run_batch, _, _ = self.invoke(resume=True)
        run_batch.assert_not_called()
        current = C.S.sha(auditor)
        self.assertEqual(result['resume']['auditor_change'], dict(
            path=auditor, registered_sha256=registered, current_sha256=current))
        self.assertEqual(result['inputs'][auditor], current)
        saved = json.loads((self.out / 'verification.json').read_text())
        self.assertEqual(saved['resume']['auditor_change'], result['resume']['auditor_change'])
        self.assertEqual(saved['inputs'][auditor], current)

    def test_resume_rejects_other_stale_input_before_running_arms(self):
        self.invoke()
        offending = str((self.out / 'source/SEQONLY/logistic/command.json').resolve())
        registration = json.loads((self.out / 'registration.json').read_text())
        registration['inputs'][offending] = '0' * 64
        C.F.save(self.out / 'registration.json', registration)
        run_batch = MagicMock()
        with self.assertRaises(ValueError) as raised:
            self.invoke(resume=True, run_batch=run_batch)
        self.assertIn(offending, str(raised.exception))
        run_batch.assert_not_called()

    def test_resume_rejects_tampered_arm_command(self):
        self.invoke()
        command = self.out / 'source/SEQONLY/logistic/command.json'
        command.write_text('["tampered"]')
        with self.assertRaisesRegex(ValueError, 'inputs or outputs changed'):
            self.invoke(resume=True)

    def test_resume_requires_existing_output_and_fresh_run_refuses_one(self):
        with self.assertRaisesRegex(ValueError, 'resume output directory does not exist'):
            self.invoke(resume=True)
        self.invoke()
        with self.assertRaises(FileExistsError):
            self.invoke()


if __name__ == '__main__':
    unittest.main()
