import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import analyze_source_arm_context_models as C


def row(pid, label, split='calibration', genome='cal'):
    return dict(example_id=pid, protein_id=pid, function_type='EC', function_id='EC:1.1.1.1',
                genome_id=genome, split=split, label=label,
                label_basis='positive_truth' if label else 'closed_world_unconfirmed')


class SourceAnalysisTest(unittest.TestCase):
    def test_no_call_is_not_a_score_of_zero_at_zero_threshold(self):
        attached = C.attach([row('present',1), row('missing',1), row('negative',0)],
                           {'cal':{('present','EC','1.1.1.1'):0.}})
        cut = C.threshold(attached)
        self.assertEqual((cut['value'],cut['missing_predictions']), (0.,2))
        test = [dict(r,split='test') for r in attached]
        result = C.conditional_metrics(test, 0.)
        self.assertEqual((result['reference_matched_calls'],result['missed_reference'],
                          result['unconfirmed_or_negative_calls']), (1,1,0))
        self.assertEqual(result['brier_denominator'],1)
        with self.assertRaisesRegex(ValueError,'calibration rows only'):
            C.threshold(test)

    def test_missing_genome_and_duplicate_supervision_rejected(self):
        with self.assertRaisesRegex(ValueError,'missing genome'):
            C.attach([row('p',1)],{})
        with self.assertRaisesRegex(ValueError,'duplicate supervised'):
            C.attach([row('p',1),row('p',1)],{'cal':{}})

    def test_comparisons_include_context_source_and_native_clean_controls(self):
        pairs = C.comparisons()
        self.assertEqual(len(C.ARMS),61)
        self.assertIn(('clean_diamond_family_excluded/PATHWAY/gbdt','clean/PATHWAY/gbdt'),pairs)
        self.assertIn(('clean/ORGANISM/logistic','clean/SEQONLY/logistic'),pairs)
        self.assertIn(('diamond_homology/SEQONLY/rule','clean_native'),pairs)
        self.assertNotIn(('clean/PATHWAY/logistic','clean/UNIVERSAL/gbdt'),pairs)

    def test_ec_calibration_counts_references_never_proposed(self):
        refs={('g','p'):{'1.1.1.1'},('g','missing'):{'2.2.2.2'}}
        result=C.ec_threshold(refs,{'g':{('p','EC','1.1.1.1'):.8,('p','EC','9.9.9.9'):.9}},{'g'})
        self.assertEqual((result['value'],result['reference_pairs'],result['emitted_reference_pairs']),(.8,2,1))
        self.assertEqual(result['calibration_f1'],.5)

    def test_paired_reference_counts_gains_losses_and_unconfirmed_calls(self):
        refs={('g','a'):{'1.1.1.1'},('g','b'):{'2.2.2.2'}}
        left=C.EC.evaluate(refs,{'g':{('a','EC','1.1.1.1'):.9,('b','EC','3.3.3.3'):.8}},{'g'},.5)
        right=C.EC.evaluate(refs,{'g':{('b','EC','2.2.2.2'):.9}},{'g'},.5)
        result=C.paired_reference(left,right)
        self.assertEqual((result['top1_gained'],result['top1_lost'],result['reference_calls_gained'],
                          result['reference_calls_lost'],result['unconfirmed_call_change']),(1,1,1,1,1))
        right['records'].pop()
        with self.assertRaisesRegex(ValueError,'population differs'):
            C.paired_reference(left,right)

    def test_full_reference_denominator_and_joint_threshold_freeze(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as folder:
            out=Path(folder)
            labels=out/'labels.tsv';splits=out/'splits.tsv';table=out/'table.jsonl'
            labels.write_text('genome_id\tprotein_id\tfunction_id\tevidence\tlabel\n'
                              'cal\tpos\tEC:1.1.1.1\tECO:0000269\t1\n'
                              'test\tp\tEC:1.1.1.1\tECO:0000269\t1\n'
                              'test\tmissing\tEC:2.2.2.2\tECO:0000269\t1\n')
            splits.write_text('genome_id\tsplit\ncal\tcalibration\ntest\ttest\n')
            table.write_text(''.join(json.dumps(r)+'\n' for r in
                                    [row('pos',1),row('neg',0),row('p',1,'test','test')]))
            paths={'cal':'cal.tsv','test':'test.tsv'}
            plans={a:dict(table=str(table),paths=paths,format='claims' if a=='clean_native' else 'integrated') for a in C.ARMS}
            protocol=dict(C.POLICY,plans=plans,comparisons=C.comparisons(),genomes=['cal','test'],
                          references={'labels':str(labels),'splits':str(splits)})
            C.F.save(out/'protocol.json',protocol)
            C.F.save(out/'registration.json',dict(inputs=C.F.hashes([labels,splits,table]),
                                                  protocol_sha256=C.S.sha(out/'protocol.json')))
            def read(path):
                return {('pos','EC','1.1.1.1'):.9,('neg','EC','1.1.1.1'):.1} if path.name=='cal.tsv' else {('p','EC','1.1.1.1'):.7}
            original=C.EC.evaluate; calls=[]
            def evaluate(*args):
                self.assertEqual(set(json.loads((out/'thresholds.json').read_text())),set(C.ARMS))
                calls.append(1)
                return original(*args)
            with patch.object(C.A,'read_scores',side_effect=read),patch.object(C,'raw_clean_scores',side_effect=read), \
                 patch.object(C.EC,'evaluate',side_effect=evaluate):
                result=C.run(out)
                with self.assertRaises(FileExistsError):C.run(out)
            self.assertEqual(len(calls),184)
            for policy,summary in result.items():
                if policy=='native_emitted_calls':
                    self.assertEqual(summary['full_ec_reference']['threshold'],0.)
                    self.assertEqual(summary['full_ec_reference']['reference_matches'],1)
                    continue
                expected=0 if policy in ('calibration_selected','ec_calibration_selected') else 1
                for panel in summary['full_ec_reference'].values():
                    self.assertEqual((panel['proteins'],panel['reference_pairs'],panel['reference_matches']), (2,2,expected))
                    self.assertEqual(panel['emitted_reference_pairs'],1)
                    self.assertEqual(panel['proteins_with_candidates'],1)


if __name__=='__main__':unittest.main()
