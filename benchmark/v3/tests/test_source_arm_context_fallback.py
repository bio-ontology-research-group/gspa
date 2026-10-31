"""A complete unsupported attempt is allowed only by an explicit replay policy."""
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import verified_batch as V
import run_source_arm_context_comparison as R


class SourceFallbackTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp')))
        self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.model=self.root/'model.json';self.model.write_text('{}')
        self.out=self.root/'g.integrated.tsv'
        self.out.write_text('protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\n')
        self.manifest=self.root/'manifest.tsv'
        self.manifest.write_text(f'genome_id\tclaims_path\tout_path\ng\t{self.model}\t{self.out}\n')
        self.command=[sys.executable,'--learned-model',str(self.model)]
        self.stage=dict(name='learned_integration',status='skipped',requested=True,
            reason='no candidates supported by the model input-coverage manifest; baseline scores retained',
            candidates_scored=0,candidates_unsupported=0,model_sha256=R.S.sha(self.model))
        self.report=dict(schema='gspa.learned-stage-report/1',model={'sha256':R.S.sha(self.model)},
            counts=dict(candidates=0,unsupported=0),candidates=[],outcomes=[])

    def save(self):
        stages=[dict(name=n,status='executed') for n in ('reference_preparation','evidence_integration','feature_export')]
        stages.append(self.stage)
        Path(str(self.out)+'.stages.json').write_text(json.dumps(dict(schema_version=2,genome_id='g',stages=stages)))
        Path(str(self.out)+'.learned.json').write_text(json.dumps(self.report))

    def test_empty_genome_requires_explicit_policy_and_retains_both_artifacts(self):
        self.save()
        spec=V.describe_run(self.command,self.manifest)
        with self.assertRaisesRegex(ValueError,'did not execute'):V.output_hashes(spec)
        allowed=V.describe_run(self.command,self.manifest,allow_unsupported_learned=True)
        self.assertTrue(allowed['allow_unsupported_learned'])
        outputs=V.output_hashes(allowed)
        self.assertIn(str(self.out),outputs)
        self.assertIn(str(self.out)+'.learned.json',outputs)
        del self.report['candidates'];self.save()
        with self.assertRaisesRegex(ValueError,'did not execute'):V.output_hashes(allowed)

    def test_nonempty_unsupported_and_mutations(self):
        self.stage['candidates_unsupported']=1
        self.report.update(counts=dict(candidates=1,unsupported=1),
            candidates=[dict(function_key='p|EC|EC:1.1.1.1',score_status='unsupported',unsupported_reason='source',
                             probability=None,margin=None,above_display_threshold=None)],
            outcomes=[dict(function_key='p|EC|EC:1.1.1.1',action='unchanged_unsupported',learned_probability=None)])
        self.save();spec=V.describe_run(self.command,self.manifest,allow_unsupported_learned=True)
        V.output_hashes(spec)
        for defect in ('scored','probability','changed','wrong_model','wrong_reason','missing_outcome'):
            report,stage=copy.deepcopy(self.report),copy.deepcopy(self.stage)
            try:
                if defect=='scored':self.stage['candidates_scored']=1
                elif defect=='probability':self.report['candidates'][0]['probability']=.5
                elif defect=='changed':self.report['outcomes'][0]['action']='added'
                elif defect=='wrong_model':self.report['model']['sha256']='0'*64
                elif defect=='wrong_reason':self.stage['reason']='file missing'
                else:self.report['outcomes']=[]
                self.save()
                with self.subTest(defect=defect),self.assertRaises(ValueError):V.output_hashes(spec)
            finally:self.report,self.stage=report,stage

    def test_sequence_only_rejects_an_executed_context_stage(self):
        self.save()
        self.assertEqual(R.sequence_stage_audit(self.root,['g'])['g'],
                         ['evidence_integration','feature_export','reference_preparation'])
        self.stage.update(name='cohort_context',status='executed')
        self.save()
        with self.assertRaisesRegex(ValueError,'executed a context stage'):
            R.sequence_stage_audit(self.root,['g'])

    def test_large_learned_report_has_same_checks_without_loading_candidate_list(self):
        from unittest.mock import patch
        try:import ijson
        except ImportError:self.skipTest('optional streaming parser unavailable')
        self.stage['candidates_unsupported']=1
        self.report.update(counts=dict(candidates=1,unsupported=1),
            candidates=[dict(function_key='p|EC|EC:1.1.1.1',score_status='unsupported',unsupported_reason='source',
                             probability=None,margin=None,above_display_threshold=None,supporting_claims=[{'large':'ignored'}])],
            outcomes=[dict(function_key='p|EC|EC:1.1.1.1',action='unchanged_unsupported',learned_probability=None)])
        self.save();path=Path(str(self.out)+'.learned.json')
        small=V.learned_report_summary(path)
        with patch.object(V,'LEDGER_MEMORY_LIMIT',1):
            large=V.learned_report_summary(path)
            # Large mode retains only the model hash, rather than optional header metadata.
            self.assertEqual(small,large)
            spec=V.describe_run(self.command,self.manifest,allow_unsupported_learned=True)
            V.output_hashes(spec)
            with path.open('a') as stream:stream.write(' broken')
            with self.assertRaises((ValueError,ijson.JSONError)):V.learned_report_summary(path)


if __name__=='__main__':unittest.main()
