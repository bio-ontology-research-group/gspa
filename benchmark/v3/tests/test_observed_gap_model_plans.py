import copy
import csv
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import check_observed_gap_models as C


class ObservedGapModelPlanTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp')))
        self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.protocol=dict(genomes=['g'],plans={})
        for mode in C.G.MODES:
            p=self.root/(mode+'.tsv')
            p.write_text('genome_id\tclaims_path\tgaps_path\tout_path\n'+f'g\t{mode}.claims\tobserved.gaps\t{mode}.baseline.tsv\n')
            command=['java','--cohort','--rlc-suggester','--batch-manifest',str(p)]
            self.protocol['plans'][mode]=dict(observed_gaps=dict(manifest=str(p),command=command,features=mode+'.features'),
                no_gaps=dict(manifest='must-not-use',command=['must-not-use']))

    def test_all_four_commands_preserve_gap_inputs_and_use_matching_frozen_models(self):
        original=copy.deepcopy(self.protocol);out=self.root/'out';training=self.root/'training';source=self.root/'source.py'
        with patch.object(C.F.V,'describe_run',return_value={'inputs':{}}) as describe:
            plans,_=C.prepare_plans(self.protocol,training,out,[source])
        self.assertEqual(self.protocol,original);self.assertEqual(describe.call_count,4)
        models=[]
        for mode in C.G.MODES:
            for kind in C.F.KINDS:
                plan=plans[mode][kind];models.append(plan['model'])
                self.assertEqual(Path(plan['model']),training/mode/'default/models'/f'model_{kind}.json')
                self.assertEqual(plan['command'][plan['command'].index('--learned-model')+1],plan['model'])
                self.assertIn('--rlc-suggester',plan['command'])
                with Path(plan['manifest']).open() as f:row=next(csv.DictReader(f,delimiter='\t'))
                self.assertEqual(row['gaps_path'],'observed.gaps')
                self.assertEqual(row['claims_path'],mode+'.claims')
                self.assertEqual(Path(row['out_path']).parent,Path(plan['inference']))
                self.assertEqual(plan['original_rows'][0]['out_path'],mode+'.baseline.tsv')
        self.assertEqual(len(set(models)),4)
        for call in describe.call_args_list:self.assertEqual(call.kwargs['sources'],[source])
        with patch.object(C.F.V,'describe_run',return_value={'inputs':{}}):
            with self.assertRaises(FileExistsError):C.prepare_plans(self.protocol,training,out,[source])

    def test_incomplete_or_duplicate_baseline_genomes_reject_registration(self):
        for ids in (['other'],['g','g']):
            p=Path(self.protocol['plans']['homology']['observed_gaps']['manifest'])
            p.write_text('genome_id\tclaims_path\tgaps_path\tout_path\n'+''.join(f'{g}\tx\ty\tz\n' for g in ids))
            with self.subTest(ids=ids),self.assertRaisesRegex(ValueError,'coverage differs'):
                C.prepare_plans(self.protocol,self.root/'training',self.root/'out',[])


if __name__=='__main__':unittest.main()
