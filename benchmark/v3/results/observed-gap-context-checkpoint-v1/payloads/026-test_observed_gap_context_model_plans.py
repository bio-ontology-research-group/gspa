import copy,csv,os
from pathlib import Path
import sys,tempfile,unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import check_observed_gap_context_models as C

class ObservedGapContextModelPlansTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp')))
        self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name);self.protocol=dict(genomes=['g'],plans={})
        for transfer in C.R.G.MODES:
            self.protocol['plans'][transfer]={}
            for mode in C.R.MODES:
                name=transfer+'-'+mode;manifest=self.root/(name+'.tsv')
                manifest.write_text('genome_id\tclaims_path\tgaps_path\tout_path\ng\t'+transfer+'.claims\tactual.gaps\t'+name+'.baseline.tsv\n')
                self.protocol['plans'][transfer][mode]=dict(manifest=str(manifest),features=name+'.features',
                    command=['java','--cohort','--rlc-suggester','--cohort-context-mode',mode])
    def test_all_twelve_plans_keep_context_baseline_and_matching_model(self):
        old=copy.deepcopy(self.protocol);training=self.root/'training';source=self.root/'source.py'
        with patch.object(C.F.V,'describe_run',return_value={'inputs':{}}) as describe:
            plans,_=C.make_plans(self.protocol,training,self.root/'out',[source])
        self.assertEqual(old,self.protocol);self.assertEqual(describe.call_count,12);models=set()
        for transfer,modes in plans.items():
            for mode,kinds in modes.items():
                for kind,p in kinds.items():
                    self.assertEqual(Path(p['model']),training/transfer/mode/'models'/f'model_{kind}.json');models.add(p['model'])
                    self.assertEqual(p['baseline'],self.protocol['plans'][transfer][mode])
                    self.assertEqual(p['command'][p['command'].index('--cohort-context-mode')+1],mode)
                    self.assertIn('--rlc-suggester',p['command'])
                    with Path(p['manifest']).open() as f:row=next(csv.DictReader(f,delimiter='\t'))
                    self.assertEqual(row['claims_path'],transfer+'.claims');self.assertEqual(row['gaps_path'],'actual.gaps')
                    self.assertEqual(Path(row['out_path']).parent,Path(p['inference']))
        self.assertEqual(len(models),12)
        for call in describe.call_args_list:self.assertEqual(call.kwargs['sources'],[source])
        with self.assertRaises(FileExistsError):C.make_plans(self.protocol,training,self.root/'out',[])
    def test_missing_or_duplicate_genomes_fail(self):
        for ids in [('other',),('g','g')]:
            p=Path(self.protocol['plans']['homology']['UNIVERSAL']['manifest'])
            p.write_text('genome_id\tclaims_path\tgaps_path\tout_path\n'+''.join(f'{g}\tc\tgap\tout\n' for g in ids))
            with self.subTest(ids=ids),self.assertRaisesRegex(ValueError,'coverage differs'):
                C.make_plans(self.protocol,self.root/'training',self.root/'out',[])

if __name__=='__main__':unittest.main()
