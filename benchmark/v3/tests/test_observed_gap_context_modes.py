import copy,csv,json,os
from contextlib import ExitStack
from pathlib import Path
import sys,tempfile,unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import run_observed_gap_context_modes as C

class ObservedGapContextModesTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp')))
        self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.protocol=dict(genomes=['g'],gaps={'g':{'path':'actual.gaps'}},plans={})
        for transfer in C.G.MODES:
            manifest=self.root/(transfer+'.tsv')
            manifest.write_text('genome_id\tclaims_path\tgaps_path\tout_path\ng\t'+transfer+'.claims\tactual.gaps\toriginal.tsv\n')
            self.protocol['plans'][transfer]={'observed_gaps':dict(manifest=str(manifest),command=['java','-jar','old.jar','--cohort','--rlc-suggester'],features='original.features',inference='original.inference')}
    def test_all_eight_plans_preserve_claims_and_gaps_with_matched_mode_inputs(self):
        original=copy.deepcopy(self.protocol);context=self.root/'context';out=self.root/'out';source=self.root/'source.py'
        with patch.object(C.F.V,'describe_run',return_value={'inputs':{}}) as describe:
            plans,_=C.make_plans(self.protocol,context,self.root/'current.jar',out,[source])
        self.assertEqual(self.protocol,original);self.assertEqual(describe.call_count,8)
        for transfer in C.G.MODES:
            self.assertEqual(list(plans[transfer]),list(C.CONTROLS))
            for mode,plan in plans[transfer].items():
                command=plan['command'];self.assertIn('--rlc-suggester',command);self.assertNotIn('--learned-model',command)
                with Path(plan['manifest']).open() as f:row=next(csv.DictReader(f,delimiter='\t'))
                self.assertEqual(row['claims_path'],transfer+'.claims');self.assertEqual(row['gaps_path'],'actual.gaps')
                if mode=='default':self.assertNotIn('--cohort-context-mode',command)
                else:
                    self.assertEqual(command[command.index('--cohort-context-mode')+1],mode)
                    self.assertEqual(command[command.index('--cohort-pathways')+1],str(context/'pathway-groups.tsv'))
                    self.assertEqual(command[command.index('--cohort-metadata')+1],str(context/'organism-metadata.json'))
        for call in describe.call_args_list:self.assertEqual(call.kwargs['sources'],[source])
        with self.assertRaises(FileExistsError):C.make_plans(self.protocol,context,self.root/'j.jar',out,[])
    def test_invalid_coverage_gap_or_scored_baseline_cannot_register(self):
        for defect in ['coverage','gap','model']:
            p=copy.deepcopy(self.protocol);baseline=p['plans']['homology']['observed_gaps']
            if defect=='coverage':p['genomes']=['other']
            elif defect=='gap':p['gaps']['g']['path']='different.gaps'
            else:baseline['command']+=['--learned-model','model.json']
            with self.subTest(defect=defect),self.assertRaises(ValueError):
                C.make_plans(p,self.root/'context',self.root/'j.jar',self.root/defect,[])
    def pair(self,old,new):
        a=self.root/'before.jsonl';b=self.root/'after.jsonl'
        a.write_text(json.dumps(old)+'\n');b.write_text(json.dumps(new)+'\n');return a,b
    def test_context_comparison_preserves_non_context_features_and_provenance(self):
        old=dict(candidate_key='p|EC|1|CATALYTIC',genome_id='g',sequence_sources=['s'],sequence_contracts=[{'source':'s'}],
            quality_evidence=[],cohort_context_mode='UNIVERSAL',features={'cohort.w.ec_direct':1.,'cohort.w.local':.2,'ctx.gap.available':1.,'seq.x':.3})
        new=copy.deepcopy(old);new['cohort_context_mode']='PATHWAY';new['features']['cohort.w.local']=.8
        self.assertEqual(C.compare_features(*self.pair(old,new),'PATHWAY')['changed_columns'],{'cohort.w.local':1})
        with self.assertRaises(ValueError):C.compare_features(*self.pair(old,new),'PATHWAY',exact=True)
        for defect in ('sequence_contracts','quality_evidence','gap','sequence','mode','nonfinite'):
            bad=copy.deepcopy(new)
            if defect in ('sequence_contracts','quality_evidence'):bad[defect]=['changed']
            elif defect=='gap':bad['features']['ctx.gap.available']=0
            elif defect=='sequence':bad['features']['seq.x']=.5
            elif defect=='nonfinite':bad['features']['cohort.w.local']=float('nan')
            else:bad['cohort_context_mode']='ORGANISM'
            with self.subTest(defect=defect),self.assertRaises(ValueError):C.compare_features(*self.pair(old,bad),'PATHWAY')
    def test_supervision_allows_context_changes_only(self):
        old=dict(example_id='e',label=1,split='test',role='CATALYTIC',cohort_context_mode='UNIVERSAL',features={'cohort.w.ec_direct':1})
        new=copy.deepcopy(old);new['cohort_context_mode']='ORGANISM'
        self.assertEqual(C.compare_supervision(*self.pair(old,new),'ORGANISM'),1)
        new['label']=0
        with self.assertRaisesRegex(ValueError,'supervision identity'):C.compare_supervision(*self.pair(old,new),'ORGANISM')
    def test_all_supervision_dependencies_are_bound_and_checked(self):
        paths={}
        for key in ['labels','splits','label_policy','registry','go_obo','example_clusters','locus_groups']:
            path=self.root/key;path.write_text('frozen '+key);paths[key]=str(path)
        recipes={t:dict(supervision=paths) for t in C.G.MODES}
        reference=C.F.hashes([Path(p) for p in paths.values()])
        self.assertEqual(C.supervision_inputs(recipes,reference),reference)
        for key,path in paths.items():
            with self.subTest(key=key):
                incomplete={p:h for p,h in reference.items() if p!=path}
                with self.assertRaisesRegex(ValueError,'unregistered supervision'):
                    C.supervision_inputs(recipes,incomplete)
                original=Path(path).read_bytes();Path(path).write_text('changed')
                with self.assertRaises(ValueError):C.supervision_inputs(recipes,reference)
                Path(path).write_bytes(original)
        Path(paths['registry']).unlink()
        with self.assertRaises((ValueError,FileNotFoundError)):C.supervision_inputs(recipes,reference)
    def test_all_twelve_models_freeze_before_test_evaluation(self):
        out=self.root/'fit';out.mkdir();(out/'export.json').write_text(json.dumps(dict(status='verified',outputs={})))
        recipe=dict(supervision={k:'input' for k in ['labels','splits','label_policy','registry','go_obo','example_clusters']},
            original_table='original',disjoint_keys=['genome_id'],settings={k:{} for k in C.F.KINDS})
        p=dict(recipes={t:copy.deepcopy(recipe) for t in C.G.MODES})
        for t in C.G.MODES:
            for m in C.MODES:
                f=out/t/m/'features/g';f.mkdir(parents=True);(f/'feature_contract.json').write_text('{}')
        rows=[{'split':s} for s in C.F.T.SPLITS]
        def build(*args,**kwargs):
            dest=args[5];dest.mkdir();(dest/'training_table.jsonl').write_text('{}\n');(dest/'training_table_manifest.json').write_text('{}')
            return {'table':str(dest/'training_table.jsonl')}
        def fit(train,calibration,*args,**kwargs):
            self.assertTrue(all(r['split']=='train' for r in train));self.assertTrue(all(r['split']=='calibration' for r in calibration))
            self.assertFalse((out/'training.json').exists());return {'stub':True}
        def evaluate(model,contract,test,*args):
            self.assertEqual(len(list(out.glob('*/*/models/model_*.json'))),12)
            self.assertEqual(json.loads((out/'training.json').read_text())['status'],'verified')
            self.assertEqual(test,[{'split':'test'}]);return {'n':1}
        with ExitStack() as stack:
            for obj,name,value in [(C,'verify',p),(C,'compare_supervision',3),(C.F.T,'load_contract',{}),
                (C.F.T,'read_table',rows),(C.F.T,'check_disjoint',None),(C.F.T,'split_identity',{})]:
                stack.enter_context(patch.object(obj,name,return_value=value))
            stack.enter_context(patch.object(C.F.P,'build_table',side_effect=build))
            train=stack.enter_context(patch.object(C.F.T,'train',side_effect=fit))
            stack.enter_context(patch.object(C.F.T,'write_frozen',side_effect=lambda a,p:p.write_text(json.dumps(a))))
            test=stack.enter_context(patch.object(C.F.T,'evaluate_frozen',side_effect=evaluate))
            C.train(out);self.assertEqual(train.call_count,12);test.assert_not_called()
            C.evaluate(out);self.assertEqual(test.call_count,12)
            with self.assertRaises(FileExistsError):C.evaluate(out)
    def test_incomplete_freeze_rejected_before_test_rows_are_read(self):
        out=self.root/'incomplete';out.mkdir();(out/'training.json').write_text(json.dumps({'status':'verified','outputs':{}}))
        with patch.object(C,'verify',return_value={}),patch.object(C.F.T,'read_table') as read:
            with self.assertRaisesRegex(ValueError,'missing frozen model'):C.evaluate(out)
            read.assert_not_called()

if __name__=='__main__':unittest.main()
