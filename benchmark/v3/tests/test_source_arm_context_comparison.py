import copy,csv,json,os
from contextlib import ExitStack
from pathlib import Path
import sys,tempfile,unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import run_source_arm_context_comparison as R
A=R.A;S=R.S;F=R.F;G=R.G
SMALL=dict(genomes=1,proteins=3,native=2,unsupported=1,exact_known=2)


def jsonl(path,rows):Path(path).write_text(''.join(json.dumps(r)+'\n' for r in rows))


class Base(unittest.TestCase):
    def setUp(self):
        tmp=Path(os.environ.get('TMPDIR',str(Path.home()/'tmp')));tmp.mkdir(exist_ok=True)
        self.tmp=tempfile.TemporaryDirectory(dir=tmp);self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)


class CleanAccountingTest(Base):
    def setUp(self):
        super().setUp();root=self.root/'clean';self.clean=root
        self.contract=dict(schema='gspa.predictor-contract/1',tool='CLEAN',tool_version='fixture',tool_sha256='c'*64,
            reference_sha256='c'*64,parameters_sha256='c'*64,adapter_sha256='c'*64,score_model='fixture')
        genome=dict(genome_id='g',proteins=3,sha256='d'*64,proteome='g.faa')
        self.unsupported=dict(protein_id='u',status='unsupported_input',prediction='abstain',reason='residue_not_in_native_esm_alphabet')
        folder=root/'attempts/g/job-1';folder.mkdir(parents=True)
        self.claims=[dict(source='clean',protein_id=p,function_type='EC',function_id='EC:1.1.1.1',raw_score=.5,
            evidence_type='SEQUENCE_DEEPLEARNING',metadata=dict(predictor_contract=self.contract,genome_id='g')) for p in 'pq']
        jsonl(folder/'claims.jsonl',self.claims);(folder/'completed.json').write_text('{}')
        overlap=dict(rows=[],by_genome={'g':dict(exact_proteins=2,truncated_input_matches=2)},exact_proteins=2,truncated_input_matches=2,interpretation='fixture')
        (root/'reference-overlap.json').write_text(json.dumps(overlap))
        self.plan=dict(schema='gspa.clean-cohort/2',genomes=[genome],unsupported_inputs={'g':[self.unsupported]},predictor_contract=self.contract,
            reference_overlap_sha256=S.sha(root/'reference-overlap.json'))
        (root/'plan.json').write_text(json.dumps(self.plan))
        self.report=dict(schema='gspa.clean-cohort-audit/2',plan_sha256=S.sha(root/'plan.json'),auditor_sha256=S.sha(A.__file__),genomes=1,
            proteins=3,original_genomes=[genome],native_proteins=2,unsupported_proteins=1,
            unsupported_outcomes=[dict(genome_id='g',original_fasta_sha256='d'*64,original_proteins=3,**self.unsupported)],calls=2,
            panels=[dict(genome_id='g',proteins=2,original_proteins=3,attempt=str(folder),claims_sha256=S.sha(folder/'claims.jsonl'),
                completed_sha256=S.sha(folder/'completed.json'))],failed_attempts=[],reference_overlap=overlap)
        self.write()

    def write(self):
        self.audit=self.root/'audit.json';self.audit.write_text(json.dumps(self.report))

    def verify(self):return R.verify_clean(self.clean,self.audit,['g'],SMALL)

    def test_complete_audit_is_bound_without_rerunning_it(self):
        with patch.object(A,'audit',side_effect=AssertionError('audit must not rerun')):r=self.verify()
        s=r['summary'];self.assertEqual((s['proteins'],s['native_proteins'],s['unsupported_proteins']),(3,2,1))
        self.assertEqual(s['unsupported_outcomes'][0]['prediction'],'abstain');self.assertEqual(s['reference_overlap']['exact_proteins'],2)
        self.assertNotIn('rows',s['reference_overlap']);self.assertEqual(s['original_genomes'],self.plan['genomes'])
        self.assertIn(str((self.clean/'plan.json').resolve()),r['inputs'])

    def test_every_defect_fails_closed(self):
        def tamper_claims(t):(self.clean/'attempts/g/job-1/claims.jsonl').write_text('{}\n')
        def invented(t):
            rows=self.claims+[dict(self.claims[0],protein_id='u')];jsonl(self.clean/'attempts/g/job-1/claims.jsonl',rows)
            t.report['panels'][0]['claims_sha256']=S.sha(self.clean/'attempts/g/job-1/claims.jsonl')
        cases={'schema':lambda t:t.report.update(schema='gspa.clean-cohort-audit/1'),
            'plan':lambda t:t.report.update(plan_sha256='0'*64),'auditor':lambda t:t.report.update(auditor_sha256='0'*64),
            'missing_genome':lambda t:t.report.update(panels=[],genomes=0),'native':lambda t:t.report.update(native_proteins=3),
            'dropped_outcome':lambda t:t.report.update(unsupported_outcomes=[]),'overlap':lambda t:t.report['reference_overlap'].update(exact_proteins=1),
            'claims_hash':tamper_claims,'invented_negative':invented,
            'genomes':lambda t:t.report.update(original_genomes=[dict(t.plan['genomes'][0],genome_id='h')])}
        for name,change in cases.items():
            with self.subTest(name),self.assertRaises((ValueError,KeyError)):
                saved=copy.deepcopy(self.report);files={p:p.read_bytes() for p in self.clean.rglob('*') if p.is_file()}
                try:change(self);self.write();self.verify()
                finally:
                    self.report=saved
                    for p,b in files.items():p.write_bytes(b)
                    self.write()
        self.verify()


class PlanTest(Base):
    def setUp(self):
        super().setUp()
        def contract(c):return dict(schema='gspa.predictor-contract/1',tool='diamond',tool_version='f',tool_sha256=c*64,
            reference_sha256=c*64,parameters_sha256=c*64,adapter_sha256=c*64,score_model='fixture')
        self.contracts={'homology':contract('1'),'family_excluded':contract('2')}
        self.cleanc=dict(contract('3'),tool='CLEAN')
        self.base=dict(genomes=['g'],gaps={'g':{'path':'actual.gaps'}},plans={},claims={},contracts={})
        for t in G.MODES:
            claim=self.root/(t+'.claims');jsonl(claim,[dict(source='swissprot_diamond_baseline',protein_id='p',function_type='EC',
                function_id='EC:1.1.1.1',raw_score=.4,evidence_type='SEQUENCE_SIMILARITY',metadata=dict(predictor_contract=self.contracts[t]))])
            self.base['claims'][t]={'g':dict(path=str(claim),sha256=S.sha(claim))}
            self.base['contracts'][t]=[dict(source='swissprot_diamond_baseline',contract=self.contracts[t])]
            manifest=self.root/(t+'.tsv');manifest.write_text('\t'.join(R.MANIFEST_COLUMNS[:7])+'\ng\told.claims\told.out\tbacteria\tbacteria\toperons\tactual.gaps\n')
            command=['java','-Xmx1g','-jar','old.jar','integrate','--batch-manifest','old.tsv','--go-owl','go-'+t+'.owl','--ec2go','e.txt',
                '--pathways','p.tsv','--cohort','--cohort-metadata','m.json','--cohort-orthogroups','og.tsv','--cohort-holdout','CLADE',
                '--rlc-suggester','--rlc-alpha','0.5','--operons','op.tsv','--gaps','x','--genome-layout','l','--orthogroups','o',
                '--cluster-consensus','c','--enable-fba-refiner','--fba-medium','m','--experiment-manifest','e.json','--community-manifest','cm.json',
                '--external-context','x.json','--enable-priors','a,b','--enforce-stage2','--infer-taxon','--features-out','f','--suggestions-out','s',
                '--feature-export-dir','fe']
            self.base['plans'][t]={'observed_gaps':dict(manifest=str(manifest),command=command)}
        cclaim=self.root/'clean.claims';jsonl(cclaim,[dict(source='clean',protein_id='p',function_type='EC',function_id='EC:1.1.1.1',raw_score=.9,
            evidence_type='SEQUENCE_DEEPLEARNING',metadata=dict(predictor_contract=self.cleanc,genome_id='g'))])
        self.clean=dict(allowed=dict(source='clean',contract=self.cleanc),claims={'g':dict(path=str(cclaim),sha256=S.sha(cclaim))})

    def test_five_arms_use_declared_sources_contracts_and_recipes(self):
        out=self.root/'out';wanted={'diamond_homology':{'swissprot_diamond_baseline'},'diamond_family_excluded':{'swissprot_diamond_baseline'},
            'clean':{'clean'},'clean_diamond_homology':{'clean','swissprot_diamond_baseline'},'clean_diamond_family_excluded':{'clean','swissprot_diamond_baseline'}}
        self.assertEqual(set(wanted),set(R.ARMS))
        for arm,sources in wanted.items():
            files,allowed=R.arm_claims(arm,self.base,self.clean,out,['g']);spec=R.ARM_SOURCES[arm]
            self.assertEqual(set(files['g']['by_source']),sources)
            parts=([Path(self.base['claims'][spec['diamond']]['g']['path'])] if spec['diamond'] else [])+([Path(self.clean['claims']['g']['path'])] if spec['clean'] else [])
            self.assertEqual(Path(files['g']['path']).read_bytes(),b''.join(p.read_bytes() for p in parts))
            self.assertEqual([e['contract'] for e in allowed if e['source']!='clean'],[self.contracts[spec['diamond']]] if spec['diamond'] else [])
            self.assertEqual(S.audit_export.__module__,'cohort_sequence_inputs')
        self.assertEqual(R.ARM_SOURCES['clean']['recipe'],'homology');self.assertEqual(R.ARM_SOURCES['clean_diamond_family_excluded']['recipe'],'family_excluded')

    def plans(self,out):
        plans={}
        with patch.object(F.V,'describe_run',return_value={'inputs':{}}):
            for arm in R.ARMS:
                claims,_=R.arm_claims(arm,self.base,self.clean,out,['g']);plans[arm]={}
                for ctx in R.CONTEXTS:plans[arm][ctx],_=R.make_plan(arm,ctx,self.base,claims,self.root/'new.jar',self.root/'context',out,[])
        return plans

    def table(self,plan):
        with Path(plan['manifest']).open() as f:reader=csv.reader(f,delimiter='\t');return list(reader)

    def test_twenty_plans_sequence_only_removes_every_context_route(self):
        original=copy.deepcopy(self.base);plans=self.plans(self.root/'out');self.assertEqual(self.base,original)
        self.assertEqual(sum(len(v) for v in plans.values()),20)
        for arm,by in plans.items():
            recipe=R.ARM_SOURCES[arm]['recipe']
            for ctx,plan in by.items():
                command=plan['command'];value=lambda flag:command[command.index(flag)+1]
                self.assertEqual(value('--go-owl'),'go-'+recipe+'.owl');self.assertEqual(value('-jar'),str(self.root/'new.jar'))
                self.assertEqual(value('--batch-manifest'),plan['manifest']);header,row=self.table(plan)
                self.assertTrue(row[1].endswith('/claims/'+arm+'/g.claims.jsonl'));self.assertEqual(header[:7],R.MANIFEST_COLUMNS[:7])
                if ctx=='SEQONLY':
                    self.assertEqual(header,R.MANIFEST_COLUMNS)
                    for c in R.MANIFEST_COLUMNS:
                        got=row[header.index(c)]
                        if c in R.SEQONLY_BLANK:self.assertEqual(got,'-',c)
                    self.assertEqual((row[3],row[4]),('bacteria','bacteria'))
                    self.assertEqual((row[header.index('enforce_stage2')],row[header.index('enable_priors')]),('false','none'))
                    self.assertEqual(value('--enable-priors'),'none')
                    for t in command:
                        self.assertFalse(t.startswith(R.CONTEXT_PREFIXES),t)
                        self.assertNotIn(t,R.CONTEXT_FLAGS-{'--enable-priors'})
                    for kept in ('--go-owl','--ec2go','--pathways','--feature-export-dir','--batch-stats-out','--batch-manifest'):self.assertIn(kept,command)
                    self.assertEqual(command.count('--enable-priors'),1)
                else:
                    self.assertEqual(value('--cohort-context-mode'),ctx);self.assertIn('--cohort',command);self.assertIn('--rlc-suggester',command)
                    self.assertEqual(row[6],'actual.gaps');self.assertEqual(value('--cohort-pathways'),str(self.root/'context/pathway-groups.tsv'))

    def test_strip_keeps_unrelated_values_and_equals_forms(self):
        command=['integrate','--cohort-holdout=CLADE','--cohort','--go-owl','go.owl','--no-enforce-stage2','--fba-gram','neg','--taxon','2']
        self.assertEqual(R.strip_context(command),['integrate','--go-owl','go.owl','--taxon','2','--enable-priors','none'])

    def test_input_mutations_and_missing_genomes_cannot_register(self):
        out=self.root/'out';plans=self.plans(out)
        with self.assertRaises(FileExistsError),patch.object(F.V,'describe_run',return_value={'inputs':{}}):
            R.arm_claims('clean',self.base,self.clean,out,['g'])
        for defect in ('model','gap','genome','columns','claims'):
            b=copy.deepcopy(self.base);cl=copy.deepcopy(self.clean);baseline=b['plans']['homology']['observed_gaps']
            files={p:p.read_bytes() for p in self.root.glob('*') if p.is_file()}
            if defect=='model':baseline['command']+=['--learned-model','m.json']
            elif defect=='gap':b['gaps']['g']['path']='other.gaps'
            elif defect=='genome':b['genomes']=['g','h']
            elif defect=='columns':
                Path(baseline['manifest']).write_text('claims_path\tgenome_id\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\n'
                    'c\tg\to\tb\tb\t-\tactual.gaps\n')
            else:Path(cl['claims']['g']['path']).write_text('changed\n')
            try:
                with self.subTest(defect),patch.object(F.V,'describe_run',return_value={'inputs':{}}),self.assertRaises((ValueError,KeyError)) as caught:
                    claims,_=R.arm_claims('clean_diamond_homology',b,cl,self.root/('bad-'+defect),b['genomes'])
                    R.make_plan('clean_diamond_homology','SEQONLY',b,claims,self.root/'j.jar',self.root/'c',self.root/('bad-'+defect),[])
                if defect=='columns':self.assertIn('positional',str(caught.exception))
            finally:
                for p,data in files.items():p.write_bytes(data)

    def test_supervision_agreement_requires_shared_truth_sources(self):
        paths={}
        for k in ('labels','label_policy','go_obo','splits'):
            for t in G.MODES:
                p=self.root/f'{k}-{t}';p.write_text('shared' if k!='splits' else t);paths[k,t]=str(p)
        recipes={t:dict(supervision={k:paths[k,t] for k in ('labels','label_policy','go_obo','splits')}) for t in G.MODES}
        with self.assertRaisesRegex(ValueError,'splits'):R.supervision_agreement(recipes)
        Path(paths['splits','family_excluded']).write_bytes(Path(paths['splits','homology']).read_bytes())
        agree=R.supervision_agreement(recipes);self.assertEqual((agree['labels'],agree['splits']),(True,True))
        Path(paths['labels','family_excluded']).write_text('different')
        with self.assertRaisesRegex(ValueError,'labels'):R.supervision_agreement(recipes)


class ExportAuditTest(Base):
    def row(self,key,**features):
        return dict(candidate_key=key,genome_id='g',protein_id='p',function_type='EC',function_id=key,role='CATALYTIC',sequence_sources=['clean'],
            sequence_contracts=[],quality_evidence=[],cohort_context_mode='UNIVERSAL',
            features={'seq.x':.5,'ctx.gap.available':0,'cohort.w.ec_direct':0,'asm.n':0,'fba.n':0,**features})

    def pair(self,old,new):
        a=self.root/'a.jsonl';b=self.root/'b.jsonl';jsonl(a,old);jsonl(b,new);return a,b

    def test_comparison_uses_candidate_keys_not_line_order_and_fails_closed(self):
        rows=[self.row('a'),self.row('b')];moved=copy.deepcopy(rows[::-1]);moved[0]['features']['cohort.w.ec_direct']=1
        moved[0]['cohort_context_mode']='PATHWAY';moved[0]['features']['cohort.w.ec_direct']=1
        r=R.compare_candidates(*self.pair(rows,moved),'PATHWAY',lambda n:n.startswith('cohort.'));self.assertEqual(r['candidates'],2)
        for name,change in dict(extra=lambda t:t.append(self.row('c')),seq=lambda t:t[0]['features'].update({'seq.x':.6}),
            provenance=lambda t:t[0].update(sequence_sources=['x']),mode=lambda t:t[1].update(cohort_context_mode='ORGANISM'),
            nonfinite=lambda t:t[0]['features'].update({'cohort.w.ec_direct':float('nan')})).items():
            bad=copy.deepcopy(moved);change(bad)
            with self.subTest(name),self.assertRaises(ValueError):R.compare_candidates(*self.pair(rows,bad),'PATHWAY',lambda n:n.startswith('cohort.'))

    def test_zero_audit_requires_every_context_column_zero(self):
        d=self.root/'features/g';d.mkdir(parents=True);jsonl(d/'candidate_features.jsonl',[self.row('a')])
        self.assertEqual(R.zero_audit(self.root/'features'),4)
        for change in (dict(features={'ctx.gap.available':1}),dict(features={'cohort.x':1e-9}),dict(cohort_context_mode='PATHWAY')):
            row=self.row('a');row['features'].update(change.get('features',{}));row.update({k:v for k,v in change.items() if k!='features'})
            jsonl(d/'candidate_features.jsonl',[row])
            with self.assertRaises(ValueError):R.zero_audit(self.root/'features')
        jsonl(d/'candidate_features.jsonl',[dict(self.row('a'),features={'seq.x':1})])
        with self.assertRaisesRegex(ValueError,'no context columns'):R.zero_audit(self.root/'features')

    def test_supervision_comparison_preserves_shared_labels_and_splits(self):
        def ex(i,g,split,label,mode='UNIVERSAL'):
            key=i+'|EC|EC:1.1.1.1|CATALYTIC'
            return dict(example_id=g+'|'+key,candidate_key=key,protein_id=i,function_type='EC',function_id='EC:1.1.1.1',role='CATALYTIC',genome_id=g,group='x',split=split,label=label,label_basis='b',dataset='d',cohort_context_mode=mode,features={'cohort.w.ec_direct':0})
        old=[ex('1','g','train',1),ex('2','h','test',0),ex('3','h','test',1)]
        a=self.root/'o.jsonl';b=self.root/'n.jsonl';jsonl(a,old)
        jsonl(b,[ex('1','g','train',1),ex('2','h','test',0),ex('9','h','test',0)])
        self.assertEqual(R.compare_supervision(a,b,'SEQONLY'),dict(shared_examples=2,original_only=1,new_only=1,representation_changes=0))
        jsonl(b,[ex('1','g','train',1)])
        self.assertEqual(R.compare_supervision(a,b,'SEQONLY'),dict(shared_examples=1,original_only=2,new_only=0,representation_changes=0))
        for bad in ([ex('1','g','train',0),ex('2','h','test',0)],[ex('1','g','test',1),ex('2','h','test',0)],
                    [ex('1','g','train',1),ex('2','h','test',0,'PATHWAY')],[ex('9','g','test',1),ex('2','h','test',0)]):
            jsonl(b,bad)
            with self.subTest(bad),self.assertRaises(ValueError):R.compare_supervision(a,b,'SEQONLY')

    def test_new_source_can_label_a_genome_absent_from_old_candidate_table(self):
        a=self.root/'old.jsonl';b=self.root/'new.jsonl';splits=self.root/'splits.tsv'
        def ex(i,g):
            key=i+'|EC|EC:1.1.1.1|CATALYTIC'
            return dict(example_id=g+'|'+key,candidate_key=key,protein_id=i,function_type='EC',function_id='EC:1.1.1.1',role='CATALYTIC',genome_id=g,group='f',split='train',label=1,label_basis='positive_truth',
                        dataset='d',cohort_context_mode='UNIVERSAL',features={'cohort.w.ec_direct':0})
        jsonl(a,[ex('one','old')]);jsonl(b,[ex('two','new')])
        splits.write_text('genome_id\tsplit\nold\ttrain\nnew\ttrain\n')
        self.assertEqual(R.compare_supervision(a,b,'SEQONLY',splits),dict(shared_examples=0,original_only=1,new_only=1,representation_changes=0))
        splits.write_text('genome_id\tsplit\nold\ttrain\nnew\ttest\n')
        with self.assertRaisesRegex(ValueError,'genome split changed'):R.compare_supervision(a,b,'SEQONLY',splits)


class SupervisionIdentityMigrationTest(Base):
    def test_prefix_migration_preserves_truth_and_detects_collisions(self):
        key='custom|p|EC|1.1.1.1|CATALYTIC'
        old=dict(example_id='g|'+key,candidate_key=key,genome_id='g',protein_id='custom|p',function_type='EC',
                 function_id='EC:1.1.1.1',role='CATALYTIC',group='family',split='train',label=1,
                 label_basis='positive_truth',label_references=['PMID:123'],dataset='d',
                 cohort_context_mode='UNIVERSAL',features={'cohort.w.ec_direct':0})
        new=dict(old,candidate_key='custom|p|EC|EC:1.1.1.1|CATALYTIC',example_id='g|custom|p|EC|EC:1.1.1.1|CATALYTIC')
        a=self.root/'old.jsonl';b=self.root/'new.jsonl';jsonl(a,[old]);jsonl(b,[new])
        self.assertEqual(R.compare_supervision(a,b,'SEQONLY'),
                         dict(shared_examples=1,original_only=0,new_only=0,representation_changes=1))
        jsonl(a,[old,new])
        with self.assertRaisesRegex(ValueError,'duplicate canonical supervised example'):
            R.compare_supervision(a,b,'SEQONLY')
        jsonl(a,[old])
        for bad in ([old,new],[dict(new,label=0)],[dict(new,label_references=['PMID:456'])],
                    [dict(new,role='ACCESSORY')],[dict(new,example_id='invented')]):
            jsonl(b,bad)
            with self.subTest(bad=bad),self.assertRaises(ValueError):
                R.compare_supervision(a,b,'SEQONLY')


class StageTest(Base):
    def setUp(self):
        super().setUp();self.out=self.root/'run'
        self.recipe=dict(supervision={k:'input' for k in ('labels','splits','label_policy','registry','go_obo','example_clusters')},
            original_table='original',disjoint_keys=['genome_id'],settings={k:{} for k in F.KINDS})
        self.protocol=dict(recipes={t:copy.deepcopy(self.recipe) for t in G.MODES},pending_stages=['full_reference_ec_endpoint','production_replay'],
            limitations='x',genomes=['g'],gaps={'g':{'path':'gaps'}},claims={a:{'g':{}} for a in R.ARMS},contracts={a:[] for a in R.ARMS},
            plans={a:{c:dict(command=['x'],manifest='m',inference='i',features=str(self.out/a/c/'features')) for c in R.CONTEXTS} for a in R.ARMS},
            zero_audit=dict(prefixes=list(R.ZERO_PREFIXES)))

    def test_export_audits_each_arm_with_candidate_specific_audit_and_requires_zero_context(self):
        for a in R.ARMS:(self.out/a).mkdir(parents=True)
        with ExitStack() as stack:
            stack.enter_context(patch.object(R,'verify',return_value=self.protocol))
            batch=stack.enter_context(patch.object(F.V,'run_batch'))
            audit=stack.enter_context(patch.object(S,'audit_export',return_value={'g':2}))
            legacy=stack.enter_context(patch.object(F,'audit_export',side_effect=AssertionError('legacy audit used')))
            zero=stack.enter_context(patch.object(R,'zero_audit',return_value=4))
            stack.enter_context(patch.object(R,'sequence_stage_audit',return_value={}))
            stack.enter_context(patch.object(G.R,'jsonl',return_value=[]))
            stack.enter_context(patch.object(R,'compare_candidates',return_value={}))
            stack.enter_context(patch.object(F,'hashes',return_value={}))
            R.export(self.out)
            self.assertEqual(batch.call_count,20);self.assertEqual(audit.call_count,20);self.assertEqual(zero.call_count,5)
            self.assertEqual({c.args[1].__class__ for c in audit.call_args_list},{dict})
            with self.assertRaises(FileExistsError):R.export(self.out)
            legacy.assert_not_called()

    def test_export_rejects_missing_genome_and_unverified_extras(self):
        self.out.mkdir()
        for result in ({},{'h':1},ValueError('candidate/claim identities differ')):
            with self.subTest(result),patch.object(R,'verify',return_value=self.protocol),patch.object(F.V,'run_batch'),\
                    patch.object(S,'audit_export',**({'side_effect':result} if isinstance(result,Exception) else {'return_value':result})),\
                    self.assertRaises(ValueError):R.export(self.out)

    def test_forty_models_freeze_before_any_test_evaluation(self):
        self.out.mkdir();(self.out/'export.json').write_text(json.dumps(dict(status='verified',outputs={})))
        for a in R.ARMS:
            for c in R.CONTEXTS:
                f=self.out/a/c/'features/g';f.mkdir(parents=True);(f/'feature_contract.json').write_text('{}')
        rows=[{'split':s} for s in F.T.SPLITS]
        def build(*args,**kwargs):
            dest=args[5];dest.mkdir();(dest/'training_table.jsonl').write_text('{}\n');(dest/'training_table_manifest.json').write_text('{}')
            return {'table':str(dest/'training_table.jsonl')}
        def fit(train,calibration,*args,**kwargs):
            self.assertTrue(all(r['split']=='train' for r in train));self.assertTrue(all(r['split']=='calibration' for r in calibration))
            self.assertFalse((self.out/'training.json').exists());return {'stub':True}
        def score(model,contract,test,*args):
            self.assertEqual(len(list(self.out.glob('*/*/models/model_*.json'))),40);self.assertEqual(test,[{'split':'test'}]);return {'n':1}
        with ExitStack() as stack:
            for obj,name,value in [(R,'verify',self.protocol),(R,'compare_supervision',dict(original_only=0,new_only=0)),(F.T,'load_contract',{}),(F.T,'read_table',rows),
                    (F.T,'check_disjoint',None),(F.T,'split_identity',{})]:stack.enter_context(patch.object(obj,name,return_value=value))
            stack.enter_context(patch.object(F.P,'build_table',side_effect=build))
            train=stack.enter_context(patch.object(F.T,'train',side_effect=fit))
            stack.enter_context(patch.object(F.T,'write_frozen',side_effect=lambda artifact,p:p.write_text(json.dumps(artifact))))
            test=stack.enter_context(patch.object(F.T,'evaluate_frozen',side_effect=score))
            R.train(self.out);self.assertEqual(train.call_count,40);test.assert_not_called()
            R.evaluate(self.out);self.assertEqual(test.call_count,40)
            report=json.loads((self.out/'evaluation.json').read_text())
            self.assertEqual(report['pending_stages'],{'full_reference_ec_endpoint':'not_run','production_replay':'not_run'})
            with self.assertRaises(FileExistsError):R.evaluate(self.out)

    def test_incomplete_freeze_rejected_before_test_rows_are_read(self):
        self.out.mkdir();(self.out/'training.json').write_text(json.dumps({'status':'verified','outputs':{}}))
        with patch.object(R,'verify',return_value=self.protocol),patch.object(F.T,'read_table') as read:
            with self.assertRaisesRegex(ValueError,'missing frozen model'):R.evaluate(self.out)
            read.assert_not_called()


if __name__=='__main__':unittest.main()
