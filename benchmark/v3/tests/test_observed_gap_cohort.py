import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import run_observed_gap_cohort as G


class ObservedGapCohortTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp')))
        self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)

    def write(self,path,data):
        path.parent.mkdir(parents=True,exist_ok=True);G.F.save(path,data);return path

    def replay(self):
        root=self.root/'replay';annotation=root/'g/annotation'
        gap=dict(pathway_id='gapsmith:p',reaction_id='r',gapseq_guessed=False)
        path=self.write(annotation/'g_integrated.tsv.gaps.jsonl',gap)
        obs=self.write(annotation/'g_gapsmith_observations.json',dict(
            reaction_rows=[dict(pathway='p',rxn='r',status='no_blast',spont='FALSE',is_complex='FALSE')],
            pathway_rows=[dict(ID='p',Prediction='TRUE')]))
        # JSONL gap records must be single-line, unlike the JSON receipts.
        path.write_text(json.dumps(gap)+'\n')
        receipt=dict(status='verified',acquisition_status='verified',missing_genomes=[],
            replayed_genomes=1,expected_genomes=1,summaries={'g':dict(gap_records=1)},inputs={},
            outputs={str(p.relative_to(root)):G.F.T.sha256_file(p) for p in [path,obs]})
        self.write(root/'verification.json',receipt)
        return root,path,receipt

    def test_observed_gaps_require_complete_receipt_and_unchanged_outputs(self):
        root,path,receipt=self.replay()
        gaps,_=G.verified_gaps(root);self.assertEqual(gaps['g']['records'],1)
        path.write_text(path.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'changed'):G.verified_gaps(root)
        receipt['status']='partial_verified';self.write(root/'verification.json',receipt)
        with self.assertRaisesRegex(ValueError,'complete verified'):G.verified_gaps(root)

    def test_resealed_gap_must_still_match_complete_search_eligibility(self):
        root,path,receipt=self.replay()
        row=json.loads(path.read_text());row['reaction_id']='invented';path.write_text(json.dumps(row)+'\n')
        receipt['outputs'][str(path.relative_to(root))]=G.F.T.sha256_file(path)
        self.write(root/'verification.json',receipt)
        with self.assertRaisesRegex(ValueError,'lost or invented'):G.verified_gaps(root)

    def feature_pair(self):
        x=dict(candidate_key='p|GO|GO:1|UNSPECIFIED',genome_id='g',sequence_sources=['s'],
            sequence_contracts=[{'source':'s','contract':{'tool':'fixture'}}],
            cohort_context_mode='UNIVERSAL',quality_evidence=None,
            features={'seq.score':.25,'cohort.context':.5,'ctx.gap.available':0.,'ctx.gap.function_listed':0.})
        y=copy.deepcopy(x);y['features']['ctx.gap.available']=1.;y['features']['ctx.gap.function_listed']=1.
        a=self.root/'a.jsonl';b=self.root/'b.jsonl'
        a.write_text(json.dumps(x)+'\n');b.write_text(json.dumps(y)+'\n');return a,b,x,y

    def test_paired_features_retain_all_candidates_and_only_gap_columns_change(self):
        a,b,x,y=self.feature_pair()
        self.assertEqual(G.paired_features(a,b),dict(candidates=1,changed_columns={
            'ctx.gap.available':1,'ctx.gap.function_listed':1}))
        for key in ('seq.score','cohort.context'):
            altered=copy.deepcopy(y);altered['features'][key]+=.1;b.write_text(json.dumps(altered)+'\n')
            with self.assertRaisesRegex(ValueError,'non-gap'):G.paired_features(a,b)
        b.write_text('')
        with self.assertRaisesRegex(ValueError,'count changed'):G.paired_features(a,b)

    def test_predictor_identity_and_gap_availability_cannot_be_swapped(self):
        a,b,x,y=self.feature_pair()
        y['sequence_contracts'][0]['contract']['tool']='other';b.write_text(json.dumps(y)+'\n')
        with self.assertRaisesRegex(ValueError,'provenance'):G.paired_features(a,b)
        y=copy.deepcopy(x);b.write_text(json.dumps(y)+'\n')
        with self.assertRaisesRegex(ValueError,'availability'):G.paired_features(a,b)

    def test_registration_pairs_modes_without_models_or_labels_and_binds_gap_bytes(self):
        replay,gap,receipt=self.replay()
        claims=self.root/'claims';claims.mkdir()
        contracts={mode:{'identity':{'tool':mode}} for mode in G.MODES}
        self.write(claims/'contracts.json',contracts);self.write(claims/'inputs.json',{})
        files={}
        for mode in G.MODES:
            p=claims/(mode+'.jsonl');p.write_text('{}\n')
            files[mode]={'g':dict(path=str(p),sha256=G.F.T.sha256_file(p))}
        self.write(claims/'run.json',dict(status='verified',enriched_claims_reversible=True,output_claims=files,
            contract_manifest_sha256=G.F.T.sha256_file(claims/'contracts.json'),
            inputs_sha256=G.F.T.sha256_file(claims/'inputs.json')))
        jar=self.root/'engine.jar';jar.write_text('fixture')
        templates={}
        for mode in G.MODES:
            p=self.root/mode;p.mkdir();templates[mode]=p
            (p/'manifest.tsv').write_text('genome_id\tclaims_path\tout_path\tgaps_path\n'+'g\tunused\tunused\t-\n')
            self.write(p/'command.json',['java','-jar',str(jar),'integrate','--cohort'])
        out=self.root/'prepared';protocol=G.prepare(templates,claims,replay,jar,out)
        self.assertFalse(protocol['labels_read']);self.assertFalse(protocol['models_fitted'])
        for mode in G.MODES:
            for arm in G.ARMS:
                plan=protocol['plans'][mode][arm]
                with Path(plan['manifest']).open() as f:row=next(G.csv.DictReader(f,delimiter='\t'))
                self.assertEqual(row['claims_path'],files[mode]['g']['path'])
                self.assertEqual(row['gaps_path'],str(gap) if arm=='observed_gaps' else '-')
                self.assertEqual('--rlc-suggester' in plan['command'],arm=='observed_gaps')
                self.assertNotIn('--learned-model',plan['command'])
        gap.write_text(gap.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'changed'):G.verify(out)


if __name__=='__main__':unittest.main()
