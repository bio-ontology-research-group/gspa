"""Identity/range coverage mechanics; synthetic contracts are not calibration evidence."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from test_model_input_support import T, row, tiny_contract, contracts


class PredictorContractSupportTests(unittest.TestCase):
    def setUp(self):
        self.c = tiny_contract(('fn.type_ec', 'role.catalytic', 'seq.src.other.reported', 'seq.src.other.max_raw'))
        self.sources = ['tool_a']
        self.contracts = contracts(self.sources)
        self.rows = [dict(row(i, split, split, split, label, [1,1,1,score]),
                          sequence_sources=self.sources, sequence_contracts=copy.deepcopy(self.contracts))
                     for i,(split,label,score) in enumerate([('train',0,.2),('train',1,.8),('calibration',0,.1),('calibration',1,.9)])]
        self.art = dict(intended_use='production',contract={'features':self.c['names']},
                        input_support=T.build_input_support(self.rows[:2],self.rows[2:],self.c['names']))

    def reason(self, x=None, contracts=None):
        return T.unsupported_input_reason(self.art,self.rows[0]['x'] if x is None else x,self.sources,
                                           sequence_contracts=self.contracts if contracts is None else contracts)

    def test_identity_mutations_with_identical_vectors_are_unsupported(self):
        self.assertIsNone(self.reason())
        for key in ('tool_version','tool_sha256','reference_sha256','parameters_sha256','adapter_sha256','score_model'):
            changed=copy.deepcopy(self.contracts)
            changed[0]['contract'][key]='e'*64 if key.endswith('_sha256') else 'changed'
            self.assertEqual(self.reason(contracts=changed),'unseen_predictor_contract')

    def test_ranges_include_calibration_endpoints_but_no_unseen_extremes(self):
        profile=self.art['input_support']['profiles'][0]
        self.assertEqual((profile['lower'],profile['upper']),([.1],[.9]))
        for score,expected in [(0.,'sequence_score_outside_training_calibration_range'),(.1,None),(.9,None),(1.,'sequence_score_outside_training_calibration_range')]:
            x=self.rows[0]['x'].copy();x[-1]=score
            self.assertEqual(self.reason(x),expected)

    def test_missing_contracts_and_legacy_profiles_do_not_grant_coverage(self):
        self.assertEqual(T.unsupported_input_reason(self.art,self.rows[0]['x'],self.sources),'candidate_predictor_contract_unavailable')
        self.assertEqual(self.reason(contracts=[dict(source='tool_a',contract=None)]),'candidate_predictor_contract_unavailable')
        for schema in ('gspa.model-input-support/1','gspa.model-input-support/2'):
            art=copy.deepcopy(self.art);art['input_support']['schema']=schema
            self.assertEqual(T.unsupported_input_reason(art,self.rows[0]['x'],self.sources,sequence_contracts=self.contracts),
                             'model_predictor_contract_support_unavailable')

    def test_new_version_cannot_borrow_other_versions_class_coverage(self):
        rows=copy.deepcopy(self.rows)
        rows[-1]['sequence_contracts'][0]['contract']['tool_version']='2'
        art=dict(self.art,input_support=T.build_input_support(rows[:2],rows[2:],self.c['names']))
        self.assertEqual(len(art['input_support']['profiles']),2)
        for r in (rows[0],rows[-1]):
            self.assertEqual(T.unsupported_input_reason(art,r['x'],r['sequence_sources'],sequence_contracts=r['sequence_contracts']),
                             'profile_lacks_both_classes_in_train_and_calibration')

    def test_contracts_roundtrip_and_bind_training_digest_without_mutable_aliases(self):
        before=T.training_data_digest(self.rows[:2],self.rows[2:])
        changed=copy.deepcopy(self.rows);changed[0]['sequence_contracts'][0]['contract']['tool_version']='2'
        self.assertNotEqual(before,T.training_data_digest(changed[:2],changed[2:]))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'rows.jsonl';T.write_rows_jsonl(self.rows,self.c['names'],p)
            self.assertEqual([r['sequence_contracts'] for r in T.read_table(p,self.c)],[self.contracts]*4)
        self.rows[0]['sequence_contracts'][0]['contract']['tool_version']='altered'
        self.assertEqual(self.art['input_support']['profiles'][0]['sequence_contracts'][0]['contract']['tool_version'],'1')

    def test_corrupt_profile_shapes_ranges_and_duplicate_identities_fail(self):
        for kind in ('missing_bounds','nan','reversed','duplicate','float_duplicate','omitted_feature','incomplete_contract','availability'):
            art=copy.deepcopy(self.art);p=art['input_support']['profiles'][0]
            if kind=='missing_bounds':p['lower']=[]
            if kind=='nan':p['lower'][0]=float('nan')
            if kind=='reversed':p['lower'][0]=2
            if kind=='duplicate':art['input_support']['profiles'].append(copy.deepcopy(p))
            if kind=='float_duplicate':art['input_support']['profiles'].append(dict(copy.deepcopy(p),flags=list(map(float,p['flags']))))
            if kind=='availability':art['input_support']['available']='true'
            if kind=='omitted_feature':art['input_support']['score_features']=[]
            if kind=='incomplete_contract':p['sequence_contracts'][0]['contract']=None
            with self.subTest(kind=kind),self.assertRaises(T.ProvenanceError):
                T.unsupported_input_reason(art,self.rows[0]['x'],self.sources,sequence_contracts=self.contracts)

    def test_one_unsupported_member_prevents_partial_locus_probability(self):
        xs=[self.rows[0]['x'].copy(),self.rows[0]['x'].copy()];xs[1][-1]=1.
        with self.assertRaisesRegex(T.ProvenanceError,'outside_training'):
            T.assign_locus(self.art,xs,[self.sources]*2,sequence_contracts=[self.contracts]*2)

    def test_unknown_training_identities_cannot_be_silently_certified(self):
        rows=copy.deepcopy(self.rows);rows[0].pop('sequence_contracts')
        self.assertFalse(T.build_input_support(rows[:2],rows[2:],self.c['names'])['available'])
        malformed=copy.deepcopy(self.contracts);malformed[0]['contract']['reference_sha256']='not-a-hash'
        with self.assertRaises(T.ProvenanceError):T.validate_sequence_contracts(malformed,self.sources)


if __name__=='__main__':unittest.main()
