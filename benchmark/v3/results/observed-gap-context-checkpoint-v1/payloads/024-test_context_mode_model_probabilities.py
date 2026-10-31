import json,os
from pathlib import Path
import sys,tempfile,unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import check_context_mode_models as C

class ContextModeModelProbabilityTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp')))
        self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.exports=self.root/'exports';arm=self.exports/'UNIVERSAL';(arm/'models').mkdir(parents=True);(arm/'table').mkdir()
        for kind in ('logistic','gbdt'):(arm/'models'/f'model_{kind}.json').write_text(json.dumps({'contract':{'features':['x']}}))
        (arm/'table/training_table.jsonl').write_text(json.dumps(dict(example_id='e',features={'x':.2},sequence_sources=[]))+'\n')
        self.jar=self.root/'engine.jar';self.jar.write_bytes(b'test jar; external reader is mocked')
        self.counter=0
    def compare(self,expected,observed):
        def invoke(command,**kwargs):
            Path(command[-1]).write_text(json.dumps(dict(example_id='e',mode='UNIVERSAL',unsupported_reason=None,probability=observed))+'\n')
        self.counter+=1
        with patch.object(C.subprocess,'run',side_effect=invoke),patch.object(C.T,'unsupported_input_reason',return_value=None),patch.object(C.T,'predict_proba',return_value=expected):
            return C.run(self.exports,self.jar,self.root/f'result-{self.counter}',('UNIVERSAL',))
    def test_nonfinite_boolean_and_out_of_range_probabilities_fail_even_if_difference_looks_zero(self):
        for expected,observed in [(float('nan'),.4),(.4,float('nan')),(1.,True),(-.1,-.1),(1.1,1.1)]:
            with self.subTest(expected=expected,observed=observed),self.assertRaises(ValueError):
                self.compare(expected,observed)
    def test_valid_probability_parity_remains_accepted(self):
        result=self.compare(.4,.4)
        self.assertEqual(result['status'],'verified')
        self.assertTrue(all(r['max_probability_error']==0 for r in result['results'].values()))

if __name__=='__main__':unittest.main()
