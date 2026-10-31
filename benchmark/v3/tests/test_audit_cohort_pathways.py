import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import audit_cohort_pathways as A


class CohortObservationAudit(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='gspa-pathway-audit-')
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)/'search';self.root.mkdir()
        self.row=dict(genome_id='g',taxonomy='Archaea',proteins=1,proteome='/remote/proteomes/g.faa',sha256='a'*64)
        self.manifest=dict(genomes=[self.row],threads=2,
                           sources={'/remote/observe_gapsmith.py':'b'*64,
                                    '/remote/observe_gapsmith_cohort.py':'c'*64,
                                    '/remote/binary':'d'*64,'/remote/diamond':'e'*64})
        self.write(self.root/'manifest.json',self.manifest)
        self.write(self.root/'staging-hashes.json',{'manifest.json':A.digest(self.root/'manifest.json')})
        self.directory=self.root/'observations/g';self.directory.mkdir(parents=True)
        output=self.directory/'outputs';output.mkdir()
        (output/'g-all-Reactions.tbl').write_text('status\ngood_blast\n')
        (output/'g-all-Pathways.tbl').write_text('ID\tPrediction\nP\tTRUE\n')
        self.refs={name:'f'*64 for name in A.TABLES}
        for category in ['rev','unrev','rxn']:
            self.refs[f'seq/Archaea/{category}/r.fasta']='f'*64
        self.write(self.directory/'reference-hashes.json',self.refs)
        self.observation=dict(status='observed',inputs={p:h for p,h in self.manifest['sources'].items()
                             if not p.endswith('observe_gapsmith_cohort.py')},
                             taxonomy='Archaea',threads=2,exit_code=0,job_id='fixture',
                             output_hashes={str(p.relative_to(self.directory)):A.digest(p) for p in output.iterdir()},
                             reaction_columns=['status'],reaction_statuses={'good_blast':1},pathways=1,
                             predicted_pathways=1,reference_files=len(self.refs))
        self.observation['inputs'][self.row['proteome']]=self.row['sha256']
        self.rebind()

    @staticmethod
    def write(path,data): path.write_text(json.dumps(data,sort_keys=True)+'\n')

    def rebind(self):
        self.write(self.directory/'observation.json',self.observation)
        self.write(self.directory/'execution.json',dict(status='observed',genome_id='g',
                   manifest_sha256=A.digest(self.root/'manifest.json'),sources=self.manifest['sources'],
                   observation_sha256=A.digest(self.directory/'observation.json')))

    def test_complete_bound_search_is_verified_and_failed_attempt_is_missing(self):
        result=A.audit([self.root]);self.assertEqual(result['status'],'verified')
        self.assertEqual(result['observations']['g']['reaction_rows'],1)
        self.write(self.directory/'execution.json',{'status':'failed'})
        result=A.audit([self.root]);self.assertEqual(result['status'],'incomplete')
        self.assertEqual(result['missing_genomes'],['g'])
        self.assertEqual(result['attempts'][0]['status'],'failed')

    def test_changed_output_is_rejected(self):
        (self.directory/'outputs/g-all-Reactions.tbl').write_text('status\nno_blast\n')
        with self.assertRaisesRegex(ValueError,'output bytes differ'): A.audit([self.root])

    def test_rehashed_receipt_cannot_substitute_another_proteome(self):
        self.observation['inputs'][self.row['proteome']]='0'*64;self.rebind()
        with self.assertRaisesRegex(ValueError,'bindings differ'): A.audit([self.root])

    def test_rehashed_receipt_cannot_invent_counts(self):
        self.observation['pathways']=100;self.rebind()
        with self.assertRaisesRegex(ValueError,'summary differs'): A.audit([self.root])

    def test_wrong_taxonomy_reference_is_rejected(self):
        self.refs.pop('seq/Archaea/rev/r.fasta');self.refs['seq/Bacteria/rev/r.fasta']='f'*64
        self.write(self.directory/'reference-hashes.json',self.refs)
        with self.assertRaisesRegex(ValueError,'wrong taxonomy'): A.audit([self.root])

    def test_search_must_match_the_verified_acquired_reference(self):
        acquisition=self.root.parent/'acquisition';acquisition.mkdir()
        hashes=acquisition/'reference-hashes.json';self.write(hashes,self.refs)
        receipt=acquisition/'acquisition.json'
        self.write(receipt,dict(status='verified',reference_manifest_sha256=A.digest(hashes)))
        self.manifest['reference_acquisition']='/remote/acquisition/acquisition.json'
        self.write(self.root/'manifest.json',self.manifest)
        self.write(self.root/'staging-hashes.json',{'manifest.json':A.digest(self.root/'manifest.json')})
        self.rebind()
        self.assertEqual(A.audit([self.root])['status'],'verified')
        changed=dict(self.refs);changed['seq/Archaea/rev/r.fasta']='0'*64
        self.write(hashes,changed)
        self.write(receipt,dict(status='verified',reference_manifest_sha256=A.digest(hashes)))
        with self.assertRaisesRegex(ValueError,'reference differs'): A.audit([self.root])


if __name__=='__main__': unittest.main()
