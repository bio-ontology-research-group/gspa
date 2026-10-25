"""A browser must not invent measured QC or fail when the optional report is absent."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / 'gspa-cli/src/main/resources/visualize/make_viz.py'


class OptionalQualityTest(unittest.TestCase):
    def render(self, quality=None, fasta=None, sidecar=False, analyses=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'prokka_out').mkdir()
            (root/'gspa_out').mkdir()
            (root/'prokka_out/sample.tsv').write_text('locus_tag\tftype\tlength_bp\tproduct\np1\tCDS\t300\tTest protein\n')
            (root/'prokka_out/sample.gff').write_text('##gff-version 3\nctg\ttest\tCDS\t1\t300\t.\t+\t0\tID=p1;locus_tag=p1\n')
            (root/'gspa_out/integrated.tsv').write_text('protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\tn_supporting\np1\tGO\tGO:0003824\tMF\t0.8\t1\n')
            (root/'gspa_out/provenance.json').write_text('[{"supporting_sources":["diamond"]}]')
            (root/'go.obo').write_text('format-version: 1.2\n')
            (root/'ec2go').write_text('! no mappings\n')
            if analyses is not None:
                (root/'gspa_out/analyses.json').write_text(json.dumps(analyses))
            if quality is not None:
                (root/'gspa_out/quality_gspa.json').write_text(json.dumps(quality))
            if fasta is not None:
                (root/'genome.fna').write_text(fasta)
            if sidecar:
                (root/'clean_out').mkdir()
                (root/'clean_out/sample.clean.tsv').write_text('protein_id\tterm\tscore\tannotation_type\np1\tEC:1.2.3.4\t0.99\tEC\n')
            output = root/'browser.html'
            subprocess.run([sys.executable,str(SCRIPT)], check=True,capture_output=True,text=True,env={
                **os.environ,'GSPA_WORKDIR':str(root),'GSPA_GENOME_ID':'sample','GSPA_OUT':str(output),
                'GSPA_RUN_DIR':'gspa_out','GSPA_FASTA':str(root/'genome.fna'),'GSPA_GO_OBO':str(root/'go.obo'),'GSPA_EC2GO':str(root/'ec2go')})
            return json.loads(re.search(r'const DATA = (.*);',output.read_text()).group(1))

    def test_absent_quality_and_external_qc_are_null(self):
        data = self.render()
        self.assertFalse(data['quality']['available'])
        for key in ('composite','annotation_coverage','consistent','completeness_score'):
            self.assertIsNone(data['quality'][key])
        for key in ('size_bp','completeness','contamination','busco'):
            self.assertIsNone(data['genome'][key])
        self.assertEqual(len(data['annotations']),1)
        self.assertEqual(data['build']['tutorial'], 'sample')
        self.assertEqual(data['build']['predictors_run'], ['diamond'])

    def test_sequence_statistics_do_not_imply_measured_completeness(self):
        data = self.render(fasta='>c1\nGCAA\n>c2\nGCTT\n')
        self.assertEqual(data['genome']['size_bp'], 8)
        self.assertEqual(data['genome']['n_contigs'], 2)
        self.assertEqual(data['genome']['gc'], 50)
        self.assertIsNone(data['genome']['completeness'])

    def test_supplied_measurements_survive(self):
        quality = {'summary':{'annotated_proteins':1,'annotation_coverage':1,'composite_score':.8},
                   'completeness':{'score':.5,'present_count':1,'missing_count':1,'missing_functions':['GO:0000001']},
                   'coherence':{'process_coherence':1,'pathway_coherence':.4,'complex_coherence':.5},
                   'consistency':{'consistent':False,'violation_count':2}}
        data = self.render(quality)
        self.assertTrue(data['quality']['available'])
        self.assertEqual(data['quality']['composite'],.8)
        self.assertFalse(data['quality']['consistent'])
        self.assertEqual(data['quality']['missing_essentials'],['GO:0000001'])

    def test_quality_opt_out_does_not_hide_independently_requested_pathways(self):
        quality = {'summary':{'annotated_proteins':1,'annotation_coverage':1,'composite_score':.8},
                   'completeness':{'score':.5,'present_count':1,'missing_count':1,'missing_functions':[]},
                   'coherence':{'process_coherence':1,'pathway_coherence':.4,'complex_coherence':.5,'pathway_detail':[]},
                   'consistency':{'consistent':True,'violation_count':0}}
        data=self.render(quality,analyses={'quality':{'status':'skipped'},'pathways':{'status':'complete'}})
        self.assertFalse(data['quality']['available'])
        self.assertIsNone(data['quality']['composite'])
        self.assertEqual(data['quality']['pathway_coherence'],.4)

    def test_raw_sidecar_scores_do_not_inflate_integrated_summaries(self):
        data = self.render(sidecar=True)
        self.assertEqual(len(data['annotations']), 2)
        self.assertEqual(data['stats']['n_high'], 1)
        self.assertEqual(sum(data['stats']['posterior_hist']), 1)
        self.assertEqual([term['term'] for term in data['stats']['top_terms']], ['GO:0003824'])
        self.assertEqual(data['proteins'][0]['nt'], 1)


if __name__ == '__main__': unittest.main()
