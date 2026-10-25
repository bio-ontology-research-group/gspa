"""Operon labels and coverage must reflect functional evidence, not placeholders or duplicate genes."""
import importlib.util
from pathlib import Path
import unittest

path=Path(__file__).resolve().parents[2]/'gspa-cli/src/main/resources/visualize/make_viz.py'
spec=importlib.util.spec_from_file_location('viz',path)
viz=importlib.util.module_from_spec(spec);spec.loader.exec_module(viz)

class OperonEvidenceTest(unittest.TestCase):
    def test_placeholder_products_are_not_operon_names(self):
        for product in ('Predicted CDS','hypothetical protein','Uncharacterized protein'):
            ops=[{'members':['p1','p2']}]
            viz.name_operons(ops,[],{}, {p:{'product':product} for p in ('p1','p2')},[])
            self.assertIsNone(ops[0]['name'])

    def test_duplicate_proteins_do_not_inflate_pathway_coverage_and_ec_counts(self):
        ops=[{'members':['p1','p2','p3']}]
        annos=[{'protein_id':p,'term':'EC:1.1.1.1','aspect':'F','post':.9} for p in ops[0]['members']]
        pathway={'id':'KEGG:M1','name':'Test','required':2,'present_terms':['GO:1'],'missing_terms':['GO:2']}
        viz.name_operons(ops,annos,{}, {},[pathway],{'EC:1.1.1.1':'GO:1'})
        self.assertEqual(ops[0]['pathways'][0]['coverage'],.5)
        self.assertEqual(ops[0]['pathways'][0]['n_members'],3)

if __name__=='__main__':unittest.main()
