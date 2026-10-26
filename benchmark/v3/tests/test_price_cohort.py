import json
import tempfile
import unittest
from pathlib import Path

from benchmark.v3.prepare_price_cohort import CohortError, prepare, split_namespace
from benchmark.v3.prepare_price_layout import prepare as prepare_layout

PROTEINS = {
    # two strains of one species, a second species of the same genus, and a genome whose taxon is not in the XML
    'GCF_1.1': (100, ['p1', 'p2', 'p3']),
    'GCF_2.1': (101, ['p1', 'q2']),          # p1 is the same accession string as in GCF_1.1: a different protein
    'GCF_3.1': (200, ['r1', 'r2']),
    'GCA_4.1': (300, ['s__1', 's2']),         # '__' inside a protein ID must survive the namespace split
}

TAXON = '''<Taxon><TaxId>{tid}</TaxId><ScientificName>{name}</ScientificName><Rank>{rank}</Rank>
<Lineage>x</Lineage><LineageEx>
<Taxon><TaxId>2</TaxId><ScientificName>Bacteria</ScientificName><Rank>domain</Rank></Taxon>
<Taxon><TaxId>{genus}</TaxId><ScientificName>Genus</ScientificName><Rank>genus</Rank></Taxon>{species}
</LineageEx></Taxon>'''


def taxon(tid, rank, genus, species):
    sp = f'<Taxon><TaxId>{species}</TaxId><ScientificName>Sp</ScientificName><Rank>species</Rank></Taxon>' if species else ''
    return TAXON.format(tid=tid, name=f'taxon{tid}', rank=rank, genus=genus, species=sp)


class CohortBase(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        genomes = self.root / 'genomes'
        for asm, (taxid, ids) in PROTEINS.items():
            folder = genomes / asm
            folder.mkdir(parents=True)
            gff = f'#!genome-build-accession NCBI_Assembly:{asm}\n##sequence-region chr 1 100000\n'
            for i, pid in enumerate(ids):
                gff += f'chr\tNCBI\tCDS\t{1 + 300 * i}\t{300 + 300 * i}\t.\t+\t0\tprotein_id={pid};product=SECRET_FUNCTION\n'
            (folder / 'genomic.gff').write_text(gff)
            (folder / 'protein.faa').write_text(''.join(f'>{pid} SECRET_FUNCTION\nMKV{"A" * (i + 3)}\n' for i, pid in enumerate(ids)))
            (folder / 'assembly_data_report.json').write_text(json.dumps(
                {'accession': asm, 'organism': {'taxId': taxid}, 'assemblyInfo': {'assemblyLevel': 'Complete Genome'}}))
        self.layout = self.root / 'layout'
        prepare_layout(genomes, self.layout)
        # genus 10 holds species 11 (taxids 100, 101 are two strains of species 11) and species 200 itself
        self.taxonomy = self.write_taxonomy([taxon(100, 'strain', 10, 11), taxon(101, 'strain', 10, 11), taxon(200, 'species', 10, None)])
        self.log = self.root / 'mmseqs.log'
        self.log.write_text('easy-cluster /x/all.faa /x/out /x/tmp --min-seq-id 0.3 -c 0.8 --cov-mode 0\n\nMMseqs Version:  \t15-6f452\n')

    def write_taxonomy(self, taxa):
        path = self.root / f'taxonomy{len(list(self.root.glob("taxonomy*")))}.xml'
        path.write_text('<?xml version="1.0" ?><TaxaSet>' + ''.join(taxa) + '</TaxaSet>')
        return path

    def all_ids(self):
        return [f'{asm}__{p}' for asm, (_, ids) in PROTEINS.items() for p in ids]

    def write_clusters(self, groups, name='clusters.tsv'):
        path = self.root / name
        path.write_text(''.join(f'{g[0]}\t{m}\n' for g in groups for m in g))
        return path

    def default_groups(self):
        # GCF_1.1__p1 and p2 are in-genome paralogs of one family; GCF_2.1__p1 is NOT their family despite the same accession
        return [['GCF_1.1__p1', 'GCF_1.1__p2', 'GCF_2.1__q2'], ['GCF_1.1__p3'], ['GCF_2.1__p1'],
                ['GCF_3.1__r1', 'GCF_3.1__r2'], ['GCA_4.1__s__1'], ['GCA_4.1__s2']]

    def run_prepare(self, groups=None, clusters=None, **kw):
        clusters = clusters or self.write_clusters(groups or self.default_groups())
        kw.setdefault('mmseqs_log', self.log)
        return prepare(self.layout, clusters, kw.pop('taxonomy', self.taxonomy), self.root / f'out{len(list(self.root.glob("out*")))}', **kw)

    def last_out(self):
        return sorted(self.root.glob('out*'))[-1]


class MembershipTest(CohortBase):
    def test_paralogs_share_a_family_and_stay_distinct_proteins(self):
        report = self.run_prepare()
        rows = [l.split('\t') for l in (self.last_out() / 'cohort_orthogroups.tsv').read_text().splitlines()]
        self.assertEqual(rows[0], ['genome_id', 'protein_id', 'orthogroup_id', 'subfamily_id', 'weight'])
        fam = {(r[0], r[1]): r[2] for r in rows[1:]}
        self.assertEqual(fam[('GCF_1.1', 'p1')], fam[('GCF_1.1', 'p2')])
        self.assertEqual(fam[('GCF_1.1', 'p1')], fam[('GCF_2.1', 'q2')])
        self.assertNotEqual(fam[('GCF_1.1', 'p1')], fam[('GCF_2.1', 'p1')])
        self.assertEqual(report['family_stats']['families_with_in_genome_paralogs'], 2)  # p1/p2 and r1/r2
        self.assertTrue(all(len(r) == 5 and r[3] == '-' and r[4] == '1' for r in rows[1:]))
        self.assertTrue(all(':' not in r[2] for r in rows[1:]))  # normalizeFamily would strip a namespace prefix
        self.assertEqual(len(rows) - 1, 9)

    def test_namespace_splits_once_and_unnamespaced_ids_are_rejected(self):
        self.assertEqual(split_namespace('GCA_4.1__s__1'), ('GCA_4.1', 's__1'))
        with self.assertRaises(CohortError):
            split_namespace('p1')
        rows = (self.run_prepare() and (self.last_out() / 'cohort_orthogroups.tsv').read_text())
        self.assertIn('GCA_4.1\ts__1\t', rows)

    def test_complete_exactly_once_membership_certifies_every_genome(self):
        report = self.run_prepare()
        meta = json.loads((self.last_out() / 'cohort_metadata.json').read_text())
        self.assertTrue(report['certified_family_search'])
        for g in meta['genomes'].values():
            self.assertIs(g['familySearchComplete'], True)
            self.assertEqual(g['familySearchCoverage'], 1.0)
        self.assertIn('not established orthologs', report['interpretation'])

    def test_partial_coverage_is_rejected_and_never_certified(self):
        groups = [g for g in self.default_groups() if g != ['GCF_3.1__r1', 'GCF_3.1__r2']]
        with self.assertRaisesRegex(CohortError, 'refusing to certify'):
            self.run_prepare(groups)
        report = self.run_prepare(groups, allow_partial=True)
        meta = json.loads((self.last_out() / 'cohort_metadata.json').read_text())['genomes']
        self.assertFalse(report['certified_family_search'])
        self.assertFalse(meta['GCF_3.1']['familySearchComplete'])
        self.assertEqual(meta['GCF_3.1']['familySearchCoverage'], 0.0)
        self.assertFalse(meta['GCF_1.1']['familySearchComplete'])  # no certification from an incomplete cluster file
        self.assertEqual(report['cluster_audit']['missing_proteins'], 2)

    def test_wrong_ids_and_repeated_members_are_rejected(self):
        with self.assertRaisesRegex(CohortError, 'foreign'):
            self.run_prepare(self.default_groups() + [['GCF_9.9__p1']])
        with self.assertRaisesRegex(CohortError, 'foreign'):
            self.run_prepare(self.default_groups() + [['GCF_1.1__p9']])      # genome known, protein never acquired
        with self.assertRaisesRegex(CohortError, 'foreign'):
            self.run_prepare(self.default_groups() + [['p1']])               # namespace stripped
        with self.assertRaisesRegex(CohortError, 'repeated'):
            self.run_prepare(self.default_groups() + [['GCF_3.1__r1']])      # one protein in two families
        with self.assertRaisesRegex(CohortError, 'foreign'):
            self.run_prepare(self.default_groups() + [['GCF_1.1__p1.1']])    # a different protein version is never joined


class TaxonomyTest(CohortBase):
    def metadata(self, **kw):
        self.run_prepare(**kw)
        return json.loads((self.last_out() / 'cohort_metadata.json').read_text())['genomes']

    def test_strains_share_species_block_and_species_share_clade_not_cluster(self):
        meta = self.metadata(taxonomy=self.taxonomy)
        self.assertEqual(meta['GCF_1.1']['strainCluster'], meta['GCF_2.1']['strainCluster'])
        self.assertEqual(meta['GCF_1.1']['clade'], meta['GCF_3.1']['clade'])
        self.assertNotEqual(meta['GCF_1.1']['strainCluster'], meta['GCF_3.1']['strainCluster'])
        self.assertEqual(meta['GCF_1.1']['clade'], 'ncbi_genus_10')
        self.assertEqual(meta['GCF_3.1']['strainCluster'], 'ncbi_species_200')

    def test_missing_taxonomy_stays_unresolved_not_guessed(self):
        meta = self.metadata()
        self.assertNotIn('clade', meta['GCA_4.1'])
        self.assertNotIn('strainCluster', meta['GCA_4.1'])
        report = json.loads((self.last_out() / 'report.json').read_text())
        self.assertIn('GCA_4.1', report['taxonomy']['unresolved_genomes'])
        self.assertEqual(report['taxonomy']['absent_from_xml'], [300])

    def test_missing_rank_stays_unknown(self):
        tax = self.write_taxonomy([taxon(100, 'strain', 10, None), taxon(101, 'strain', 10, 11), taxon(200, 'species', 10, None)])
        meta = self.metadata(taxonomy=tax)
        self.assertEqual(meta['GCF_1.1']['clade'], 'ncbi_genus_10')
        self.assertNotIn('strainCluster', meta['GCF_1.1'])

    def test_top_level_taxid_must_be_the_requested_one(self):
        tax = self.write_taxonomy([taxon(100, 'strain', 10, 11), taxon(101, 'strain', 10, 11), taxon(999, 'species', 10, None)])
        with self.assertRaisesRegex(CohortError, 'was not requested'):
            self.run_prepare(taxonomy=tax)

    def test_metadata_has_only_schema_keys_and_no_functional_text(self):
        self.run_prepare()
        text = (self.last_out() / 'cohort_metadata.json').read_text()
        allowed = {'clade', 'strainCluster', 'assemblyStatus', 'completenessFraction', 'familySearchComplete',
                   'familySearchCoverage', 'contigs'}
        for g in json.loads(text)['genomes'].values():
            self.assertLessEqual(set(g), allowed)
            self.assertNotIn('completenessFraction', g)
        for name in ('cohort_metadata.json', 'cohort_orthogroups.tsv', 'families.tsv', 'report.json'):
            self.assertNotIn('SECRET_FUNCTION', (self.last_out() / name).read_text())
        self.assertEqual(json.loads(text)['genomes']['GCF_1.1']['assemblyStatus'], 'COMPLETE')


class ProvenanceAndSubfamilyTest(CohortBase):
    def test_hashes_protocol_and_coverage_are_recorded(self):
        report = self.run_prepare()
        self.assertEqual(report['mmseqs']['version'], '15-6f452')
        self.assertTrue(report['mmseqs']['command'].startswith('mmseqs easy-cluster'))
        self.assertEqual(report['mmseqs']['parameters']['min-seq-id'], '0.3')
        self.assertEqual(len(report['inputs']['clusters_sha256']), 64)
        self.assertEqual(len(report['inputs']['taxonomy_sha256']), 64)
        self.assertEqual(report['cluster_audit']['expected_proteins'], 9)
        self.assertFalse(report['subfamilies']['available'])

    def test_protocol_is_required(self):
        with self.assertRaisesRegex(CohortError, 'version and command'):
            self.run_prepare(mmseqs_log=None)

    def test_tampered_layout_is_refused(self):
        with (self.layout / 'all_proteins.faa').open('a') as fh:
            fh.write('>GCF_1.1__p1\nMKV\n')
        with self.assertRaisesRegex(CohortError, 'digest'):
            self.run_prepare()

    def test_nested_subfamilies_are_hashed_and_non_nested_ones_need_explicit_policy(self):
        sub_ok = self.write_clusters([['GCF_1.1__p1', 'GCF_1.1__p2'], ['GCF_2.1__q2'], ['GCF_1.1__p3'], ['GCF_2.1__p1'],
                                      ['GCF_3.1__r1'], ['GCF_3.1__r2'], ['GCA_4.1__s__1'], ['GCA_4.1__s2']], 'sub70.tsv')
        report = self.run_prepare(subclusters=sub_ok)
        self.assertEqual(report['subfamilies']['nesting_violations'], 0)
        rows = [l.split('\t') for l in (self.last_out() / 'cohort_orthogroups.tsv').read_text().splitlines()[1:]]
        sub = {(r[0], r[1]): r[3] for r in rows}
        self.assertEqual(sub[('GCF_1.1', 'p1')], sub[('GCF_1.1', 'p2')])
        self.assertNotEqual(sub[('GCF_1.1', 'p1')], sub[('GCF_2.1', 'q2')])
        # a 70% subcluster that straddles two 30% families
        sub_bad = self.write_clusters([['GCF_1.1__p3', 'GCF_2.1__p1'], ['GCF_1.1__p1', 'GCF_1.1__p2'], ['GCF_2.1__q2'],
                                       ['GCF_3.1__r1'], ['GCF_3.1__r2'], ['GCA_4.1__s__1'], ['GCA_4.1__s2']], 'sub70bad.tsv')
        with self.assertRaisesRegex(CohortError, 'nesting'):
            self.run_prepare(subclusters=sub_bad)
        report = self.run_prepare(subclusters=sub_bad, subfamily_policy='combined')
        self.assertEqual(report['subfamilies']['nesting_violations'], 1)
        self.assertTrue(report['subfamilies']['warnings'])
        rows = [l.split('\t') for l in (self.last_out() / 'cohort_orthogroups.tsv').read_text().splitlines()[1:]]
        sub = {(r[0], r[1]): r[3] for r in rows}
        self.assertNotEqual(sub[('GCF_1.1', 'p3')], sub[('GCF_2.1', 'p1')])  # combined signature keeps families apart

    def test_incomplete_subcluster_file_is_rejected(self):
        sub = self.write_clusters([['GCF_1.1__p1']], 'sub70short.tsv')
        with self.assertRaisesRegex(CohortError, 'exact-once'):
            self.run_prepare(subclusters=sub)


if __name__ == '__main__':
    unittest.main()
