package gspa.integration.cohort

import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import gspa.model.Strand
import spock.lang.Specification
import static gspa.integration.cohort.CohortFixtures.*

class CohortEvidenceIdentitySpec extends Specification {
    private static Map identities(CohortInput input) {
        def byGenome = input.claims.groupBy { it.genomeId }.collectEntries { g, cs ->
            [(g): cs.groupBy { it.proteinId }]
        }
        CohortEvidenceIdentity.fingerprints(input, byGenome)
    }

    private static CohortInput pair() {
        def input = newInput()
        ['A', 'B'].each { g ->
            addGenome(input, [id: g, clade: 'c', cluster: 's', status: 'COMPLETE', closed: true,
                              family: 'OG1', candidate: 'T', anchors: ['N1', 'N2']])
        }
        input
    }

    def "complete inputs match across arbitrary genome protein and contig renaming and input order"() {
        given:
        def input = pair()
        def g = input.genomes.B
        def layout = new GenomeLayout('B')
        g.layout.byProtein.values().toList().reverseEach { l ->
            layout.add(new ProteinLocus(genomeId: 'B', proteinId: 'new-' + l.proteinId,
                contig: 'renamed', start: l.start, end: l.end, strand: l.strand))
        }
        layout.finishLoading()
        g.layout = layout
        g.assembly.contigs.renamed = g.assembly.contigs.remove('c1')
        input.membership.B = input.membership.B.collectEntries { p, aa ->
            aa.each { it.proteinId = 'new-' + p }
            [('new-' + p): aa]
        }
        input.claims.findAll { it.genomeId == 'B' }.each { it.proteinId = 'new-' + it.proteinId }
        Collections.reverse(input.claims)

        expect:
        identities(input).A == identities(input).B
    }

    def "distinct observed inputs remain distinct: #field"() {
        given:
        def input = pair()
        change(input)

        expect:
        identities(input).A != identities(input).B

        where:
        field | change
        'confidence' | { i -> i.claims.find { it.genomeId == 'B' }.baselineConfidence = 0.8d }
        'provenance' | { i -> i.claims.find { it.genomeId == 'B' }.provenance = 'independent-source' }
        'EC mapping' | { i -> i.claims.find { it.genomeId == 'B' }.ecNumber = '1.1.1.1' }
        'EC direct' | { i -> i.claims.find { it.genomeId == 'B' }.ecDirect = true }
        'EC multiplicity' | { i -> i.claims.find { it.genomeId == 'B' }.ecReactionCount = 2 }
        'membership weight' | { i -> i.membership.B.B_20[0].weight = 0.5d }
        'subfamily' | { i -> i.membership.B.B_20[0].subfamilyId = 'sub' }
        'paralog' | { i -> i.assign('B', new CohortInput.FamilyAssignment(proteinId: 'B_0', orthogroupId: 'OG1')) }
        'coordinate' | { i -> i.genomes.B.layout.get('B_0').end = 950 }
        'strand' | { i -> i.genomes.B.layout.get('B_0').strand = Strand.MINUS }
        'assembly status' | { i -> i.genomes.B.assembly.status = AssemblyStatus.DRAFT }
        'assembly completeness' | { i -> i.genomes.B.assembly.completenessFraction = 0.9d }
        'contig length' | { i -> i.genomes.B.assembly.contigs.c1.length = 50000 }
        'contig closure' | { i -> i.genomes.B.assembly.contigs.c1.closed = null }
        'search certification' | { i -> i.genomes.B.familySearchComplete = null }
        'search coverage' | { i -> i.genomes.B.familySearchCoverage = 1d }
        'known empty organism reactions' | { i -> i.setOrganismReactions('B', []) }
        'role evidence' | { i -> i.addRoleEvidence(new CohortInput.RoleEvidence(genomeId: 'B', proteinId: 'B_20', source: 'domain', scores: [(Role.CATALYTIC): 0.9d])) }
        'unplaced protein' | { i -> i.genomes.B.unplacedProteins.add('unplaced') }
    }

    def "equal evidence is collapsed only inside a strain unit and subsets keep their representative"() {
        given:
        def input = pair()
        addGenome(input, [id: 'C', clade: 'c', cluster: 'other', status: 'COMPLETE', closed: true,
                          family: 'OG1', candidate: 'T', anchors: ['N1', 'N2']])
        def clades = new CladeModel(input.genomes.values(), null, 0.1d, 1d)
        def hashes = identities(input)

        expect:
        hashes.values().toSet().size() == 1
        clades.units(['A', 'B', 'C'], hashes)*.genomes.flatten().size() == 2
        clades.units(['B'], hashes)[0].genomes == ['B']
        clades.units(['A', 'B', 'C'])*.genomes.flatten().size() == 3
    }
}
