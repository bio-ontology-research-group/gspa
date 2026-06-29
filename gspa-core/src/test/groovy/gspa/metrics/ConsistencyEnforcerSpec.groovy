package gspa.metrics

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * ConsistencyEnforcer: acts on taxon-inconsistent GO annotations as a
 * post-annotation pass. Uses a minimal hand-built taxon constraint set so the
 * test needs no GO ontology: GO:0000001 is only_in T1, GO:0000002 only_in T2,
 * and T1/T2 are disjoint sibling taxa — so a protein carrying both is UNSAT.
 */
class ConsistencyEnforcerSpec extends Specification {

    private SatConsistencyChecker conflictingChecker() {
        def tc = new TaxonConstraints()
        tc.onlyInTaxon['GO:0000001'] << 'NCBITaxon:T1'
        tc.onlyInTaxon['GO:0000002'] << 'NCBITaxon:T2'
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyHierarchy([
            'NCBITaxon:T1': 'NCBITaxon:Root',
            'NCBITaxon:T2': 'NCBITaxon:Root',
        ])
        checker
    }

    private Genome genomeWith(List<String> terms) {
        def genome = new Genome(id: 'g')
        def contig = new Contig(id: 'c1')
        def protein = new Protein(id: 'p1', sequence: 'M')
        terms.each { t ->
            protein.annotations.add(new Annotation(type: AnnotationType.GO, value: t,
                score: 0.9d, source: 'test', evidence: 'IEA'))
        }
        contig.addProtein(protein)
        genome.addContig(contig)
        genome
    }

    def 'sanity: the two terms really are SAT-inconsistent together'() {
        expect:
        !conflictingChecker().check(['GO:0000001', 'GO:0000002'] as Set).consistent
        conflictingChecker().check(['GO:0000001'] as Set).consistent
    }

    def 'remove mode drops the violating annotations'() {
        given:
        def genome = genomeWith(['GO:0000001', 'GO:0000002'])
        def enforcer = new ConsistencyEnforcer(satChecker: conflictingChecker(), mode: 'remove')

        when:
        def result = enforcer.enforce(genome)

        then:
        result.violations >= 1
        result.proteinsAffected == 1
        result.annotationsAffected == 2
        genome.proteins[0].annotations.goAnnotations().isEmpty()
    }

    def 'downrank mode scales the score and flags evidence'() {
        given:
        def genome = genomeWith(['GO:0000001', 'GO:0000002'])
        def enforcer = new ConsistencyEnforcer(satChecker: conflictingChecker(),
            mode: 'downrank', downrankFactor: 0.5d)

        when:
        enforcer.enforce(genome)
        def anns = genome.proteins[0].annotations.goAnnotations()

        then:
        anns.size() == 2
        anns.every { Math.abs(it.score - 0.45d) < 1e-9 }
        anns.every { it.evidence.contains('taxon-inconsistent') }
    }

    def 'flag mode keeps annotations but marks them'() {
        given:
        def genome = genomeWith(['GO:0000001', 'GO:0000002'])
        def enforcer = new ConsistencyEnforcer(satChecker: conflictingChecker(), mode: 'flag')

        when:
        enforcer.enforce(genome)
        def anns = genome.proteins[0].annotations.goAnnotations()

        then:
        anns.size() == 2
        anns.every { it.evidence.contains('taxon-inconsistent') }
    }

    def 'a consistent protein is left untouched'() {
        given:
        def genome = genomeWith(['GO:0000001'])
        def enforcer = new ConsistencyEnforcer(satChecker: conflictingChecker(), mode: 'remove')

        when:
        def result = enforcer.enforce(genome)

        then:
        result.annotationsAffected == 0
        genome.proteins[0].annotations.goAnnotations().size() == 1
    }

    @TempDir
    Path tmp

    def 'minimal-flip resolves a co-annotation conflict by dropping the lower-weight term'() {
        given: 'two mutually-disjoint terms on one protein, different scores, no organism'
        def tc = new TaxonConstraints()
        tc.onlyInTaxon['GO:euk'] << 'NCBITaxon_2759'
        tc.onlyInTaxon['GO:bact'] << 'NCBITaxon_2'
        def hier = tmp.resolve('h2.tsv').toFile()
        hier.text = "Term\tRelationship\tParent\nNCBITaxon_2\tdisjoint_from\tNCBITaxon_2759\n"
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(hier)

        def genome = new Genome(id: 'g')
        def c = new Contig(id: 'c1')
        def p = new Protein(id: 'p1', sequence: 'M')
        p.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:euk', score: 0.9d, source: 'x', evidence: 'IEA'))
        p.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:bact', score: 0.1d, source: 'x', evidence: 'IEA'))
        c.addProtein(p)
        genome.addContig(c)
        def enforcer = new ConsistencyEnforcer(satChecker: checker, mode: 'minimal-flip')

        when:
        def result = enforcer.enforce(genome)

        then: 'only the lower-weight (bact, 0.1) term is dropped, not both'
        result.annotationsAffected == 1
        p.annotations.goAnnotations()*.value == ['GO:euk']
    }

    def 'enforcement records provenance actions and trail lines'() {
        given:
        def report = new EnforcementReport(provenance: true)
        def genome = genomeWith(['GO:0000001', 'GO:0000002'])
        def enforcer = new ConsistencyEnforcer(satChecker: conflictingChecker(), mode: 'flag', report: report)

        when:
        enforcer.enforce(genome)

        then:
        report.actions.size() == 2
        report.actions.every { it.dimension == 'consistency' && it.action == 'flag' }
        report.actions.every { it.basis.startsWith('test@') }
        genome.proteins[0].annotations.goAnnotations().every { it.provenance.any { line -> line.startsWith('consistency:flag') } }
    }

    def 'organism mode removes only the conflicting term, not a consistent co-annotation'() {
        given: 'organism = Bacteria; one eukaryote-only term, one cellular-organism term'
        def tc = new TaxonConstraints()
        tc.onlyInTaxon['GO:euk'] << 'NCBITaxon_2759'
        tc.onlyInTaxon['GO:cell'] << 'NCBITaxon_131567'
        def hier = tmp.resolve('h.tsv').toFile()
        hier.text = "Term\tRelationship\tParent\n" +
            "NCBITaxon_2\tis_a\tNCBITaxon_131567\n" +
            "NCBITaxon_2\tdisjoint_from\tNCBITaxon_2759\n"
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(hier)
        checker.organismTaxon = 'NCBITaxon_2'

        def genome = genomeWith(['GO:euk', 'GO:cell'])
        def enforcer = new ConsistencyEnforcer(satChecker: checker, mode: 'remove')

        when:
        def result = enforcer.enforce(genome)
        def remaining = genome.proteins[0].annotations.goAnnotations()*.value

        then: 'the eukaryote-only term is removed; the cellular-organism term survives'
        result.annotationsAffected == 1
        remaining == ['GO:cell']
    }
}
