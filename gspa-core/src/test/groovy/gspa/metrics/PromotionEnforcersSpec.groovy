package gspa.metrics

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import gspa.ontology.GoOntology
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * CompletenessEnforcer + CoherenceEnforcer promotion/demotion, over a tiny
 * inline ontology (GO:0000010 is_a GO:0000001).
 */
class PromotionEnforcersSpec extends Specification {

    @TempDir
    Path tmp

    @Shared
    GoOntology go

    def setupSpec() {
        // loaded per-spec in setup() because @TempDir is per-iteration
    }

    private GoOntology ontology() {
        if (go != null) return go
        def obo = tmp.resolve('mini.obo').toFile()
        obo.text = '''format-version: 1.2

[Term]
id: GO:0000001
name: root

[Term]
id: GO:0000010
name: child
is_a: GO:0000001
'''
        go = new GoOntology()
        go.loadOwl(obo)
        go
    }

    private Protein protein(String id, Map<String, Double> terms) {
        def p = new Protein(id: id, sequence: 'M')
        terms.each { v, s ->
            p.annotations.add(new Annotation(type: AnnotationType.GO, value: v, score: s, source: 'pred', evidence: 'IEA'))
        }
        p
    }

    private Genome genomeOf(Protein... proteins) {
        def g = new Genome(id: 'g')
        def c = new Contig(id: 'c1')
        proteins.each { c.addProtein(it) }
        g.addContig(c)
        g
    }

    def 'completeness promotes a missing essential onto the protein with the best near-ancestor'() {
        given: 'GO:0000010 is the missing essential; a protein carries its parent GO:0000001'
        def p = protein('p1', ['GO:0000001': 0.8d])
        def genome = genomeOf(p)
        def report = new EnforcementReport(provenance: true)
        def enf = new CompletenessEnforcer(goOntology: ontology(), report: report, promoteMinScore: 0.01d)

        when:
        def r = enf.enforce(genome, ['GO:0000010'])

        then: 'the essential is imputed onto p1, with provenance + ISC evidence'
        r.promoted == 1
        def added = p.annotations.goAnnotations().find { it.value == 'GO:0000010' }
        added != null
        added.evidence == 'ISC'
        added.source == 'gspa-enforce:completeness'
        added.provenance.any { it.startsWith('completeness:promote') }
        report.actions.size() == 1
    }

    def 'completeness leaves a gap unfilled when there is no evidence'() {
        given: 'no protein carries anything near GO:0000010'
        def p = protein('p1', ['GO:0000099': 0.9d])
        def enf = new CompletenessEnforcer(goOntology: ontology(), promoteMinScore: 0.01d)

        when:
        def r = enf.enforce(genomeOf(p), ['GO:0000010'])

        then:
        r.promoted == 0
        r.unfillable == 1
        p.annotations.goAnnotations().every { it.value != 'GO:0000010' }
    }

    def 'coherence demotes an obligate-complex singleton (demote-only)'() {
        given:
        def p = protein('p1', ['GO:0000010': 0.9d])
        def enf = new CoherenceEnforcer(goOntology: ontology(), promotePartner: false,
            enforceProcess: false, complexClassification: ['GO:0000010': 'n'])

        when:
        def r = enf.enforce(genomeOf(p))

        then:
        r.complexDemoted == 1
        p.annotations.goAnnotations().every { it.value != 'GO:0000010' }
    }

    def 'coherence promotes a complex partner when one is plausible'() {
        given: 'p1 carries the singleton complex; p2 carries its parent'
        def p1 = protein('p1', ['GO:0000010': 0.9d])
        def p2 = protein('p2', ['GO:0000001': 0.8d])
        def report = new EnforcementReport(provenance: true)
        def enf = new CoherenceEnforcer(goOntology: ontology(), report: report, promotePartner: true,
            enforceProcess: false, promoteMinScore: 0.01d, complexClassification: ['GO:0000010': 'n'])

        when:
        def r = enf.enforce(genomeOf(p1, p2))

        then: 'the complex term is promoted onto p2 (the partner), p1 keeps it'
        r.complexPromoted == 1
        r.complexDemoted == 0
        p2.annotations.goAnnotations().any { it.value == 'GO:0000010' }
        p1.annotations.goAnnotations().any { it.value == 'GO:0000010' }
        report.actions.any { it.dimension == 'coherence' && it.action == 'promote' }
    }

    def 'coherence promotes a missing has_part process partner'() {
        given: 'process pair (C present, F=GO:0000010 missing); a protein carries the parent of F'
        def p = protein('p1', ['GO:0000001': 0.7d])
        def enf = new CoherenceEnforcer(goOntology: ontology(), promotePartner: false,
            enforceProcess: true, promoteMinScore: 0.01d)
        def pairs = [new AbstractMap.SimpleEntry('GO:0000050', 'GO:0000010')] as List<Map.Entry<String, String>>

        when:
        def r = enf.enforce(genomeOf(p), pairs)

        then:
        r.processPromoted == 1
        p.annotations.goAnnotations().any { it.value == 'GO:0000010' }
    }
}
