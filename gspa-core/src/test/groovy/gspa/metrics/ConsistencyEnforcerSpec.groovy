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

    def 'enforcer and SymbolicEnforcementStage flag and remove forbidden EC claims for bacterial strain via lineage while keeping unmapped EC as not checkable'() {
        given: 'a bacterial strain taxon (NCBITaxon_640081) resolved via lineage under NCBITaxon_2'
        def tc = new TaxonConstraints()
        tc.neverInTaxon['GO:0004164'] << 'NCBITaxon_2'   // diphthine synthase (EC:2.1.1.98)
        tc.onlyInTaxon['GO:0046905'] << 'NCBITaxon_2759' // phytoene synthase (EC:2.5.1.32)
        def hier = tmp.resolve('h_ec.tsv').toFile()
        hier.text = "Term\tRelationship\tParent\n" +
            "NCBITaxon_2\tis_a\tNCBITaxon_131567\n" +
            "NCBITaxon_2\tdisjoint_from\tNCBITaxon_2759\n"
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(hier)
        checker.ecToGoAlternatives = [
            'EC:2.1.1.98': ['GO:0004164'],
            'EC:2.5.1.32': ['GO:0046905'],
            'EC:2.7.1.1' : ['GO:0004396'],
        ]
        checker.organismTaxon = 'NCBITaxon_640081'
        checker.applyTaxonLineage('NCBITaxon_131567,NCBITaxon_2,NCBITaxon_1224', 'NCBITaxon_640081')

        def genome = new Genome(id: 'g')
        def c = new Contig(id: 'c1')
        def p = new Protein(id: 'p1', sequence: 'M')
        p.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:2.1.1.98', score: 0.91d, source: 'topec', evidence: 'IEA'))
        p.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:2.5.1.32', score: 0.85d, source: 'topec', evidence: 'IEA'))
        p.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:2.7.1.1', score: 0.95d, source: 'topec', evidence: 'IEA'))
        p.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:1.97.1.12', score: 0.80d, source: 'topec', evidence: 'IEA')) // unmapped
        c.addProtein(p)
        genome.addContig(c)

        def report = new EnforcementReport(provenance: true)
        def enforcer = new ConsistencyEnforcer(satChecker: checker, ecToGoAlternatives: checker.ecToGoAlternatives,
            mode: 'relax', report: report)

        when:
        def result = enforcer.enforce(genome)

        then: 'forbidden ECs are removed with GO cause in provenance; allowed and unmapped ECs survive; unmapped EC is counted as not checkable'
        result.status == 'applied'
        result.taxonInHierarchy
        result.constraintsApplied
        result.ecClaimsChecked == 3
        result.ecClaimsNotCheckable == 1
        result.notCheckableEcTerms == ['EC:1.97.1.12'] as Set
        result.ecViolations == 2
        result.annotationsAffected == 2
        p.annotations.ecAnnotations()*.value == ['EC:2.7.1.1', 'EC:1.97.1.12']
        report.actions.size() == 2
        report.actions.find { it.term == 'EC:2.1.1.98' }.reason.contains('GO:0004164')
        report.actions.find { it.term == 'EC:2.5.1.32' }.reason.contains('GO:0046905')
    }

    def 'enforcer reports explicit taxon not in hierarchy: constraints not applied when strain has no lineage'() {
        given: 'NCBITaxon_640081 not in hierarchy and no lineage supplied'
        def tc = new TaxonConstraints()
        tc.neverInTaxon['GO:0004164'] << 'NCBITaxon_2'
        def hier = tmp.resolve('h_nolin.tsv').toFile()
        hier.text = "Term\tRelationship\tParent\nNCBITaxon_2\tis_a\tNCBITaxon_131567\n"
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(hier)
        checker.ecToGoAlternatives = ['EC:2.1.1.98': ['GO:0004164']]
        checker.organismTaxon = 'NCBITaxon_640081'

        def genome = new Genome(id: 'g')
        def c = new Contig(id: 'c1')
        def p = new Protein(id: 'p1', sequence: 'M')
        p.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:2.1.1.98', score: 0.91d, source: 'topec', evidence: 'IEA'))
        p.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:1.97.1.12', score: 0.80d, source: 'topec', evidence: 'IEA'))
        c.addProtein(p)
        genome.addContig(c)

        def enforcer = new ConsistencyEnforcer(satChecker: checker, ecToGoAlternatives: checker.ecToGoAlternatives, mode: 'relax')

        when:
        def result = enforcer.enforce(genome)

        then:
        result.status == 'taxon not in hierarchy: constraints not applied'
        !result.taxonInHierarchy
        !result.constraintsApplied
        result.violations == 0
        result.annotationsAffected == 0
        result.ecClaimsNotCheckable == 1
        p.annotations.ecAnnotations().size() == 2
    }

    def 'IntegrationWorkflow and SymbolicEnforcementStage remove forbidden EC claims with strain lineage and record explicit ledger state without lineage'() {
        given:
        File ec2goFile = tmp.resolve('ec2go.txt').toFile()
        ec2goFile.text = "EC:2.1.1.98 > GO:diphthine synthase activity ; GO:0004164\n" +
            "EC:2.5.1.32 > GO:phytoene synthase activity ; GO:0046905\n" +
            "EC:2.7.1.1 > GO:hexokinase activity ; GO:0004396\n"
        def claims = [
            new gspa.integration.EvidenceClaim(proteinId: 'protA', functionId: 'EC:2.1.1.98',
                functionType: AnnotationType.EC, source: 'topec', evidenceType: gspa.integration.EvidenceType.SEQUENCE_DEEPLEARNING,
                rawScore: 0.92d, calibratedProb: 0.92d),
            new gspa.integration.EvidenceClaim(proteinId: 'protB', functionId: 'EC:2.5.1.32',
                functionType: AnnotationType.EC, source: 'topec', evidenceType: gspa.integration.EvidenceType.SEQUENCE_DEEPLEARNING,
                rawScore: 0.88d, calibratedProb: 0.88d),
            new gspa.integration.EvidenceClaim(proteinId: 'protC', functionId: 'EC:2.7.1.1',
                functionType: AnnotationType.EC, source: 'topec', evidenceType: gspa.integration.EvidenceType.SEQUENCE_DEEPLEARNING,
                rawScore: 0.95d, calibratedProb: 0.95d),
            new gspa.integration.EvidenceClaim(proteinId: 'protD', functionId: 'EC:1.97.1.12',
                functionType: AnnotationType.EC, source: 'topec', evidenceType: gspa.integration.EvidenceType.SEQUENCE_DEEPLEARNING,
                rawScore: 0.85d, calibratedProb: 0.85d),
        ]

        when: 'run without lineage for unknown strain NCBITaxon_640081'
        def cfgNoLin = new gspa.integration.workflow.IntegrationWorkflowConfig(
            lite: true,
            enabledPriors: 'consistency',
            ec2goFile: ec2goFile,
            organismTaxon: 'NCBITaxon_640081',
            enforceStage2: true,
            enforceConsistencyMode: 'relax'
        )
        def resNoLin = new gspa.integration.workflow.IntegrationWorkflow(cfgNoLin).run(claims)
        def stagesNoLin = resNoLin.audit.report().stages.collectEntries { [it.name, it] }

        and: 'run with lineage for NCBITaxon_640081'
        def cfgWithLin = new gspa.integration.workflow.IntegrationWorkflowConfig(
            lite: true,
            enabledPriors: 'consistency',
            ec2goFile: ec2goFile,
            organismTaxon: 'NCBITaxon_640081',
            taxonLineage: 'NCBITaxon_131567,NCBITaxon_2,NCBITaxon_1224,NCBITaxon_28216',
            enforceStage2: true,
            enforceConsistencyMode: 'relax'
        )
        def resWithLin = new gspa.integration.workflow.IntegrationWorkflow(cfgWithLin).run(claims)
        def stagesWithLin = resWithLin.audit.report().stages.collectEntries { [it.name, it] }

        then: 'without lineage, constraints are explicitly recorded as not applied and all 4 ECs remain'
        stagesNoLin.reference_preparation.taxon_in_hierarchy == false
        stagesNoLin.reference_preparation.taxon_consistency_state == 'taxon not in hierarchy: constraints not applied'
        stagesNoLin.reference_preparation.taxon_not_in_hierarchy_count == 1
        stagesNoLin.symbolic_enforcement.taxon_consistency_state == 'taxon not in hierarchy: constraints not applied'
        stagesNoLin.symbolic_enforcement.constraints_applied == false
        resNoLin.integrated.provenance.keySet().size() == 4

        and: 'with lineage, constraints are applied, the 2 forbidden ECs are removed, and the allowed + unmapped ECs remain'
        stagesWithLin.reference_preparation.taxon_in_hierarchy == true
        stagesWithLin.reference_preparation.taxon_consistency_state == 'applied'
        stagesWithLin.reference_preparation.ec_claims_checkable == 3
        stagesWithLin.reference_preparation.ec_claims_not_checkable == 1
        stagesWithLin.symbolic_enforcement.taxon_consistency_state == 'applied'
        stagesWithLin.symbolic_enforcement.constraints_applied == true
        stagesWithLin.symbolic_enforcement.ec_violations == 2
        stagesWithLin.symbolic_enforcement.ec_claims_checked == 3
        stagesWithLin.symbolic_enforcement.ec_claims_not_checkable == 1
        resWithLin.integrated.provenance.keySet() == ['protC|EC|EC:2.7.1.1', 'protD|EC|EC:1.97.1.12'] as Set
    }
}
