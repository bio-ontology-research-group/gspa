package gspa.ontology

import gspa.integration.CalibrationTable
import gspa.integration.IntegrationState
import gspa.integration.prior.AnatomicalPrior
import gspa.integration.prior.CoherencePrior
import gspa.metrics.ConsistencyEnforcer
import gspa.metrics.EnforcementReport
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import spock.lang.Shared
import spock.lang.Specification

class AnatomicalAndComplexConstraintsSpec extends Specification {

    @Shared AnatomicalConstraints anatomical
    @Shared TaxonConstraints taxonConstraints
    @Shared SatConsistencyChecker checker

    def setupSpec() {
        anatomical = AnatomicalConstraints.loadBundled()
        taxonConstraints = new TaxonConstraints()
        File tcFile = new File(getClass().getResource('/taxon-constraints/go-taxon-constraints.tsv').toURI())
        taxonConstraints.loadFromTsv(tcFile)
        anatomical.augmentTaxonConstraints(taxonConstraints)
        checker = new SatConsistencyChecker(taxonConstraints)
        File hierFile = new File(getClass().getResource('/taxon-constraints/ncbi-taxon-hierarchy.tsv').toURI())
        checker.loadTaxonomyTsv(hierFile)
    }

    def "bundled anatomical constraints load GO-Uberon-PO-CL links and taxon rules"() {
        expect:
        anatomical.constrainedGoTermCount() >= 2600
        anatomical.constrainedAnatomyCount() >= 1200
    }

    def "anatomical constraints accurately distinguish clade-specific developmental processes across eukaryotes"() {
        expect: "heart development (GO:0007507 -> UBERON:0000948) is forbidden in yeast, plants, and bacteria, but valid in fly and mouse"
        !anatomical.checkTerm('GO:0007507', 'NCBITaxon_4932', checker).isEmpty()   // S. cerevisiae
        !anatomical.checkTerm('GO:0007507', 'NCBITaxon_3702', checker).isEmpty()   // A. thaliana
        !anatomical.checkTerm('GO:0007507', 'NCBITaxon_83333', checker).isEmpty()  // E. coli K-12
        anatomical.checkTerm('GO:0007507', 'NCBITaxon_7227', checker).isEmpty()    // D. melanogaster
        anatomical.checkTerm('GO:0007507', 'NCBITaxon_10090', checker).isEmpty()   // M. musculus

        and: "flower development (GO:0009908 -> PO:0009046) is forbidden in yeast, fly, and non-flowering moss P. patens, but valid in A. thaliana"
        !anatomical.checkTerm('GO:0009908', 'NCBITaxon_4932', checker).isEmpty()   // S. cerevisiae
        !anatomical.checkTerm('GO:0009908', 'NCBITaxon_7227', checker).isEmpty()   // D. melanogaster
        !anatomical.checkTerm('GO:0009908', 'NCBITaxon_3218', checker).isEmpty()   // P. patens (Bryophyta)
        anatomical.checkTerm('GO:0009908', 'NCBITaxon_3702', checker).isEmpty()    // A. thaliana (Magnoliopsida)

        and: "neural crest cell migration (GO:0001755 -> UBERON:0002342) is forbidden in invertebrates (C. elegans, D. melanogaster) but valid in vertebrates (M. musculus)"
        !anatomical.checkTerm('GO:0001755', 'NCBITaxon_6239', checker).isEmpty()   // C. elegans
        !anatomical.checkTerm('GO:0001755', 'NCBITaxon_7227', checker).isEmpty()   // D. melanogaster
        anatomical.checkTerm('GO:0001755', 'NCBITaxon_10090', checker).isEmpty()   // M. musculus
    }

    def "AnatomicalPrior penalizes impossible anatomy and boosts active anatomical structures"() {
        given:
        def localChecker = new SatConsistencyChecker(taxonConstraints)
        File hierFile = new File(getClass().getResource('/taxon-constraints/ncbi-taxon-hierarchy.tsv').toURI())
        localChecker.loadTaxonomyTsv(hierFile)
        localChecker.organismTaxon = 'NCBITaxon_7227' // D. melanogaster

        def state = new IntegrationState(new Genome(id: 'dmel'))
        state.anatomicalConstraints = anatomical
        state.satConsistencyChecker = localChecker
        // Two fly proteins already annotated with heart development (GO:0007507 -> UBERON:0000948)
        state.set('P1|GO|GO:0007507', CalibrationTable.probToLogOdds(0.85d))
        state.set('P2|GO|GO:0007507', CalibrationTable.probToLogOdds(0.80d))
        // Third protein has a sub-threshold candidate for heart development, plus a false flower development hit
        state.set('P3|GO|GO:0007507', CalibrationTable.probToLogOdds(0.35d))
        state.set('P3|GO|GO:0009908', CalibrationTable.probToLogOdds(0.70d))

        def prior = new AnatomicalPrior()
        prior.beginIteration(state)

        expect: "flower development (GO:0009908) is penalized on fly, while sub-threshold heart development (GO:0007507) is boosted"
        prior.logOddsBoost('P3', 'P3|GO|GO:0009908', state) == -prior.alphaAnatPenalty
        prior.logOddsBoost('P3', 'P3|GO|GO:0007507', state) == prior.alphaAnatBoost
    }

    def "CoherencePrior boosts missing partner of obligate heteromeric complex and penalizes isolated singleton"() {
        given:
        def state = new IntegrationState(new Genome(id: 'yeast'))
        state.loadBundledComplexTerms()
        // GO:0000124 (SAGA complex) is heteromeric ('n'), GO:0000176 (exosome) is heteromeric ('n')
        // Case 1: GO:0000124 has P1 at 0.85 and candidate partner P2 at 0.35
        state.set('P1|GO|GO:0000124', CalibrationTable.probToLogOdds(0.85d))
        state.set('P2|GO|GO:0000124', CalibrationTable.probToLogOdds(0.35d))
        // Case 2: GO:0000176 has P3 at 0.65 and no other candidate protein in the genome
        state.set('P3|GO|GO:0000176', CalibrationTable.probToLogOdds(0.65d))

        def prior = new CoherencePrior()
        prior.beginIteration(state)

        expect:
        prior.logOddsBoost('P2', 'P2|GO|GO:0000124', state) == prior.alphaCoh * prior.complexPartnerWeight
        prior.logOddsBoost('P3', 'P3|GO|GO:0000176', state) == -prior.alphaCoh * prior.complexSingletonPenaltyWeight
    }

    def "ConsistencyEnforcer in relax mode replaces anatomically invalid GO term with closest valid ancestor"() {
        given: "a mock GoOntology where GO:0009908 (flower development) -> GO:0048856 (anatomical structure development) -> GO:0008150 (root)"
        GoOntology mockGo = GroovyMock(GoOntology) {
            getDirectParents('GO:0009908') >> (['GO:0048856'] as Set)
            getPartOfParents('GO:0009908') >> ([] as Set)
            getDirectParents('GO:0048856') >> (['GO:0008150'] as Set)
            getPartOfParents('GO:0048856') >> ([] as Set)
            getAncestors('GO:0009908') >> (['GO:0048856', 'GO:0008150'] as Set)
            getAncestors('GO:0048856') >> (['GO:0008150'] as Set)
            propagateAnnotations(['GO:0009908'] as Set) >> (['GO:0009908', 'GO:0048856', 'GO:0008150'] as Set)
            propagateAnnotations(['GO:0048856'] as Set) >> (['GO:0048856', 'GO:0008150'] as Set)
            getAspect('GO:0048856') >> 'BP'
            getLabel('GO:0048856') >> 'anatomical structure development'
        }

        def localChecker = new SatConsistencyChecker(taxonConstraints)
        File hierFile = new File(getClass().getResource('/taxon-constraints/ncbi-taxon-hierarchy.tsv').toURI())
        localChecker.loadTaxonomyTsv(hierFile)
        localChecker.organismTaxon = 'NCBITaxon_3218' // Physcomitrium patens (moss, non-flowering plant)

        def genome = new Genome(id: 'ppatens')
        def contig = new Contig(id: 'c1', genome: genome)
        genome.contigs << contig
        def prot = new Protein(id: 'Pp1', contig: contig)
        prot.annotations.add(new Annotation(
            type: AnnotationType.GO,
            value: 'GO:0009908',
            score: 0.82d,
            source: 'predictor',
            evidence: 'IEA',
            goAspect: 'BP'
        ))
        contig.proteins << prot

        def report = new EnforcementReport(provenance: true)
        def enforcer = new ConsistencyEnforcer(
            goOntology: mockGo,
            satChecker: localChecker,
            anatomicalConstraints: anatomical,
            mode: 'relax',
            report: report
        )

        when:
        def res = enforcer.enforce(genome)

        then: "GO:0009908 is relaxed to GO:0048856 with score preserved"
        res.annotationsAffected == 1
        res.annotationsRelaxed == 1
        prot.annotations.goTermIds() == ['GO:0048856'] as Set
        prot.annotations.goAnnotations()[0].score == 0.82d
        prot.annotations.goAnnotations()[0].evidence.contains('relaxed-from:GO:0009908')
    }
}
