package gspa.ontology

import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class SatConsistencyCheckerSpec extends Specification {

    @TempDir
    Path tmp

    private File hierarchyFile(String body) {
        def f = tmp.resolve('hier.tsv').toFile()
        f.text = "Term\tRelationship\tParent\n" + body
        f
    }

    private File resource(String path) {
        def f = tmp.resolve(path.replaceAll('/', '_')).toFile()
        f.withOutputStream { out -> getClass().getResourceAsStream(path).withCloseable { it.transferTo(out) } }
        f
    }

    def "should detect consistent annotations"() {
        given: "A simple taxonomy: Bacteria -> Proteobacteria -> Gammaproteobacteria -> E.coli"
        def constraints = new TaxonConstraints()
        constraints.onlyInTaxon['GO:0006412'] = ['NCBITaxon:131567'] as Set // cellular organisms
        constraints.onlyInTaxon['GO:0006260'] = ['NCBITaxon:131567'] as Set

        def checker = new SatConsistencyChecker(constraints)
        checker.loadTaxonomyHierarchy([
            'NCBITaxon:562'   : 'NCBITaxon:561',      // E.coli -> Escherichia
            'NCBITaxon:561'   : 'NCBITaxon:543',      // Escherichia -> Enterobacteriaceae
            'NCBITaxon:543'   : 'NCBITaxon:91347',    // -> Enterobacterales
            'NCBITaxon:91347' : 'NCBITaxon:1236',     // -> Gammaproteobacteria
            'NCBITaxon:1236'  : 'NCBITaxon:1224',     // -> Proteobacteria
            'NCBITaxon:1224'  : 'NCBITaxon:2',        // -> Bacteria
            'NCBITaxon:2'     : 'NCBITaxon:131567',   // -> cellular organisms
        ])

        when:
        def result = checker.check(['GO:0006412', 'GO:0006260'] as Set)

        then:
        result.consistent
        result.violations.isEmpty()
    }

    def "should detect inconsistent annotations (positive + negative on same taxon)"() {
        given: "One GO term requires Bacteria, another says never in Bacteria"
        def constraints = new TaxonConstraints()
        constraints.onlyInTaxon['GO:0009306'] = ['NCBITaxon:2'] as Set    // must be in Bacteria
        constraints.neverInTaxon['GO:0015979'] = ['NCBITaxon:2'] as Set   // must NOT be in Bacteria

        def checker = new SatConsistencyChecker(constraints)
        checker.loadTaxonomyHierarchy([
            'NCBITaxon:2' : 'NCBITaxon:131567',
        ])

        when:
        def result = checker.check(['GO:0009306', 'GO:0015979'] as Set)

        then:
        !result.consistent
        result.violations.size() >= 1
        result.violations[0].type == gspa.model.ConsistencyViolation.ViolationType.TAXON_CONFLICT
    }

    def "should detect inconsistent annotations via sibling disjointness"() {
        given: "Annotations requiring both Bacteria-only and Plant-only terms (siblings under cellular organisms)"
        def constraints = new TaxonConstraints()
        constraints.onlyInTaxon['GO:0009306'] = ['NCBITaxon:2'] as Set      // only in Bacteria
        constraints.onlyInTaxon['GO:0015979'] = ['NCBITaxon:33090'] as Set   // only in Viridiplantae

        def checker = new SatConsistencyChecker(constraints)
        // Bacteria and Viridiplantae are siblings under cellular organisms
        checker.loadTaxonomyHierarchy([
            'NCBITaxon:2'     : 'NCBITaxon:131567',   // Bacteria -> cellular organisms
            'NCBITaxon:33090' : 'NCBITaxon:131567',   // Plants -> cellular organisms
            'NCBITaxon:2759'  : 'NCBITaxon:131567',   // Eukaryota -> cellular organisms
        ])

        when:
        def result = checker.check(['GO:0009306', 'GO:0015979'] as Set)

        then: "Should be inconsistent because an organism cannot be both Bacteria and Plant"
        !result.consistent
        result.violations.size() >= 1
    }

    def "should allow consistent annotations within same lineage"() {
        given: "Annotations requiring Bacteria and Proteobacteria (child of Bacteria)"
        def constraints = new TaxonConstraints()
        constraints.onlyInTaxon['GO:0009306'] = ['NCBITaxon:2'] as Set      // only in Bacteria
        constraints.onlyInTaxon['GO:0042597'] = ['NCBITaxon:1224'] as Set    // only in Proteobacteria

        def checker = new SatConsistencyChecker(constraints)
        checker.loadTaxonomyHierarchy([
            'NCBITaxon:1224'  : 'NCBITaxon:2',        // Proteobacteria -> Bacteria
            'NCBITaxon:2'     : 'NCBITaxon:131567',   // Bacteria -> cellular organisms
        ])

        when:
        def result = checker.check(['GO:0009306', 'GO:0042597'] as Set)

        then: "Should be consistent - Proteobacteria is within Bacteria"
        result.consistent
    }

    def "should detect never_in_taxon violations"() {
        given: "A bacterial function annotated to a genome with a plant-only term"
        def constraints = new TaxonConstraints()
        constraints.onlyInTaxon['GO:0006412'] = ['NCBITaxon:2'] as Set     // only in Bacteria
        constraints.neverInTaxon['GO:0006412_fake'] = ['NCBITaxon:2'] as Set // never in Bacteria (fake conflict)
        // force both terms to require and exclude Bacteria
        constraints.onlyInTaxon['GO:0006412_fake'] = ['NCBITaxon:2'] as Set

        def checker = new SatConsistencyChecker(constraints)
        checker.loadTaxonomyHierarchy([:])

        when:
        def result = checker.check(['GO:0006412', 'GO:0006412_fake'] as Set)

        then:
        !result.consistent
    }

    def "should handle empty annotations"() {
        given:
        def constraints = new TaxonConstraints()
        def checker = new SatConsistencyChecker(constraints)

        when:
        def result = checker.check([] as Set)

        then:
        result.consistent
    }

    // --- explicit disjointness + organism assertion (Asaad model) ------------

    def 'explicit disjoint_from makes two only_in taxa conflict'() {
        given:
        def tc = new TaxonConstraints()
        tc.onlyInTaxon['GO:0000010'] << 'NCBITaxon_2'      // Bacteria
        tc.onlyInTaxon['GO:0000020'] << 'NCBITaxon_2759'   // Eukaryota
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(hierarchyFile("NCBITaxon_2\tdisjoint_from\tNCBITaxon_2759\n"))

        expect:
        !checker.check(['GO:0000010', 'GO:0000020'] as Set).consistent
        checker.check(['GO:0000010'] as Set).consistent
    }

    def 'asserting the organism taxon flags a lone term it cannot carry'() {
        given: 'organism is a bacterium; one term is eukaryote-only, another bacteria-only'
        def tc = new TaxonConstraints()
        tc.onlyInTaxon['GO:euk'] << 'NCBITaxon_2759'
        tc.onlyInTaxon['GO:bact'] << 'NCBITaxon_2'
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(hierarchyFile("NCBITaxon_2\tdisjoint_from\tNCBITaxon_2759\n"))
        checker.organismTaxon = 'NCBITaxon_2'

        expect:
        !checker.check(['GO:euk'] as Set).consistent
        checker.check(['GO:bact'] as Set).consistent
    }

    def 'never_in the organism taxon is a violation'() {
        given:
        def tc = new TaxonConstraints()
        tc.neverInTaxon['GO:nb'] << 'NCBITaxon_2'
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(hierarchyFile("NCBITaxon_2\tis_a\tNCBITaxon_131567\n"))
        checker.organismTaxon = 'NCBITaxon_2'

        expect:
        !checker.check(['GO:nb'] as Set).consistent
    }

    def 'union_of members are treated as pairwise disjoint'() {
        given: 'a disjoint-union node groups Bacteria and Archaea'
        def tc = new TaxonConstraints()
        tc.onlyInTaxon['GO:b'] << 'NCBITaxon_2'
        tc.onlyInTaxon['GO:a'] << 'NCBITaxon_2157'
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(hierarchyFile(
            "NCBITaxon_U\tunion_of\tNCBITaxon_2\n" +
            "NCBITaxon_U\tunion_of\tNCBITaxon_2157\n"))

        expect:
        !checker.check(['GO:b', 'GO:a'] as Set).consistent
    }

    def 'bundled GO taxon constraints + backbone flag a eukaryote-only term on a bacterium'() {
        given: 'the real vendored constraint + hierarchy data, organism = Bacteria'
        def tc = new TaxonConstraints()
        tc.loadFromTsv(resource('/taxon-constraints/go-taxon-constraints.tsv'))
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(resource('/taxon-constraints/ncbi-taxon-hierarchy.tsv'))
        checker.organismTaxon = 'NCBITaxon_2'

        expect: 'GO:0000001 is only_in Eukaryota (2759) in the bundled data'
        tc.onlyInTaxon['GO:0000001'].contains('NCBITaxon_2759')

        and: 'so it is inconsistent on a bacterium, while a cellular-organism term is fine'
        !checker.check(['GO:0000001'] as Set).consistent
    }
}
