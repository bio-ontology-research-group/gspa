package gspa.ontology

import spock.lang.Specification

class SatConsistencyCheckerSpec extends Specification {

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
}
