package gspa.metrics

import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import spock.lang.Shared
import spock.lang.Specification

/**
 * Locks {@link TaxonInference}'s fast direct-hierarchy violation count to the
 * authoritative {@link SatConsistencyChecker}: for every sampled (candidate
 * taxon, constrained term) pair, the cheap {@code violates()} predicate must
 * agree with a full SAT solve asserting that organism. Uses the real bundled
 * constraint + hierarchy resources so the equivalence holds on production data.
 */
class TaxonInferenceSatEquivalenceSpec extends Specification {

    @Shared SatConsistencyChecker checker
    @Shared TaxonInference inference
    @Shared List<String> sampleTerms

    def setupSpec() {
        def tc = new TaxonConstraints()
        tc.loadFromTsv(resource('/taxon-constraints/go-taxon-constraints.tsv'))
        checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(resource('/taxon-constraints/ncbi-taxon-hierarchy.tsv'))
        inference = new TaxonInference(checker: checker)
        // A spread of constrained terms: every Nth only_in and never_in term.
        def only = tc.onlyInTaxon.keySet().toList()
        def never = tc.neverInTaxon.keySet().toList()
        sampleTerms = (stride(only, 90) + stride(never, 90)).unique()
    }

    private static File resource(String path) {
        new File(TaxonInferenceSatEquivalenceSpec.getResource(path).toURI())
    }

    private static List<String> stride(List<String> xs, int n) {
        def step = Math.max(1, (int) (xs.size() / n))
        (0..<xs.size()).step(step).collect { xs[it] }
    }

    def "direct violation count agrees with SAT for #organism"() {
        expect:
        sampleTerms.every { term ->
            checker.organismTaxon = organism
            boolean satConsistent = checker.check([term] as Set).consistent
            boolean directConsistent = !invokeViolates(organism, term)
            assert directConsistent == satConsistent,
                "mismatch organism=${organism} term=${term}: sat=${satConsistent} direct=${directConsistent}"
            true
        }

        cleanup:
        checker.organismTaxon = null

        where:
        organism << [
            'NCBITaxon_2',      // Bacteria
            'NCBITaxon_2157',   // Archaea
            'NCBITaxon_2759',   // Eukaryota
            'NCBITaxon_10239',  // Viruses
            'NCBITaxon_33208',  // Metazoa
            'NCBITaxon_7742',   // Vertebrata
            'NCBITaxon_40674',  // Mammalia
            'NCBITaxon_9606',   // Homo sapiens
            'NCBITaxon_50557',  // Insecta
            'NCBITaxon_4751',   // Fungi
            'NCBITaxon_4932',   // Saccharomyces cerevisiae
            'NCBITaxon_33090',  // Viridiplantae
            'NCBITaxon_131567', // cellular organisms
        ]
    }

    /** {@code violates} is private; Groovy lets us invoke it directly for the lock. */
    private boolean invokeViolates(String organism, String term) {
        inference.invokeMethod('violates', [organism, term] as Object[]) as boolean
    }
}
