package gspa.metrics

import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * TaxonInference: read the organism's domain off the predicted functions via the
 * GO taxon constraints (Asaad-style, in reverse). Minimal hand-built constraints,
 * no GO ontology needed: GO:bact is only_in Bacteria, GO:euk only_in Eukaryota,
 * GO:cell only_in cellular-organisms (shared ancestor), with the four domain
 * taxa mutually disjoint (as in the real NCBI backbone).
 */
class TaxonInferenceSpec extends Specification {

    @TempDir
    Path tmp

    private SatConsistencyChecker domainChecker() {
        def tc = new TaxonConstraints()
        tc.onlyInTaxon['GO:bact'] << 'NCBITaxon_2'
        tc.onlyInTaxon['GO:euk'] << 'NCBITaxon_2759'
        tc.onlyInTaxon['GO:cell'] << 'NCBITaxon_131567'
        def hier = tmp.resolve('h.tsv').toFile()
        hier.text = "Term\tRelationship\tParent\n" +
            "NCBITaxon_2\tis_a\tNCBITaxon_131567\n" +
            "NCBITaxon_2157\tis_a\tNCBITaxon_131567\n" +
            "NCBITaxon_2759\tis_a\tNCBITaxon_131567\n" +
            "NCBITaxon_2\tdisjoint_from\tNCBITaxon_2157\n" +
            "NCBITaxon_2\tdisjoint_from\tNCBITaxon_2759\n" +
            "NCBITaxon_2\tdisjoint_from\tNCBITaxon_10239\n" +
            "NCBITaxon_2157\tdisjoint_from\tNCBITaxon_2759\n" +
            "NCBITaxon_2157\tdisjoint_from\tNCBITaxon_10239\n" +
            "NCBITaxon_2759\tdisjoint_from\tNCBITaxon_10239\n"
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(hier)
        checker
    }

    def 'infers Bacteria from a bacterial term set, with zero violations'() {
        given:
        def inf = new TaxonInference(checker: domainChecker(), goOntology: null, minScore: 0.5d)

        when: 'a bacteria-only term plus a shared cellular-organism term'
        def res = inf.infer(['GO:bact': 0.95d, 'GO:cell': 0.95d])

        then:
        res.taxon == 'NCBITaxon_2'
        res.label == 'Bacteria'
        res.confident
        res.candidates.find { it.taxon == 'NCBITaxon_2' }.forbidden == 0
        res.candidates.find { it.taxon == 'NCBITaxon_2759' }.forbidden >= 1
    }

    def 'infers Eukaryota from a eukaryotic term set'() {
        given:
        def inf = new TaxonInference(checker: domainChecker(), goOntology: null, minScore: 0.5d)

        when:
        def res = inf.infer(['GO:euk': 0.9d, 'GO:cell': 0.9d])

        then:
        res.taxon == 'NCBITaxon_2759'
        res.label == 'Eukaryota'
    }

    def 'ignores the low-score tail: nothing above the threshold is undecidable'() {
        given: 'a high threshold so the only term falls below it'
        def inf = new TaxonInference(checker: domainChecker(), goOntology: null, minScore: 0.9d)

        when:
        def res = inf.infer(['GO:euk': 0.2d])

        then: 'no constraint-bearing term survives -> no call'
        res.taxon == null
        res.label == 'Unknown'
        res.constrainedPresent == 0
    }
}
