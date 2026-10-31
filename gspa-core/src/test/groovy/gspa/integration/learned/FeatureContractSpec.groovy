package gspa.integration.learned

import spock.lang.Specification

/** The versioned feature contract: uniqueness, documentation, fingerprint, and sync with the trainer's copy. */
class FeatureContractSpec extends Specification {

    def "version 3 excludes descriptive cohort counts and their unweighted ratios"() {
        expect:
        FeatureContract.CURRENT.id == 'gspa.context-features/3'
        FeatureContract.CURRENT.names().findAll { it.startsWith('cohort.') }.every { it.startsWith('cohort.w.') }
        !FeatureContract.CURRENT.names().contains('cohort.w.n_peer_genomes')
        !FeatureContract.CURRENT.names().contains('cohort.w.bg_windows')
    }

    def "a genuine version 2 fixture cannot be silently reinterpreted under version 3"() {
        when:
        LearnedTestSupport.mechanismFixture('legacy_model_v2.json')

        then:
        def e = thrown(ModelUnavailableException)
        e.message.contains('feature contract gspa.context-features/2')
        e.message.contains('does not match gspa.context-features/3')
    }

    def "every feature is documented, uniquely named and typed by group"() {
        expect:
        FeatureContract.CURRENT.names().toSet().size() == FeatureContract.CURRENT.size()
        FeatureContract.CURRENT.features.every { it.description?.length() > 10 && it.group != null }
        FeatureContract.CURRENT.features*.group.toSet() == FeatureGroup.values().toList().toSet()
    }

    def "baseline sequence features carry no posterior / calibrated-probability / bayes-factor column"() {
        expect:
        FeatureContract.CURRENT.names().every { n ->
            !(n.toLowerCase() =~ /posterior|likelihood|bayes|calibrated|logodds|floor|consensus|rlgc/)
        }
    }

    def "fingerprint tracks names and groups but not descriptions"() {
        given:
        def c = FeatureContract.CURRENT
        def reworded = new FeatureContract(c.id, c.features.collect { new FeatureSpec(it.name, it.group, 'reworded text') })
        def renamed = new FeatureContract(c.id, c.features.collect { new FeatureSpec(it.name + 'x', it.group, it.description) })

        expect:
        reworded.fingerprint() == c.fingerprint()
        renamed.fingerprint() != c.fingerprint()
    }

    def "contract JSON round-trips"() {
        expect:
        FeatureContract.fromJson(FeatureContract.CURRENT.toJson()).fingerprint() == FeatureContract.CURRENT.fingerprint()
        FeatureContract.fromJson(FeatureContract.LOCUS_GROUP.toJson()).names() == FeatureContract.LOCUS_GROUP.names()
    }

    def "writes the contract for the Python trainer, and the committed copy matches the code"() {
        given: "build output used by `train_context_model.py make-mechanism-fixture --contract ...`"
        def out = new File('build/learned-fixtures')
        out.mkdirs()
        new File(out, 'feature_contract_v1.json').text = FeatureContract.CURRENT.toJson()
        def committed = getClass().getResourceAsStream('/learned/feature_contract_v1.json')

        expect: "the fixture copy in test resources must be regenerated whenever the contract changes"
        committed != null
        committed.getText('UTF-8') == FeatureContract.CURRENT.toJson()
    }
}
