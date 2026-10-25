package gspa.ontology

import spock.lang.Specification
import spock.lang.TempDir

class GoReasonerCacheSpec extends Specification {
    @TempDir File cache

    def "a persisted inference cache is reusable only for the same ontology"() {
        given:
        def go = new GoOntology()
        go.loadOwl(new File(getClass().getResource('/test-ontology/go-tiny.obo').toURI()))
        def first = new GoReasoner(go)
        first.cacheDir = cache
        first.initialize()

        when:
        def expected = first.getHasPartPairs()
        def second = new GoReasoner(go)
        second.cacheDir = cache

        then:
        second.loadHasPartCache()
        second.getHasPartPairs() == expected
        new File(cache, 'has_part_pairs.tsv').readLines()[0] == "# ontology_sha256=${go.sourceSha256}"

        when:
        go.sourceSha256 = 'different-ontology'
        def changed = new GoReasoner(go)
        changed.cacheDir = cache

        then:
        !changed.loadHasPartCache()

        cleanup:
        first.dispose()
    }
}
