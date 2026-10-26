package gspa.ontology

import org.semanticweb.owlapi.model.IRI
import org.semanticweb.owlapi.reasoner.OWLReasoner
import spock.lang.Specification
import spock.lang.TempDir

import java.lang.reflect.InvocationTargetException
import java.lang.reflect.Proxy
import java.util.concurrent.CountDownLatch
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

/**
 * Concurrent has_part extraction must do the ELK work once, return the full
 * entailed relation, and never serve stale or damaged cache content.
 *
 * Fixture biology (only the first two has_part facts are asserted):
 *   DNA replication (GO:0006260) has_part some DNA replication initiation (GO:0006270)
 *   DNA replication initiation (GO:0006270) is_a DNA metabolic process (GO:0006259) is_a cellular process (GO:0009987)
 *   DNA-templated DNA replication (GO:0006261) is_a DNA replication
 *   translation (GO:0006412) has_part some tRNA aminoacylation for protein translation (GO:0006418)
 * Entailed: has_part is inherited by subclasses (GO:0006261) and widens along
 * the filler hierarchy (has_part some initiation => has_part some DNA metabolic process).
 */
class GoReasonerConcurrencySpec extends Specification {
    @TempDir File cache
    @TempDir File work

    static final String FIXTURE = '''Prefix(:=<http://purl.obolibrary.org/obo/>)
Ontology(<http://purl.obolibrary.org/obo/go-hp-test.owl>
Declaration(ObjectProperty(:BFO_0000051))
Declaration(Class(:GO_0009987))
Declaration(Class(:GO_0006259))
Declaration(Class(:GO_0006270))
Declaration(Class(:GO_0006260))
Declaration(Class(:GO_0006261))
Declaration(Class(:GO_0006412))
Declaration(Class(:GO_0006418))
SubClassOf(:GO_0006259 :GO_0009987)
SubClassOf(:GO_0006270 :GO_0006259)
SubClassOf(:GO_0006261 :GO_0006260)
SubClassOf(:GO_0006260 ObjectSomeValuesFrom(:BFO_0000051 :GO_0006270))
SubClassOf(:GO_0006412 ObjectSomeValuesFrom(:BFO_0000051 :GO_0006418))
)
'''

    static final Set<List<String>> ENTAILED = [
        ['GO:0006260', 'GO:0006270'], ['GO:0006260', 'GO:0006259'], ['GO:0006260', 'GO:0009987'],
        ['GO:0006261', 'GO:0006270'], ['GO:0006261', 'GO:0006259'], ['GO:0006261', 'GO:0009987'],
        ['GO:0006412', 'GO:0006418'],
    ] as Set

    /** A different ontology (own IRI): translation initiation is_a translation, so it inherits has_part aminoacylation. */
    static final String FIXTURE_EXTENDED = FIXTURE
        .replace('go-hp-test.owl', 'go-hp-extended.owl')
        .replace('Declaration(Class(:GO_0006418))\n', 'Declaration(Class(:GO_0006418))\nDeclaration(Class(:GO_0006413))\n')
        .replace('SubClassOf(:GO_0006259 :GO_0009987)\n', 'SubClassOf(:GO_0006259 :GO_0009987)\nSubClassOf(:GO_0006413 :GO_0006412)\n')

    static final Set<List<String>> ENTAILED_EXTENDED = ENTAILED + [['GO:0006413', 'GO:0006418']]

    private GoOntology load() {
        def file = new File(work, 'go-hp-test.ofn')
        file.text = FIXTURE
        def go = new GoOntology()
        go.loadOwl(file)
        go
    }

    private GoOntology loadExtended() {
        def file = new File(work, 'go-hp-extended.ofn')
        file.text = FIXTURE_EXTENDED
        def go = new GoOntology()
        go.loadOwl(file)
        go
    }

    /** Wraps the reasoner so tests can count how many class-expression queries reach ELK. */
    private static AtomicInteger countQueries(GoReasoner gr) {
        def calls = new AtomicInteger()
        OWLReasoner real = gr.reasoner
        gr.reasoner = (OWLReasoner) Proxy.newProxyInstance(OWLReasoner.classLoader, [OWLReasoner] as Class[]) { proxy, method, args ->
            if (method.name == 'getSubClasses') calls.incrementAndGet()
            try {
                method.invoke(real, args)
            } catch (InvocationTargetException e) {
                throw e.targetException
            }
        }
        calls
    }

    /** Counts calls to the named method on the reasoner currently installed in gr. */
    private static AtomicInteger spy(GoReasoner gr, String methodName) {
        def calls = new AtomicInteger()
        OWLReasoner real = gr.reasoner
        gr.reasoner = (OWLReasoner) Proxy.newProxyInstance(OWLReasoner.classLoader, [OWLReasoner] as Class[]) { proxy, method, args ->
            if (method.name == methodName) calls.incrementAndGet()
            try {
                method.invoke(real, args)
            } catch (InvocationTargetException e) {
                throw e.targetException
            }
        }
        calls
    }

    private static Set<List<String>> asSet(List<Map.Entry<String, String>> pairs) {
        pairs.collect { [it.key, it.value] } as Set
    }

    def "concurrent clients trigger one extraction and all see the full entailed has_part relation"() {
        given:
        def go = load()
        def gr = new GoReasoner(go)
        gr.initialize()
        def calls = countQueries(gr)
        int clients = 8
        def start = new CountDownLatch(1)
        def pool = Executors.newFixedThreadPool(clients)

        when:
        def futures = (1..clients).collect { pool.submit({ start.await(); gr.getHasPartPairs() } as java.util.concurrent.Callable) }
        start.countDown()
        def results = futures.collect { it.get(60, TimeUnit.SECONDS) }

        then: 'ELK was queried once per GO term, not once per client per term'
        calls.get() == go.getAllGoTerms().size()

        and: 'every client holds the same complete entailed relation, including inherited and widened pairs'
        results.every { it.is(results[0]) }
        asSet(results[0]) == ENTAILED
        results[0].size() == ENTAILED.size()

        cleanup:
        pool.shutdownNow()
        gr.dispose()
    }

    def "the shared result cannot be mutated by a client"() {
        given:
        def gr = new GoReasoner(load())
        gr.initialize()
        def first = gr.getHasPartPairs()

        when:
        first << new AbstractMap.SimpleEntry<>('GO:9999999', 'GO:0006418')

        then:
        thrown(UnsupportedOperationException)

        when:
        first[0].setValue('GO:9999999')

        then:
        thrown(UnsupportedOperationException)

        and:
        asSet(gr.getHasPartPairs()) == ENTAILED

        cleanup:
        gr.dispose()
    }

    def "the same ontology version is served from memory without more ELK queries"() {
        given:
        def gr = new GoReasoner(load())
        gr.initialize()
        def calls = countQueries(gr)
        gr.getHasPartPairs()
        int afterFirst = calls.get()

        when:
        gr.getHasPartPairs()

        then:
        calls.get() == afterFirst

        cleanup:
        gr.dispose()
    }

    def "replacing the ontology object re-initializes the reasoner and never attributes old entailments to the new ontology"() {
        given:
        def gr = new GoReasoner(load())
        gr.cacheDir = cache
        gr.initialize()
        def disposed = spy(gr, 'dispose')
        def before = asSet(gr.getHasPartPairs())
        def changed = loadExtended()

        when: 'the same GoReasoner now points at an ontology with different axioms'
        gr.goOntology = changed
        def after = asSet(gr.getHasPartPairs())

        then:
        before == ENTAILED
        after == ENTAILED_EXTENDED
        disposed.get() == 1
        new File(cache, 'has_part_pairs.tsv').readLines()[0] == "# ontology_sha256=${changed.sourceSha256}"

        cleanup:
        gr.dispose()
    }

    def "reloading different axioms into the same GoOntology is detected"() {
        given:
        def go = load()
        def gr = new GoReasoner(go)
        gr.initialize()
        assert asSet(gr.getHasPartPairs()) == ENTAILED

        when:
        def file = new File(work, 'go-hp-extended.ofn')
        file.text = FIXTURE_EXTENDED
        go.loadOwl(file)

        then:
        asSet(gr.getHasPartPairs()) == ENTAILED_EXTENDED

        cleanup:
        gr.dispose()
    }

    def "in-place axiom edits are detected even when the digest is unchanged"() {
        given:
        def go = load()
        def gr = new GoReasoner(go)
        gr.initialize()
        assert asSet(gr.getHasPartPairs()) == ENTAILED

        when: 'translation initiation is added as a subclass of translation'
        go.manager.addAxiom(go.ontology, go.factory.getOWLSubClassOfAxiom(
            go.factory.getOWLClass(IRI.create('http://purl.obolibrary.org/obo/GO_0006413')), go.getGoClass('GO:0006412')))

        then:
        asSet(gr.getHasPartPairs()) == ENTAILED + [['GO:0006413', 'GO:0006418']]

        cleanup:
        gr.dispose()
    }

    def "initializing twice disposes the previous reasoner, and a failed initialization keeps the working one"() {
        given:
        def go = load()
        def gr = new GoReasoner(go)
        gr.initialize()
        def firstDisposed = spy(gr, 'dispose')

        when:
        gr.initialize()

        then:
        firstDisposed.get() == 1

        when: 'creation fails for the next initialization'
        def secondDisposed = spy(gr, 'dispose')
        def saved = go.ontology
        go.ontology = null
        gr.initialize()

        then:
        thrown(Throwable)
        secondDisposed.get() == 0

        when:
        go.ontology = saved

        then: 'the retained reasoner still serves the unchanged ontology'
        asSet(gr.getHasPartPairs()) == ENTAILED

        cleanup:
        gr.dispose()
    }

    def "clients stay consistent while another thread re-initializes"() {
        given:
        def gr = new GoReasoner(load())
        gr.initialize()
        def pool = Executors.newFixedThreadPool(9)
        def start = new CountDownLatch(1)

        when:
        def readers = (1..8).collect {
            pool.submit({
                start.await()
                (1..20).collect { asSet(gr.getHasPartPairs()) }
            } as java.util.concurrent.Callable)
        }
        def reinit = pool.submit({ start.await(); 10.times { gr.initialize() } } as Runnable)
        start.countDown()
        reinit.get(120, TimeUnit.SECONDS)
        def seen = readers.collectMany { it.get(120, TimeUnit.SECONDS) } as Set

        then:
        seen == [ENTAILED] as Set

        cleanup:
        pool.shutdownNow()
        gr.dispose()
    }

    def "malformed or truncated disk caches are rejected and replaced by a correct recomputation"() {
        given:
        def go = load()
        def seed = new GoReasoner(go)
        seed.cacheDir = cache
        seed.initialize()
        seed.getHasPartPairs()
        def file = new File(cache, 'has_part_pairs.tsv')
        List<String> good = file.readLines()
        seed.dispose()

        when:
        file.text = corrupt(good).join('\n') + '\n'
        def gr = new GoReasoner(go)
        gr.cacheDir = cache
        gr.initialize()

        then:
        !gr.loadHasPartCache()
        asSet(gr.getHasPartPairs()) == ENTAILED
        file.readLines() as Set == good as Set  // rewritten as a valid cache

        cleanup:
        gr.dispose()

        where:
        label                        | corrupt
        'missing end marker'         | { List<String> l -> l[0..-2] }
        'truncated mid-record'       | { List<String> l -> l[0..-3] + [l[-2].take(5)] }
        'dropped record, same count' | { List<String> l -> [l[0], l[1]] + l[3..-1] }
        'count not numeric'          | { List<String> l -> [l[0], '# pairs=many'] + l[2..-1] }
        'single-column record'       | { List<String> l -> l[0..1] + ['GO:0006260'] + l[3..-1] }
        'non-GO record'              | { List<String> l -> l[0..1] + ["GO:0006260\tnot-a-term"] + l[3..-1] }
        'duplicate record'           | { List<String> l -> l[0..2] + [l[2]] + l[4..-1] }
        'pre-footer legacy format'   | { List<String> l -> l[0..-2].findAll { !it.startsWith('# pairs=') } }
        'malformed GO id'            | { List<String> l -> l[0..1] + ["GO:banana\tGO:0006270"] + l[3..-1] }
        'short GO id'                | { List<String> l -> l[0..1] + ["GO:6260\tGO:0006270"] + l[3..-1] }
        'valid-format unknown term'  | { List<String> l -> l[0..1] + ["GO:0006260\tGO:7777777"] + l[3..-1] }
        'padded GO id'               | { List<String> l -> l[0..1] + ["GO:0006260\tGO:0006270 "] + l[3..-1] }
        'empty file'               | { List<String> l -> [] }
    }

    def "a failed extraction is not cached and a later initialization succeeds"() {
        given:
        def gr = new GoReasoner(load())

        when:
        gr.getHasPartPairs()

        then:
        thrown(IllegalStateException)

        when:
        gr.initialize()

        then:
        asSet(gr.getHasPartPairs()) == ENTAILED

        cleanup:
        gr.dispose()
    }

    def "initialize and reset discard derived results; dispose keeps valid ones"() {
        given:
        def gr = new GoReasoner(load())
        gr.initialize()
        def calls = countQueries(gr)
        gr.getHasPartPairs()
        int once = calls.get()

        when: 'dispose does not invalidate results for the same ontology'
        gr.dispose()

        then:
        asSet(gr.getHasPartPairs()) == ENTAILED

        when: 'reset forces a fresh derivation, which needs a live reasoner'
        gr.resetHasPartCache()
        gr.getHasPartPairs()

        then:
        thrown(IllegalStateException)

        when:
        gr.initialize()
        def calls2 = countQueries(gr)
        gr.getHasPartPairs()

        then:
        once > 0
        calls2.get() == once

        cleanup:
        gr.dispose()
    }
}
