package gspa.predictor

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Protein
import spock.lang.Specification

class PredictorRegistrySpec extends Specification {

    def "should register and retrieve predictors"() {
        given:
        def registry = new PredictorRegistry()
        def predictor = new MockPredictor(name: 'test-predictor')

        when:
        registry.register(predictor)

        then:
        registry.get('test-predictor') == predictor
        registry.getAll().size() == 1
    }

    def "should filter available predictors"() {
        given:
        def registry = new PredictorRegistry()
        registry.register(new MockPredictor(name: 'available', available: true))
        registry.register(new MockPredictor(name: 'unavailable', available: false))

        when:
        def available = registry.getAvailable()

        then:
        available.size() == 1
        available[0].name == 'available'
    }

    def "should filter by enable/disable lists"() {
        given:
        def registry = new PredictorRegistry()
        registry.register(new MockPredictor(name: 'a', available: true))
        registry.register(new MockPredictor(name: 'b', available: true))
        registry.register(new MockPredictor(name: 'c', available: true))

        expect:
        registry.getEnabled(['a', 'b'] as Set, [] as Set)*.name == ['a', 'b']
        registry.getEnabled([] as Set, ['b'] as Set)*.name == ['a', 'c']
    }

    def "should report status"() {
        given:
        def registry = new PredictorRegistry()
        registry.register(new MockPredictor(name: 'yes', available: true))
        registry.register(new MockPredictor(name: 'no', available: false))

        when:
        def status = registry.status()

        then:
        status == [yes: true, no: false]
    }
}

class MockPredictor implements Predictor {
    String name
    boolean available = true

    String getVersion() { '1.0' }
    Set<AnnotationType> getOutputTypes() { [AnnotationType.GO] as Set }
    boolean isAvailable() { available }
    List<Annotation> predict(Protein protein) { [] }
}
