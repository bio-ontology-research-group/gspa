package gspa.io

import gspa.model.AnnotationType
import spock.lang.Specification

class GafReaderSpec extends Specification {

    def "should read GAF annotations"() {
        given:
        def gafFile = new File(getClass().getResource('/test-genomes/test_annotations.gaf').toURI())

        when:
        def annotations = GafReader.readGaf(gafFile)

        then:
        annotations.keySet().size() == 4
        annotations['protein_001'].size() == 2
        annotations['protein_003'].size() == 3

        // Check GO aspects
        def p3anns = annotations['protein_003']
        p3anns.find { it.value == 'GO:0006412' }.goAspect == 'BP'
        p3anns.find { it.value == 'GO:0003735' }.goAspect == 'MF'
        p3anns.find { it.value == 'GO:0005840' }.goAspect == 'CC'

        // All should be GO type
        annotations.values().flatten().every { it.type == AnnotationType.GO }
    }

    def "should assign evidence code scores"() {
        expect:
        GafReader.evidenceCodeScore('IDA') == 1.0
        GafReader.evidenceCodeScore('IEA') == 0.4
        GafReader.evidenceCodeScore('ISS') == 0.6
        GafReader.evidenceCodeScore('TAS') == 0.7
    }
}
