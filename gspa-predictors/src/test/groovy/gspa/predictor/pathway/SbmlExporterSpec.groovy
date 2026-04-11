package gspa.predictor.pathway

import gspa.model.*
import spock.lang.Specification

class SbmlExporterSpec extends Specification {

    def "should export valid SBML from genome annotations"() {
        given:
        def genome = new Genome(id: 'test_genome')
        def contig = new Contig(id: 'c1')
        def p1 = new Protein(id: 'p1', sequence: 'MFAKE')
        p1.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:2.7.1.1', source: 'test', score: 0.9))
        p1.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:1.2.1.12', source: 'test', score: 0.8))
        def p2 = new Protein(id: 'p2', sequence: 'MOTHER')
        p2.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:4.2.1.11', source: 'test', score: 0.7))
        contig.addProtein(p1)
        contig.addProtein(p2)
        genome.addContig(contig)

        def outputFile = File.createTempFile("model_", ".xml")

        when:
        SbmlExporter.exportSbml(genome, outputFile)

        then:
        outputFile.exists()
        def xml = outputFile.text
        xml.contains('<?xml version="1.0"')
        xml.contains('<sbml')
        xml.contains('level="3"')
        xml.contains('test_genome')
        xml.contains('EC: EC:2.7.1.1')
        xml.contains('Gene: p1')
        xml.contains('EC: EC:4.2.1.11')
        xml.contains('listOfReactions')
        xml.contains('cytoplasm')

        cleanup:
        outputFile.delete()
    }
}
