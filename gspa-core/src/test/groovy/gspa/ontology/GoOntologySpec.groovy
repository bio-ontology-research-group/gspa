package gspa.ontology

import spock.lang.Specification

class GoOntologySpec extends Specification {

    def "getLabel returns the rdfs:label literal, not Boolean.toString()"() {
        // Regression: in Groovy, OWLLiteral.literal resolves to
        // OWLAnnotationValue.isLiteral() (boolean true) instead of
        // OWLLiteral.getLiteral() (the lexical form). The fix is to call
        // getLiteral() explicitly. Without the fix, every label was the
        // literal string "true".
        given:
        def go = new GoOntology()
        def oboFile = new File(getClass().getResource('/test-ontology/go-tiny.obo').toURI())
        go.loadOwl(oboFile)

        when:
        String label = go.getLabel('GO:0006259')

        then:
        label == 'DNA metabolic process'
        label != 'true'
    }

    def "buildLabelCache populates real labels (not Boolean.toString())"() {
        given:
        def go = new GoOntology()
        def oboFile = new File(getClass().getResource('/test-ontology/go-tiny.obo').toURI())
        go.loadOwl(oboFile)

        when:
        // Force buildLabelCache by retrieving any label (it's invoked from
        // loadOwl in some configs; if not, the dynamic getLabel path is
        // exercised — both should be string-correct).
        String l1 = go.getLabel('GO:0009987')
        String l2 = go.getLabel('GO:0008150')

        then:
        l1 == 'cellular process'
        l2 == 'biological_process'
        ![l1, l2].any { it == 'true' || it == 'false' }
    }
}
