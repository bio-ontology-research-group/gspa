package gspa.config

import gspa.model.OrganismDomain
import spock.lang.Specification

class EssentialFunctionsSpec extends Specification {

    def "should load default bacteria essential functions"() {
        when:
        def ef = EssentialFunctions.getDefault(OrganismDomain.BACTERIA)

        then:
        ef.profileName == 'bacteria'
        ef.goTerms.size() > 10
        ef.goTerms.contains('GO:0006412') // translation
        ef.goTerms.contains('GO:0006260') // DNA replication
        ef.goTerms.contains('GO:0051301') // cell division
        ef.categories.contains('Core')
    }

    def "should load default eukaryote essential functions"() {
        when:
        def ef = EssentialFunctions.getDefault(OrganismDomain.EUKARYA)

        then:
        ef.profileName == 'eukaryote'
        ef.goTerms.contains('GO:0006412') // translation (universal)
        ef.goTerms.contains('GO:0007049') // cell cycle (eukaryote-specific)
        ef.goTerms.contains('GO:0008380') // RNA splicing (eukaryote-specific)
        !ef.goTerms.contains('GO:0009306') // protein secretion (removed for euk)
    }

    def "should apply runtime modifications"() {
        given:
        def ef = EssentialFunctions.getDefault(OrganismDomain.BACTERIA)
        int originalSize = ef.goTerms.size()

        when:
        def modified = ef.withModifications(['GO:9999999'], ['GO:0006412'])

        then:
        modified.goTerms.contains('GO:9999999')
        !modified.goTerms.contains('GO:0006412')
        modified.goTerms.size() == originalSize // +1 -1
        // Original unchanged
        ef.goTerms.contains('GO:0006412')
    }

    def "should get terms by category"() {
        given:
        def ef = EssentialFunctions.getDefault(OrganismDomain.BACTERIA)

        when:
        def coreTerms = ef.getTermsByCategory('Core')

        then:
        coreTerms.size() > 10
        coreTerms.contains('GO:0006412')
    }
}
