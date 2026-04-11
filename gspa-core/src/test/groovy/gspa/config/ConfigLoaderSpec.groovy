package gspa.config

import spock.lang.Specification

class ConfigLoaderSpec extends Specification {

    def "should create default config"() {
        when:
        def config = new GspaConfig()

        then:
        config.input.format == 'fasta'
        config.input.kingdom == 'auto'
        config.predictors.similarity.enabled
        config.predictors.similarity.tool == 'diamond'
        config.quality.coherence.process
        config.quality.consistency.taxonConstraints
    }

    def "should load config from YAML string"() {
        given:
        def yaml = '''
input:
  kingdom: bacteria
  type: mag
predictors:
  similarity:
    tool: mmseqs2
    evalue: 1.0e-10
  operons:
    maxIntergenicDistance: 500
quality:
  completeness:
    profile: bacteria
    addTerms:
      - GO:0009058
'''

        when:
        def config = ConfigLoader.loadFromString(yaml)

        then:
        config.input.kingdom == 'bacteria'
        config.input.type == 'mag'
        config.predictors.similarity.tool == 'mmseqs2'
        config.predictors.similarity.evalue == 1.0e-10
        config.predictors.operons.maxIntergenicDistance == 500
        config.quality.completeness.addTerms == ['GO:0009058']
    }

    def "should merge configs with override taking precedence"() {
        given:
        def base = new GspaConfig()
        base.predictors.similarity.tool = 'diamond'
        base.predictors.similarity.evalue = 1e-5

        def override = new GspaConfig()
        override.predictors.similarity.tool = 'mmseqs2'

        when:
        def merged = ConfigLoader.merge(base, override)

        then:
        merged.predictors.similarity.tool == 'mmseqs2'
    }

    def "should resolve organism domain"() {
        given:
        def config = new GspaConfig()

        when:
        config.input.kingdom = 'bacteria'

        then:
        config.resolveOrganismDomain() == gspa.model.OrganismDomain.BACTERIA
    }
}
