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

    def "domains.database defaults to null and is not required"() {
        expect:
        new GspaConfig().predictors.domains.database == null
        ConfigLoader.loadFromString('predictors:\n  domains:\n    tool: hmmer\n').predictors.domains.database == null
    }

    def "domains.database round-trips through YAML, buildConfig merge and config merge"() {
        given:
        def yaml = '''
predictors:
  domains:
    tool: hmmer
    database: /db/Pfam-A.hmm
'''
        def userFile = File.createTempFile('gspa-domains-', '.yaml')
        userFile.text = yaml

        when: "loaded directly"
        def loaded = ConfigLoader.loadFromString(yaml)

        and: "serialised and read back"
        def again = ConfigLoader.loadFromString(ConfigLoader.toYaml(loaded))

        and: "merged through the layered builder (defaults -> kingdom preset -> user file)"
        def built = ConfigLoader.buildConfig(userFile, 'bacteria')

        and: "merged over a base config that does not set it"
        def merged = ConfigLoader.merge(new GspaConfig(), loaded)

        then:
        loaded.predictors.domains.database == '/db/Pfam-A.hmm'
        again.predictors.domains.database == '/db/Pfam-A.hmm'
        again.predictors.domains.tool == 'hmmer'
        built.predictors.domains.database == '/db/Pfam-A.hmm'
        built.predictors.domains.tool == 'hmmer'
        merged.predictors.domains.database == '/db/Pfam-A.hmm'

        cleanup:
        userFile.delete()
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
