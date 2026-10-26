package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.config.ConfigLoader
import gspa.config.GspaConfig
import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Unroll

class IntegrationWorkflowConfigSpec extends Specification {
    @TempDir File dir

    def 'the published key list is exactly the typed configuration schema'() {
        given:
        def published = new ObjectMapper().readValue(
            getClass().getResourceAsStream('/workflow/integration-workflow-config.keys.json'), List) as Set

        expect: 'web and Nextflow validate against this list, so it cannot drift from the class'
        new IntegrationWorkflowConfig().toMap().keySet() == published
    }

    def 'configuration survives a JSON round trip and rejects names it does not know'() {
        given:
        def config = new IntegrationWorkflowConfig(enabledPriors: 'essentiality', fbaMedium: 'LB',
            moduleFiles: [new File('m.tsv')], theta: [dark_matter: [bf_min: 3.0d]])
        File json = new File(dir, 'c.json')
        config.writeJson(json)

        expect:
        IntegrationWorkflowConfig.fromJson(json).toMap() == config.toMap()

        when:
        IntegrationWorkflowConfig.fromMap([enablePriors: 'x'])

        then:
        thrown(Exception)
    }

    @Unroll
    def 'invalid setting #name is rejected'() {
        when:
        new IntegrationWorkflowConfig(values).validate()

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains(fragment)

        where:
        name                     | values                                          | fragment
        'unknown prior'          | [enabledPriors: 'essentiality,homology']        | "Unknown prior 'homology'"
        'all is not a prior'     | [enabledPriors: 'all']                          | "Unknown prior 'all'"
        'enforcement mode'       | [enforceConsistencyMode: 'delete']              | 'enforce-consistency-mode'
        'dark matter mode'       | [darkMatterMode: 'everywhere']                  | 'dark-matter mode'
        'outer loop w/o suggester' | [iterateGapseq: true]                         | '--iterate-gapseq requires --dark-matter'
        'reps without clustering' | [gapseqTarget: 'reps']                         | 'requires --intragenome-cluster'
        'beam width'             | [beamWidth: 0]                                  | '--beam-width'
    }

    def 'a named reference file that does not exist is an error, not an absent prior'() {
        when:
        new IntegrationWorkflowConfig(operonsFile: new File(dir, 'none.tsv')).validate()

        then:
        thrown(IllegalArgumentException)
    }

    def 'inline theta overlays the theta file section by section'() {
        given:
        File file = new File(dir, 'theta.json')
        file.text = '{"max_iter": 4, "dark_matter": {"bf_min": 5.0, "gamma_in_p": 20.0}}'
        def config = new IntegrationWorkflowConfig(thetaFile: file, theta: [dark_matter: [bf_min: 7.0d]], damping: 0.25d)

        expect:
        config.resolveTheta() == [max_iter: 4, dark_matter: [bf_min: 7.0d, gamma_in_p: 20.0d], damping: 0.25d]
    }

    def 'an annotation run implies the same workflow settings as gspa integrate'() {
        given:
        def config = new GspaConfig()
        config.integration.maxIter = 3
        config.integration.darkMatter.enabled = true
        config.integration.darkMatter.mode = 'top-q'
        config.integration.darkMatter.bfMin = 4.0d
        config.integration.outerLoop.enabled = true
        config.integration.outerLoop.maxIter = 2
        config.quality.consistency.organismTaxon = 'bacteria'
        config.quality.completeness.profile = 'bacteria'
        config.integration.workflow.enabledPriors = 'consistency'
        config.integration.workflow.fbaMedium = 'LB'

        when:
        def w = IntegrationWorkflowConfig.fromGspa(config)

        then:
        w.maxIter == 3
        w.darkMatter && w.darkMatterMode == 'top-q'
        w.resolveTheta().dark_matter.bf_min == 4.0d
        w.iterateGapseq && w.maxGapseqIter == 2
        w.organismTaxon == 'bacteria'
        w.essentialProfile == 'bacteria'
        w.enabledPriors == 'consistency'
        w.fbaMedium == 'LB'
    }

    def 'the integration block of a YAML run configuration is read strictly'() {
        when:
        def config = ConfigLoader.loadFromString('''
integration:
  enabled: true
  workflow:
    enabledPriors: essentiality
    enforceStage2: true
''')

        then:
        config.integration.enabled
        config.integration.workflow.enforceStage2
        config.integration.workflow.enabledPriors == 'essentiality'

        when:
        ConfigLoader.loadFromString('''
integration:
  workflow:
    enablePriors: essentiality
''')

        then:
        thrown(Exception)
    }
}
