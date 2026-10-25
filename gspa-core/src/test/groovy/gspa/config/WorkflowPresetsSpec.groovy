package gspa.config

import spock.lang.Specification
import spock.lang.TempDir

class WorkflowPresetsSpec extends Specification {
    @TempDir File tempDir

    def 'sparse user YAML preserves preset values and explicit false wins'() {
        given:
        File yaml = new File(tempDir, 'user.yaml')
        yaml.text = '''
quality:
  consistency:
    enforce: false
predictors:
  domains:
    applications: [Pfam]
output:
  outputDir: custom-output
'''

        when:
        def config = ConfigLoader.buildConfig(yaml, null, [:], 'eukaryote')

        then:
        config.input.kingdom == 'eukaryote'
        config.input.requireExistingGenes
        !config.predictors.operons.enabled
        !config.predictors.pathway.enabled
        !config.quality.consistency.enforce
        config.quality.consistency.organismTaxon == 'eukaryote'
        config.quality.completeness.profile == 'eukaryote'
        config.predictors.domains.applications == ['Pfam']
        config.output.outputDir == 'custom-output'
    }

    def 'legacy kingdom preset is actually loaded and survives a partial user config'() {
        given:
        File yaml = new File(tempDir, 'user.yaml')
        yaml.text = 'output: {outputDir: custom-output}'

        expect:
        !ConfigLoader.loadPreset('eukaryote').predictors.operons.enabled
        !ConfigLoader.buildConfig(yaml, 'eukaryote').predictors.operons.enabled
    }

    def 'unknown presets and missing configuration files fail instead of using defaults'() {
        when:
        ConfigLoader.buildConfig(null, null, [:], 'eukarytoe')
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('gspa presets')

        when:
        ConfigLoader.buildConfig(new File(tempDir, 'missing.yaml'))
        then:
        thrown(IllegalArgumentException)
    }

    def 'mixed communities and single bins use distinct biological assumptions'() {
        when:
        def mixed = ConfigLoader.buildConfig(null, null, [:], 'metagenomics')
        def mag = ConfigLoader.buildConfig(null, 'archaea', [:], 'mag')

        then:
        mixed.input.type == 'metagenome'
        !mixed.quality.enabled
        !mixed.quality.consistency.enforce
        mixed.quality.completeness.profile == 'skip'
        !mixed.predictors.operons.enabled
        mag.input.type == 'mag'
        mag.input.kingdom == 'archaea'
        mag.quality.completeness.profile == 'archaea'
        mag.predictors.operons.enabled
        !mag.quality.completeness.enforce
        !mag.quality.coherence.enforce
        !mag.predictors.pathway.crossfeeding
    }

    def 'CLI overrides and reusable YAML preserve settings'() {
        given:
        def config = ConfigLoader.buildConfig(null, null,
            [database: 'reference.dmnd', 'output-dir': 'run-out'], 'bacteria')
        File saved = new File(tempDir, 'resolved.yaml')
        ConfigLoader.writeConfig(config, saved)

        expect:
        config.predictors.similarity.enabled
        config.predictors.similarity.database == 'reference.dmnd'
        ConfigLoader.buildConfig(saved).output.outputDir == 'run-out'
        ConfigLoader.buildConfig(saved).quality.consistency.enforce
    }

    def 'explicit YAML null clears an inherited optional value'() {
        given:
        File yaml = new File(tempDir, 'clear.yaml')
        yaml.text = 'quality: {consistency: {organismTaxon: null}}'
        expect:
        ConfigLoader.buildConfig(yaml, null, [:], 'bacteria').quality.consistency.organismTaxon == null
    }
}
