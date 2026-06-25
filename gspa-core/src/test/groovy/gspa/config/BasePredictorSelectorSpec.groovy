package gspa.config

import spock.lang.Specification

/**
 * NeuralConfig.basePredictor selects DeepGO-PlusPlus full (GPU) or light (CPU)
 * as the workflow's base function predictor, mutually exclusively.
 */
class BasePredictorSelectorSpec extends Specification {

    def 'full enables DeepGO-PlusPlus and disables Light'() {
        given:
        def neural = new GspaConfig.NeuralConfig(basePredictor: sel)

        when:
        neural.resolveBasePredictor()

        then:
        neural.deepGoPlusPlus.enabled
        !neural.deepGoPlusPlusLight.enabled

        where:
        sel << ['full', 'FULL', 'deepgo-plusplus', 'gpu']
    }

    def 'light enables Light and disables the full predictor'() {
        given:
        def neural = new GspaConfig.NeuralConfig(basePredictor: sel)

        when:
        neural.resolveBasePredictor()

        then:
        neural.deepGoPlusPlusLight.enabled
        !neural.deepGoPlusPlus.enabled

        where:
        sel << ['light', 'cpu', 'deepgo-plusplus-light']
    }

    def 'none leaves explicit enabled flags untouched'() {
        given:
        def neural = new GspaConfig.NeuralConfig(basePredictor: 'none')
        neural.deepGoPlusPlus.enabled = true

        when:
        neural.resolveBasePredictor()

        then:
        neural.deepGoPlusPlus.enabled
        !neural.deepGoPlusPlusLight.enabled
    }

    def 'default selector is none (no-op)'() {
        given:
        def neural = new GspaConfig.NeuralConfig()

        expect:
        neural.basePredictor == 'none'

        when:
        neural.resolveBasePredictor()

        then:
        !neural.deepGoPlusPlus.enabled
        !neural.deepGoPlusPlusLight.enabled
    }

    def 'switching base from light to full flips both flags'() {
        given: 'light was selected first'
        def neural = new GspaConfig.NeuralConfig(basePredictor: 'light')
        neural.resolveBasePredictor()

        when: 'the base is changed to full and re-resolved'
        neural.basePredictor = 'full'
        neural.resolveBasePredictor()

        then:
        neural.deepGoPlusPlus.enabled
        !neural.deepGoPlusPlusLight.enabled
    }

    def 'an unknown selector is rejected'() {
        when:
        new GspaConfig.NeuralConfig(basePredictor: 'medium').resolveBasePredictor()

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('medium')
    }
}
