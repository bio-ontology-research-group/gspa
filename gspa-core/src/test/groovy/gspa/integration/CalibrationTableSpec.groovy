package gspa.integration

import spock.lang.Specification

class CalibrationTableSpec extends Specification {

    def table = new CalibrationTable()

    def "calibrated probabilities stay in (0, 1)"() {
        when:
        def p0 = table.calibrate('diamond', 0.0)
        def p1 = table.calibrate('diamond', 1.0)

        then:
        p0 > 0.0 && p0 < 1.0
        p1 > 0.0 && p1 < 1.0
    }

    def "calibration is monotone non-decreasing in raw score"() {
        expect:
        def scores = (0..10).collect { it / 10.0 }
        def calibrated = scores.collect { table.calibrate(source, it) }
        calibrated == calibrated.toSorted()

        where:
        source << ['diamond', 'hmmer', 'eggnog-mapper', 'foldseek', 'signalp', 'amrfinder', 'unknown']
    }

    def "DIAMOND calibration maps low pident to low probability"() {
        expect:
        table.calibrate('diamond', 0.3) < 0.35     // pident=30% is weak
        table.calibrate('diamond', 0.9) > 0.85     // pident=90% is strong
    }

    def "registering a custom curve overrides the default"() {
        given:
        def custom = new CalibrationTable()
        custom.register('diamond', 2.0, -1.0)

        when:
        def p = custom.calibrate('diamond', 0.5)

        then:
        // sigmoid(2*0.5 - 1) = sigmoid(0) = 0.5
        Math.abs(p - 0.5) < 1e-6
    }

    def "unknown sources fall back to default coefficients"() {
        expect:
        table.calibrate('never-heard-of-it', 0.5) > 0.0
        table.calibrate('never-heard-of-it', 0.5) < 1.0
    }
}
