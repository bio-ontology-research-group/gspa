package gspa.integration

import com.fasterxml.jackson.databind.ObjectMapper
import spock.lang.Specification

class AuditChangesSpec extends Specification {
    // Frozen previous implementation: exact record/order/JSON compatibility oracle.
    private static List<Map> legacy(Map<String, Double> before, Map<String, Double> after) {
        Set<String> keys = new LinkedHashSet<>(before.keySet())
        keys.addAll(after.keySet())
        List<Map> out = []
        keys.each { key ->
            Double old = before[key], now = after[key]
            if (old != null && now != null && Math.abs(old - now) <= 1e-12d) return
            String[] parts = IntegrationState.splitFunctionKey(key)
            out << [function_key: key, protein_id: parts[0], function_type: parts[1], function_id: parts[2],
                    change: old == null ? 'added' : (now == null ? 'removed' : (now > old ? 'increased' : 'decreased')),
                    before_probability: old, after_probability: now]
        }
        out
    }

    def 'records and pretty JSON preserve additions removals deltas tolerance and key order'() {
        given:
        Map before = ['sp|P00956|SYI_ECOLI|EC|6.1.1.9': 0.8d,
            'removed|EC|1.1.1.1': 0d, 'unchanged|EC|2.2.2.2': 0.4d,
            'increased|GO|GO:0008150': 0.2d, 'tolerance|EC|3.3.3.3': 0.5d]
        Map after = ['added|EC|1.1.1.n11': 0d, 'tolerance|EC|3.3.3.3': 0.5d + 1e-13d,
            'increased|GO|GO:0008150': 0.7d, 'unchanged|EC|2.2.2.2': 0.4d,
            'sp|P00956|SYI_ECOLI|EC|6.1.1.9': 0.3d]
        def writer = new ObjectMapper().writerWithDefaultPrettyPrinter()

        expect:
        AuditChanges.between(before, after) == legacy(before, after)
        writer.writeValueAsString(AuditChanges.between(before, after)) == writer.writeValueAsString(legacy(before, after))
        AuditChanges.between(before, after)*.change == ['decreased', 'removed', 'increased', 'added']
        AuditChanges.between(before, after)[0].protein_id == 'sp|P00956|SYI_ECOLI'
    }

    def 'captured scores survive later mutation and release of source maps'() {
        given:
        def before = ['gene|EC|1.1.1.1': 0.9d]
        def after = ['gene|EC|1.1.1.1': 0.2d]
        def expected = legacy(before, after)
        def changes = AuditChanges.between(before, after)

        when:
        before.clear()
        after['gene|EC|1.1.1.1'] = 1d

        then:
        changes == expected
        changes == expected // repeated iteration never consumes the snapshot
    }

    def 'full seeded panel matches legacy with empty and disjoint populations'() {
        given:
        Random random = new Random(83)
        Map before = [:], after = [:]
        (0..<10000).each { i ->
            String key = "protein.${i}|EC|1.1.1.${i + 1}"
            if (i % 3 != 0) before[key] = random.nextDouble()
            if (i % 5 != 0) after[key] = i % 2 == 0 && before.containsKey(key) ? before[key] : random.nextDouble()
        }
        def writer = new ObjectMapper().writerWithDefaultPrettyPrinter()

        expect:
        [[before, after], [[:], after], [before, [:]], [[:], [:]]].every { pair ->
            writer.writeValueAsString(AuditChanges.between(pair[0], pair[1])) == writer.writeValueAsString(legacy(pair[0], pair[1]))
        }
    }

    def 'presence bits preserve null signed zero and nonfinite legacy values'() {
        given:
        Map before = ['null|EC|1.1.1.1': null, 'nan|EC|1.1.1.2': Double.NaN,
            'inf|EC|1.1.1.3': Double.POSITIVE_INFINITY, 'zero|EC|1.1.1.4': -0d,
            'nanAfter|EC|1.1.1.5': 0d, 'infAfter|EC|1.1.1.6': 1d,
            'boundary|EC|1.1.1.7': 0d, 'signedZero|EC|1.1.1.8': -0d]
        Map after = ['null|EC|1.1.1.1': null, 'nan|EC|1.1.1.2': 0d,
            'inf|EC|1.1.1.3': Double.POSITIVE_INFINITY, 'zero|EC|1.1.1.4': 1d,
            'nanAfter|EC|1.1.1.5': Double.NaN, 'infAfter|EC|1.1.1.6': Double.NEGATIVE_INFINITY,
            'boundary|EC|1.1.1.7': 1e-12d, 'signedZero|EC|1.1.1.8': 0d]
        def writer = new ObjectMapper().writerWithDefaultPrettyPrinter()

        expect:
        writer.writeValueAsString(AuditChanges.between(before, after)) == writer.writeValueAsString(legacy(before, after))
    }

    def 'malformed changed identity fails during capture'() {
        when:
        AuditChanges.between([:], ['malformed': 0.5d])

        then:
        thrown(NullPointerException)
    }

    def 'removed malformed identity also fails but unchanged malformed identity is ignored'() {
        expect:
        AuditChanges.between(['malformed': 0.5d], ['malformed': 0.5d]).empty

        when:
        AuditChanges.between(['malformed': 0.5d], [:])

        then:
        thrown(NullPointerException)
    }

    def 'snapshot list and records reject mutation and check indices'() {
        given:
        def changes = AuditChanges.between([:], ['gene|EC|1.1.1.1': 0d])

        when:
        changes.add([:])
        then:
        thrown(UnsupportedOperationException)

        when:
        changes[0].after_probability = 1d
        then:
        thrown(UnsupportedOperationException)

        when:
        changes.get(1)
        then:
        thrown(IndexOutOfBoundsException)
    }
}
