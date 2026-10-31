package gspa.integration

import groovy.transform.CompileStatic

/**
 * Read-only change snapshots. Keep keys and numbers, not one seven-entry map per
 * candidate. List/Jackson readers still see the same ordered change records.
 * Records are materialized only while reading; do not use them to edit an audit.
 */
@CompileStatic
final class AuditChanges extends AbstractList<Map> implements RandomAccess {
    private final String[] keys
    private final double[] beforeValues
    private final double[] afterValues
    private final byte[] present

    private AuditChanges(int count) {
        keys = new String[count]
        beforeValues = new double[count]
        afterValues = new double[count]
        present = new byte[count]
    }

    static AuditChanges between(Map<String, Double> before, Map<String, Double> after) {
        int count = 0
        for (Map.Entry<String, Double> entry : before.entrySet()) {
            if (changed(entry.value, after.get(entry.key))) count++
        }
        for (Map.Entry<String, Double> entry : after.entrySet()) {
            if (!before.containsKey(entry.key)) count++
        }
        AuditChanges result = new AuditChanges(count)
        int index = 0
        for (Map.Entry<String, Double> entry : before.entrySet()) {
            Double now = after.get(entry.key)
            if (changed(entry.value, now)) result.store(index++, entry.key, entry.value, now)
        }
        for (Map.Entry<String, Double> entry : after.entrySet()) {
            if (!before.containsKey(entry.key)) result.store(index++, entry.key, null, entry.value)
        }
        return result
    }

    private static boolean changed(Double old, Double now) {
        return old == null || now == null || !(Math.abs(old - now) <= 1e-12d)
    }

    private void store(int index, String key, Double old, Double now) {
        // Validate at capture time, as the former eagerly materialized records did.
        if (IntegrationState.splitFunctionKey(key) == null) {
            throw new NullPointerException('Malformed audit function key: ' + key)
        }
        keys[index] = key
        if (old != null) { beforeValues[index] = old; present[index] = (byte) (present[index] | 1) }
        if (now != null) { afterValues[index] = now; present[index] = (byte) (present[index] | 2) }
    }

    @Override
    int size() { return keys.length }

    @Override
    Map get(int index) {
        Objects.checkIndex(index, keys.length)
        String key = keys[index]
        String[] parts = IntegrationState.splitFunctionKey(key)
        Double old = (present[index] & 1) != 0 ? Double.valueOf(beforeValues[index]) : null
        Double now = (present[index] & 2) != 0 ? Double.valueOf(afterValues[index]) : null
        return Collections.unmodifiableMap([function_key: key, protein_id: parts[0], function_type: parts[1], function_id: parts[2],
                change: old == null ? 'added' : (now == null ? 'removed' : (now > old ? 'increased' : 'decreased')),
                before_probability: old, after_probability: now])
    }
}
