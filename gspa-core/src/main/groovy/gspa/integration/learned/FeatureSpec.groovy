package gspa.integration.learned

import groovy.transform.Canonical

/** One documented column of a {@link FeatureContract}. */
@Canonical
class FeatureSpec {
    String name
    FeatureGroup group
    /** What the value is, its units/range, and what value stands in when it is unavailable. */
    String description
}
