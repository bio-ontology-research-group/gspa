package gspa.cli

import picocli.CommandLine.IVersionProvider

/**
 * Reads the version baked into the jar at build time from
 * {@code /version.properties}. Single source of truth so the picocli
 * version string can never drift from the Gradle project version.
 */
class VersionProvider implements IVersionProvider {

    @Override
    String[] getVersion() {
        def props = new Properties()
        def stream = getClass().getResourceAsStream('/version.properties')
        if (stream != null) {
            stream.withCloseable { props.load(it) }
        }
        def v = props.getProperty('version', 'unknown')
        ["gspa ${v}".toString()] as String[]
    }
}
