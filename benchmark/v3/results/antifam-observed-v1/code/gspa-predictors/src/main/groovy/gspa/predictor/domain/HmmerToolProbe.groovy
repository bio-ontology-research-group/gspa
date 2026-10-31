package gspa.predictor.domain

/** Shared executable probe for adapters invoking hmmsearch, which has no --version. */
final class HmmerToolProbe {
    static String version(String executable, Map<String, String> environment) {
        try {
            def builder = new ProcessBuilder([executable, '-h'])
            builder.environment().putAll(environment)
            def process = builder.redirectErrorStream(true).start()
            process.waitForOrKill(10000)
            if (process.exitValue() != 0) return null
            List<String> lines = process.text.readLines()
            if (!lines.any { it.startsWith('# hmmsearch ::') }) return null
            String version = lines.find { it.startsWith('# HMMER ') }
            version != null ? version.substring(2).trim() : null
        } catch (Exception ignored) {
            null
        }
    }
}
