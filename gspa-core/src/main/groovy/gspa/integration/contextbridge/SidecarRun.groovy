package gspa.integration.contextbridge

import com.fasterxml.jackson.databind.ObjectMapper

import java.security.MessageDigest
import java.util.concurrent.TimeUnit

/**
 * One bounded execution of a Python sidecar (metabolic experiments, community metabolism): a fresh
 * output, a hard wall-clock limit, the exit status and the output file are all checked, and a result
 * that cannot be trusted (non-zero exit, timeout, missing or stale output, unparseable output, wrong
 * schema) is a failure, never an empty success. The caller adds its own content checks on the parsed
 * result.
 *
 * <p>Credentials are never read here. A sidecar that needs an LLM endpoint reads its own environment
 * and key files; this class passes the environment through unchanged and puts no secret on a command
 * line or in a log.</p>
 */
class SidecarRun {

    static class Result {
        List<String> command
        int exitCode
        double seconds
        File log
        File output
        Map parsed
        String outputSha256
    }

    static class SidecarFailure extends IllegalStateException {
        SidecarFailure(String message) { super(message) }
    }

    private static final int LOG_TAIL = 1500

    /**
     * @param expectedSchema the {@code schema} value the parsed output must carry
     * @param output where the sidecar writes its JSON; any earlier file there is deleted first
     */
    static Result run(List<String> command, File output, File log, int timeoutSeconds, String expectedSchema) {
        output.parentFile?.mkdirs()
        log.parentFile?.mkdirs()
        output.delete()                       // fresh output: an earlier run's file can never pass for this one
        log.delete()
        long started = System.currentTimeMillis()
        Process proc
        try {
            proc = new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log).start()
        } catch (IOException e) {
            throw new SidecarFailure("cannot start ${command.take(2).join(' ')}: ${e.message}")
        }
        boolean finished = proc.waitFor(timeoutSeconds, TimeUnit.SECONDS)
        if (!finished) {
            proc.destroyForcibly()
            proc.waitFor(10, TimeUnit.SECONDS)
            output.delete()
            throw new SidecarFailure("${command[1]} exceeded its ${timeoutSeconds}s limit and was stopped; no result was used")
        }
        double seconds = (System.currentTimeMillis() - started) / 1000.0d
        int exit = proc.exitValue()
        if (exit != 0) {
            output.delete()
            throw new SidecarFailure("${command[1]} exited with status ${exit}: ${tail(log)}")
        }
        if (!output.isFile()) throw new SidecarFailure("${command[1]} exited 0 but wrote no output ${output}: ${tail(log)}")
        if (output.lastModified() < started - 2000L) {
            output.delete()
            throw new SidecarFailure("output ${output} is older than this run (stale); refusing it")
        }
        Map parsed
        try {
            parsed = new ObjectMapper().readValue(output, Map)
        } catch (Exception e) {
            output.delete()
            throw new SidecarFailure("output ${output} is not valid JSON: ${e.message}")
        }
        if (parsed.schema != expectedSchema) {
            output.delete()
            throw new SidecarFailure("output ${output} has schema '${parsed.schema}', expected '${expectedSchema}'")
        }
        new Result(command: command, exitCode: exit, seconds: seconds, log: log, output: output, parsed: parsed,
            outputSha256: sha256(output))
    }

    static String tail(File log) {
        if (!log.isFile()) return '(no log)'
        String text = log.getText('UTF-8').trim()
        text.length() <= LOG_TAIL ? text : '...' + text.substring(text.length() - LOG_TAIL)
    }

    static String sha256(File f) {
        MessageDigest md = MessageDigest.getInstance('SHA-256')
        f.withInputStream { is ->
            byte[] buf = new byte[65536]
            int n
            while ((n = is.read(buf)) > 0) md.update(buf, 0, n)
        }
        md.digest().collect { String.format('%02x', it & 0xff) }.join()
    }
}
