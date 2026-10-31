package gspa.integration.contextbridge

import com.fasterxml.jackson.databind.ObjectMapper

import java.nio.file.Files
import java.nio.file.StandardCopyOption
import java.security.MessageDigest
import java.util.concurrent.TimeUnit

/**
 * One bounded execution of a Python sidecar (metabolic experiments, community metabolism).
 *
 * <p><b>Fresh output by construction.</b> The sidecar is pointed at a PRIVATE, newly created directory, never at the
 * published path. Only after the exit status, schema and the caller's content checks have passed is the file moved
 * to the published path, replacing whatever was there; a result that fails any check never reaches it, and an
 * earlier file there that cannot be removed is an error, not something that could pass for this run.</p>
 *
 * <p><b>Bounded.</b> A hard wall-clock limit; on expiry the process AND its descendants are killed, so nothing keeps
 * running or writes late (late writes could only land in the private directory, which is removed).</p>
 *
 * <p>Credentials are never read here; the environment is passed through unchanged and no secret is put on a
 * command line or in a log.</p>
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
     * @param command the command line; the value after {@code --out} is replaced by the private path
     * @param output the published path
     * @param verify optional content check on the parsed result; throw {@link SidecarFailure} to refuse it
     */
    static Result run(List<String> command, File output, File log, int timeoutSeconds, String expectedSchema,
                      Closure verify = null) {
        output.parentFile?.mkdirs()
        log.parentFile?.mkdirs()
        if (output.exists() && !output.delete()) throw new SidecarFailure("cannot remove the earlier output ${output}; refusing to run so it cannot be mistaken for this run")
        log.delete()
        int at = command.indexOf('--out')
        if (at < 0 || at + 1 >= command.size()) throw new IllegalArgumentException('the sidecar command has no --out argument')
        File privateDir = Files.createTempDirectory(output.parentFile.toPath(), '.run-').toFile()
        File privateOut = new File(privateDir, 'result.json')
        List<String> cmd = new ArrayList<>(command)
        cmd.set(at + 1, privateOut.absolutePath)
        try {
            long started = System.currentTimeMillis()
            Process proc
            try {
                proc = new ProcessBuilder(cmd).redirectErrorStream(true).redirectOutput(log).start()
            } catch (IOException e) {
                throw new SidecarFailure("cannot start ${cmd.take(2).join(' ')}: ${e.message}")
            }
            boolean finished = proc.waitFor(timeoutSeconds, TimeUnit.SECONDS)
            if (!finished) {
                proc.descendants().forEach { it.destroyForcibly() }
                proc.destroyForcibly()
                proc.waitFor(10, TimeUnit.SECONDS)
                throw new SidecarFailure("${cmd[1]} exceeded its ${timeoutSeconds}s limit and was stopped with its child processes; no result was used")
            }
            double seconds = (System.currentTimeMillis() - started) / 1000.0d
            int exit = proc.exitValue()
            if (exit != 0) throw new SidecarFailure("${cmd[1]} exited with status ${exit}: ${tail(log)}")
            if (!privateOut.isFile()) throw new SidecarFailure("${cmd[1]} exited 0 but wrote no output: ${tail(log)}")
            Map parsed
            try {
                parsed = new ObjectMapper().readValue(privateOut, Map)
            } catch (Exception e) {
                throw new SidecarFailure("the output is not valid JSON: ${e.message}")
            }
            if (parsed.schema != expectedSchema) throw new SidecarFailure("the output has schema '${parsed.schema}', expected '${expectedSchema}'")
            if (verify != null) verify.call(parsed)
            Files.move(privateOut.toPath(), output.toPath(), StandardCopyOption.REPLACE_EXISTING)
            new Result(command: cmd, exitCode: exit, seconds: seconds, log: log, output: output, parsed: parsed,
                outputSha256: sha256(output))
        } finally {
            privateDir.deleteDir()
        }
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
