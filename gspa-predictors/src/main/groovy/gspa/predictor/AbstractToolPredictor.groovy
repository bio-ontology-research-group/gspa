package gspa.predictor

import gspa.model.Annotation
import gspa.model.Protein
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Base class for predictors that wrap external command-line tools.
 * Handles common patterns: tool availability checking, temp directory management,
 * FASTA file writing, process execution with timeout, and output parsing.
 */
abstract class AbstractToolPredictor implements Predictor {

    protected final Logger log = LoggerFactory.getLogger(getClass())

    /** Working directory for temp files */
    File workDir

    /** Timeout in minutes for external tool execution */
    int timeoutMinutes = 120

    /** Additional environment variables */
    Map<String, String> environment = [:]

    /** Path to the executable (if not on PATH) */
    String executablePath

    /** The executable name or path */
    abstract String getExecutable()

    /** Build the command line for a given input/output */
    abstract List<String> buildCommand(File inputFasta, File outputDir)

    /** Parse the tool's output into annotations */
    abstract Map<String, List<Annotation>> parseOutput(File outputDir)

    @Override
    boolean isAvailable() {
        try {
            def exec = executablePath ?: executable
            def proc = [exec, '--version'].execute()
            proc.waitForOrKill(10000)
            return proc.exitValue() == 0
        } catch (Exception e) {
            // Also try --help or -h as some tools don't have --version
            try {
                def exec = executablePath ?: executable
                def proc = [exec, '--help'].execute()
                proc.waitForOrKill(10000)
                return proc.exitValue() <= 1 // some tools return 1 for --help
            } catch (Exception e2) {
                return false
            }
        }
    }

    @Override
    String getVersion() {
        try {
            def exec = executablePath ?: executable
            def proc = [exec, '--version'].execute()
            proc.waitForOrKill(10000)
            return proc.text.trim().readLines().first()
        } catch (Exception e) {
            return 'unknown'
        }
    }

    @Override
    List<Annotation> predict(Protein protein) {
        predictBatch([protein])[protein.id] ?: []
    }

    @Override
    Map<String, List<Annotation>> predictBatch(List<Protein> proteins) {
        if (proteins.isEmpty()) return [:]

        def tmpDir = createTempDir()
        try {
            // Write input FASTA
            def inputFasta = new File(tmpDir, 'input.faa')
            writeFasta(proteins, inputFasta)

            // Build and execute command
            def command = buildCommand(inputFasta, tmpDir)
            log.info("Running ${name}: ${command.join(' ')}")
            def result = execute(command, tmpDir)

            if (result.exitCode != 0) {
                log.error("${name} failed with exit code ${result.exitCode}: ${result.stderr}")
                return [:]
            }

            // Parse output
            return parseOutput(tmpDir)
        } finally {
            tmpDir.deleteDir()
        }
    }

    /**
     * Execute an external command.
     */
    protected ProcessResult execute(List<String> command, File workingDir = null) {
        def pb = new ProcessBuilder(command)
        if (workingDir) pb.directory(workingDir)
        pb.environment().putAll(environment)
        pb.redirectErrorStream(false)

        def proc = pb.start()

        // Capture output asynchronously
        def stdout = new StringBuilder()
        def stderr = new StringBuilder()
        def stdoutThread = Thread.start { proc.inputStream.eachLine { stdout.append(it).append('\n') } }
        def stderrThread = Thread.start { proc.errorStream.eachLine { stderr.append(it).append('\n') } }

        boolean finished = proc.waitFor(timeoutMinutes, java.util.concurrent.TimeUnit.MINUTES)
        if (!finished) {
            proc.destroyForcibly()
            return new ProcessResult(exitCode: -1, stdout: stdout.toString(),
                stderr: "Process timed out after ${timeoutMinutes} minutes")
        }

        stdoutThread.join(5000)
        stderrThread.join(5000)

        new ProcessResult(exitCode: proc.exitValue(), stdout: stdout.toString(), stderr: stderr.toString())
    }

    /**
     * Write proteins to a FASTA file.
     */
    protected void writeFasta(List<Protein> proteins, File output) {
        output.withWriter { writer ->
            proteins.each { protein ->
                writer.writeLine(">${protein.id}")
                writer.writeLine(protein.sequence)
            }
        }
    }

    private File createTempDir() {
        def dir = workDir ?
            new File(workDir, "gspa_${name}_${System.currentTimeMillis()}") :
            File.createTempDir("gspa_${name}_", '')
        dir.mkdirs()
        dir
    }
}

class ProcessResult {
    int exitCode
    String stdout
    String stderr
}
