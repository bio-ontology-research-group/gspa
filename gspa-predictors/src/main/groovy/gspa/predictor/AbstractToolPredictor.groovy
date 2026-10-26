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
abstract class AbstractToolPredictor implements Predictor, StatusReportingPredictor {

    protected final Logger log = LoggerFactory.getLogger(getClass())

    /** Completion status of the most recent run (null before the first run). */
    volatile PredictorRunStatus lastRunStatus

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

    /** Generic reason; subclasses with configuration prerequisites refine it. */
    @Override
    String unavailableReason() {
        isAvailable() ? null : "executable '${executablePath ?: executable}' not found or not runnable"
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
        lastRunStatus = null   // never report a stale status for this run
        if (proteins.isEmpty()) {
            lastRunStatus = PredictorRunStatus.completed(name, 0)
            return [:]
        }

        def tmpDir = createTempDir()
        long t0 = System.currentTimeMillis()
        try {
            // Write input FASTA
            def inputFasta = new File(tmpDir, 'input.faa')
            writeFasta(proteins, inputFasta)

            def command = buildCommand(inputFasta, tmpDir)
            def result = runOrFail(command, tmpDir, t0)
            if (result == null) return [:]

            def parsed = parseGuarded(t0) { parseOutput(tmpDir) }
            recordStatus(PredictorRunStatus.completed(name, countAnnotations(parsed), result.exitCode, successNote()), t0)
            return parsed
        } catch (Exception e) {
            if (lastRunStatus == null || lastRunStatus.outcome != PredictorOutcome.FAILED) {
                recordStatus(PredictorRunStatus.failed(name, "${e.class.simpleName}: ${e.message}"), t0)
            }
            throw e
        } finally {
            tmpDir.deleteDir()
        }
    }

    // ---- shared run/contract helpers (used by every adapter, incl. genome-level overrides) ----

    /**
     * Interpretation caveat attached to SUCCESS / SUCCESS_ZERO statuses (e.g. what an
     * empty result does not mean for this tool); null when there is none.
     */
    protected String successNote() { null }

    protected void recordStatus(PredictorRunStatus status, long startMillis) {
        status.elapsedMillis = System.currentTimeMillis() - startMillis
        lastRunStatus = status
    }

    /**
     * Execute the command; a non-zero exit (including timeout, exit -1) records FAILED
     * with the stderr tail and returns null. A launch failure (IOException) propagates.
     */
    protected ProcessResult runOrFail(List<String> command, File workingDir, long startMillis) {
        log.info("Running ${name}: ${command.join(' ')}")
        def result = execute(command, workingDir)
        if (result.exitCode != 0) {
            log.error("${name} failed with exit code ${result.exitCode}: ${result.stderr}")
            recordStatus(PredictorRunStatus.failed(name,
                "exit code ${result.exitCode}: ${result.stderr?.trim()?.take(500)}", result.exitCode), startMillis)
            return null
        }
        result
    }

    /**
     * Run a parse step; any exception (schema contract violation or otherwise)
     * is recorded as FAILED with its message and rethrown.
     */
    protected <T> T parseGuarded(long startMillis, Closure<T> body) {
        try {
            body.call()
        } catch (Exception e) {
            log.error("${name} output rejected: ${e.message}")
            recordStatus(PredictorRunStatus.failed(name, "${e.class.simpleName}: ${e.message}"), startMillis)
            throw e
        }
    }

    /**
     * Return a file the tool is documented to write on every successful run, or throw
     * {@link ToolOutputContractException}. {@code allowEmpty} must be true only for tools
     * whose schema has no header/comment line and which write a zero-byte file when
     * nothing matched (DIAMOND, MMseqs2 convertalis, FoldSeek, InterProScan TSV).
     */
    protected File requireOutput(File dir, String fileName, String why, boolean allowEmpty = false) {
        def f = new File(dir, fileName)
        if (!f.isFile()) {
            throw new ToolOutputContractException(
                "${name}: expected output '${fileName}' is missing after a clean exit (${why})")
        }
        if (!allowEmpty && f.length() == 0) {
            throw new ToolOutputContractException(
                "${name}: output '${fileName}' is empty but this tool always writes a header/comment line (${why})")
        }
        f
    }

    /** Result of exporting contig sequences for a nucleotide-input tool. */
    static class ContigExport {
        List<String> written = []
        /** Contigs with no sequence loaded: not analysed. */
        List<String> skipped = []
    }

    /** Write contigs that carry sequence as FASTA; contigs without are reported, not silently dropped. */
    protected static ContigExport writeContigFasta(gspa.model.Genome genome, File out) {
        def export = new ContigExport()
        out.withWriter { w ->
            genome.contigs.each { contig ->
                if (contig.sequence) {
                    gspa.io.FastaWriter.writeEntry(w, contig.id, contig.sequence)
                    export.written << contig.id
                } else {
                    export.skipped << contig.id
                }
            }
        }
        export
    }

    /**
     * Template for nucleotide-input genome tools (MinCED, antiSMASH, ...): export contigs,
     * run, parse, and record a typed status. Contigs without sequence make the run
     * PARTIAL (never a full SUCCESS/SUCCESS_ZERO); no sequence at all is FAILED.
     *
     * @param outputDirOf maps the scratch dir to the directory passed to {@link #buildCommand}
     *                    and parsed afterwards (default: the scratch dir itself)
     */
    protected Map<String, List<Annotation>> runGenomeTool(gspa.model.Genome genome, String inputName,
                                                          Closure<File> outputDirOf = null) {
        lastRunStatus = null
        long t0 = System.currentTimeMillis()
        def tmpDir = File.createTempDir("gspa_${name}_", '')
        try {
            def input = new File(tmpDir, inputName)
            def export = writeContigFasta(genome, input)
            if (export.written.isEmpty()) {
                recordStatus(PredictorRunStatus.failed(name,
                    'no contig carries nucleotide sequence; a nucleotide tool cannot run on protein-only input'), t0)
                return [:]
            }
            def outDir = outputDirOf ? outputDirOf.call(tmpDir) : tmpDir
            def result = runOrFail(buildCommand(input, outDir), tmpDir, t0)
            if (result == null) return [:]
            def parsed = parseGuarded(t0) { parseOutput(outDir) }
            int n = countAnnotations(parsed)
            if (export.skipped) {
                recordStatus(PredictorRunStatus.partial(name, n, describeSkipped(export, genome.contigs.size()),
                    export.skipped, result.exitCode), t0)
            } else {
                recordStatus(PredictorRunStatus.completed(name, n, result.exitCode, successNote()), t0)
            }
            return parsed
        } catch (Exception e) {
            if (lastRunStatus == null || lastRunStatus.outcome != PredictorOutcome.FAILED) {
                recordStatus(PredictorRunStatus.failed(name, "${e.class.simpleName}: ${e.message}"), t0)
            }
            throw e
        } finally {
            tmpDir.deleteDir()
        }
    }

    protected static String describeSkipped(ContigExport export, int total) {
        "${export.skipped.size()} of ${total} contigs have no nucleotide sequence and were not analysed: " +
            export.skipped.take(10).join(', ') + (export.skipped.size() > 10 ? ', ...' : '')
    }

    protected static int countAnnotations(Map<String, List<Annotation>> results) {
        results ? (results.values().sum { it.size() } as Integer ?: 0) : 0
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
