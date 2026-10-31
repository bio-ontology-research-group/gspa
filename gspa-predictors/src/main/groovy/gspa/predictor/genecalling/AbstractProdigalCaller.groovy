package gspa.predictor.genecalling

import gspa.io.GeneModelReader
import gspa.model.*
import gspa.predictor.PredictorRunStatus
import gspa.predictor.StatusReportingPredictor
import org.slf4j.LoggerFactory
import java.util.concurrent.TimeUnit

/** Shared process/output contract for single-row Prodigal-family gene models.
 * Instances retain the most recent status and are used sequentially by a pipeline.
 */
abstract class AbstractProdigalCaller implements GeneCaller, StatusReportingPredictor {
    boolean metaMode = false
    int translationTable = 11
    String executablePath
    long timeoutMillis = TimeUnit.MINUTES.toMillis(30)
    long availabilityTimeoutMillis = 10000L
    PredictorRunStatus lastRunStatus
    protected String availabilityDetail

    protected abstract String versionArgument()

    @Override
    String unavailableReason() { availabilityDetail }

    @Override
    boolean isAvailable() {
        File dir = null
        try {
            dir = File.createTempDir("gspa_${name}_version_", '')
            Map result = execute([executablePath ?: name, versionArgument()], dir, availabilityTimeoutMillis)
            boolean valid = !result.unavailable && !result.timedOut && result.exitCode == 0 &&
                result.output.toLowerCase(Locale.ROOT).contains(name.toLowerCase(Locale.ROOT))
            availabilityDetail = valid ? null : "version probe failed: ${result.detail ?: result.output}"
            return valid
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            availabilityDetail = 'version probe interrupted'
            return false
        } catch (Exception e) {
            availabilityDetail = "${e.class.simpleName}: ${e.message}"
            return false
        } finally {
            dir?.deleteDir()
        }
    }

    @Override
    void callGenes(Genome genome) {
        lastRunStatus = null
        long started = System.currentTimeMillis()
        String statusName = "gene-caller:${name}"
        File dir = null
        Integer exitCode = null
        try {
            if (genome.proteinCount || genome.contigs.any { it.features.any { it.type == FeatureType.CDS } }) {
                throw new IllegalArgumentException('gene calling requires input without existing CDS or proteins')
            }
            List<String> ids = genome.contigs*.id
            if (ids.any { !it || it =~ /\s/ } || ids.toSet().size() != ids.size()) {
                throw new IllegalArgumentException('empty, duplicate or whitespace-containing contig identifier')
            }
            List<Contig> supplied = genome.contigs.findAll { it.sequence }
            List<String> skipped = genome.contigs.findAll { !it.sequence }*.id
            if (!supplied) throw new IllegalArgumentException('no contig carries nucleotide sequence; gene caller was not executed')
            dir = File.createTempDir("gspa_${name}_", '')
            File input = new File(dir, 'input.fna')
            File gff = new File(dir, 'genes.gff')
            File proteins = new File(dir, 'proteins.faa')
            input.withWriter('UTF-8') { writer ->
                supplied.each { writer.writeLine(">${it.id}"); writer.writeLine(it.sequence) }
            }
            List<String> command = [executablePath ?: name, '-i', input.absolutePath,
                '-o', gff.absolutePath, '-a', proteins.absolutePath, '-f', 'gff',
                '-g', translationTable.toString()]
            if (metaMode) command.addAll(['-p', 'meta'])
            Map result = execute(command, dir, timeoutMillis)
            exitCode = result.exitCode as Integer
            if (result.unavailable) {
                lastRunStatus = PredictorRunStatus.unavailable(statusName, result.detail as String)
                return
            }
            if (result.timedOut || exitCode != 0) {
                lastRunStatus = PredictorRunStatus.failed(statusName,
                    "${result.detail ?: 'native process failed'}; ${result.output}", exitCode)
                return
            }
            if (!gff.isFile() || !gff.length() || !proteins.isFile()) {
                throw new IllegalArgumentException('successful native exit but required GFF/protein output is missing or GFF is empty')
            }
            List<String> order = validateNativeCoverage(gff, supplied)
            Genome parsed = GeneModelReader.read(input, gff, proteins)
            Map<String, Contig> targets = genome.contigs.collectEntries { [(it.id): it] }
            Map<String, Contig> models = parsed.contigs.collectEntries { [(it.id): it] }
            int count = 0
            // Validate everything before mutation. Preserve established contig/CDS order
            // and public IDs; FASTA order must not influence the sequence-to-CDS join.
            order.each { id ->
                Contig model = models[id]
                Map<Feature, Protein> byFeature = new IdentityHashMap<>()
                model.proteins.each { byFeature[it.sourceFeature] = it }
                model.features.sort().each { feature ->
                    Protein protein = byFeature[feature]
                    String tag = "${genome.id}_${String.format('%05d', ++count)}"
                    feature.source = name
                    feature.id = tag
                    feature.attributes['locus_tag'] = [tag]
                    protein.id = tag
                    targets[id].addFeature(feature)
                    targets[id].addProtein(protein)
                }
                targets[id].sortFeatures()
            }
            lastRunStatus = skipped ? PredictorRunStatus.partial(statusName, count,
                'contigs without nucleotide sequence were not analysed', skipped, exitCode) :
                PredictorRunStatus.completed(statusName, count, exitCode,
                    count ? null : 'native caller completed with no CDS calls; this does not establish biological gene absence')
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            lastRunStatus = PredictorRunStatus.failed(statusName, 'gene calling interrupted', exitCode)
        } catch (Exception e) {
            lastRunStatus = PredictorRunStatus.failed(statusName, "${e.class.simpleName}: ${e.message}", exitCode)
        } finally {
            if (lastRunStatus != null) lastRunStatus.elapsedMillis = System.currentTimeMillis() - started
            if (lastRunStatus?.detail) LoggerFactory.getLogger(getClass()).info(lastRunStatus.detail)
            dir?.deleteDir()
        }
    }

    /** Native comments bind all output, including zero-CDS contigs, to submitted input. */
    protected static List<String> validateNativeCoverage(File gff, List<Contig> supplied) {
        Map<String, Integer> expected = supplied.collectEntries { [(it.id): it.length] }
        Map<String, Integer> observed = [:]
        boolean version = false
        gff.eachLine('UTF-8') { line ->
            if (line ==~ /##gff-version\s+3\s*/) version = true
            if (!line.startsWith('# Sequence Data:')) return
            def match = line =~ /^# Sequence Data: seqnum=\d+;seqlen=(\d+);seqhdr="([^"]*)"\s*$/
            if (!match.matches()) throw new IllegalArgumentException('malformed native contig coverage record')
            String id = match[0][2].split(/\s+/)[0]
            int length = match[0][1].toInteger()
            if (!expected.containsKey(id) || expected[id] != length || observed.containsKey(id)) {
                throw new IllegalArgumentException('duplicate, unknown or inconsistent native contig coverage')
            }
            observed[id] = length
        }
        if (!version || observed != expected) throw new IllegalArgumentException('native GFF does not cover every submitted contig')
        new ArrayList<String>(observed.keySet())
    }

    /** Disk output lets waitFor enforce the deadline even when a pipe never closes. */
    protected static Map execute(List<String> command, File directory, long timeout) {
        if (timeout < 1) throw new IllegalArgumentException('process timeout must be positive')
        File log = new File(directory, 'process.log')
        Process process
        try {
            process = new ProcessBuilder(command).directory(directory).redirectErrorStream(true).redirectOutput(log).start()
        } catch (IOException e) {
            return [unavailable: true, detail: "cannot start executable: ${e.message}", output: '']
        }
        try {
            boolean finished = process.waitFor(timeout, TimeUnit.MILLISECONDS)
            if (!finished) {
                terminate(process)
                return [timedOut: true, exitCode: null, detail: "process timed out after ${timeout} ms", output: tail(log)]
            }
            [exitCode: process.exitValue(), output: tail(log)]
        } finally {
            if (process.isAlive()) terminate(process)
        }
    }

    protected static String tail(File file) {
        if (!file.isFile()) return ''
        new RandomAccessFile(file, 'r').withCloseable { stream ->
            long start = Math.max(0L, stream.length() - 8192L)
            stream.seek(start)
            byte[] bytes = new byte[(int) (stream.length() - start)]
            stream.readFully(bytes)
            (start ? '[diagnostic tail] ' : '') + new String(bytes, java.nio.charset.StandardCharsets.UTF_8)
        }
    }

    protected static void terminate(Process process) {
        List<ProcessHandle> children
        process.toHandle().descendants().withCloseable { children = it.toList() }
        // Kill children first while their parent can still reap them.
        children.reverseEach { it.destroyForcibly() }
        boolean interrupted = false
        try { process.waitFor(200, TimeUnit.MILLISECONDS) } catch (InterruptedException ignored) { interrupted = true }
        if (process.isAlive()) process.destroyForcibly()
        try { process.waitFor(1000, TimeUnit.MILLISECONDS) } catch (InterruptedException ignored) { interrupted = true }
        if (interrupted) Thread.currentThread().interrupt()
    }
}
