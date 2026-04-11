package gspa.predictor.taxonomy

import gspa.model.AssemblyInfo
import gspa.model.Genome
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Runs CheckM2 to estimate MAG completeness and contamination.
 * Parses the quality_report.tsv output.
 */
class CheckM2Runner {

    private static final Logger log = LoggerFactory.getLogger(CheckM2Runner)

    String executablePath
    int threads = Runtime.runtime.availableProcessors()

    boolean isAvailable() {
        try {
            def proc = [executablePath ?: 'checkm2', '--help'].execute()
            proc.waitForOrKill(10000)
            return proc.exitValue() <= 1
        } catch (Exception e) {
            return false
        }
    }

    /**
     * Run CheckM2 on genome FASTA files and return quality estimates.
     */
    Map<String, AssemblyInfo> run(File inputDir, File outputDir) {
        outputDir.mkdirs()

        def cmd = [
            executablePath ?: 'checkm2',
            'predict',
            '--input', inputDir.absolutePath,
            '--output-directory', outputDir.absolutePath,
            '--threads', threads.toString(),
            '-x', 'fna',
        ]

        log.info("Running CheckM2: ${cmd.join(' ')}")
        def proc = new ProcessBuilder(cmd).start()
        proc.waitFor(120, java.util.concurrent.TimeUnit.MINUTES)

        if (proc.exitValue() != 0) {
            log.error("CheckM2 failed: ${proc.errorStream.text}")
            return [:]
        }

        parseOutput(new File(outputDir, 'quality_report.tsv'))
    }

    /**
     * Parse CheckM2 quality_report.tsv output.
     * Format: Name\tCompleteness\tContamination\t...
     */
    static Map<String, AssemblyInfo> parseOutput(File reportFile) {
        if (!reportFile.exists()) return [:]

        Map<String, AssemblyInfo> results = [:]

        reportFile.eachLine { line ->
            if (line.startsWith('Name')) return  // header
            def fields = line.split('\t')
            if (fields.length < 3) return

            String name = fields[0].trim()
            double completeness = fields[1].trim() as double
            double contamination = fields[2].trim() as double

            results[name] = new AssemblyInfo(
                completeness: completeness,
                contamination: contamination,
                assemblySource: 'MAG'
            )
        }

        log.info("CheckM2: parsed ${results.size()} genome quality estimates")
        results
    }

    /**
     * Parse and assign CheckM2 results to a single genome.
     */
    static void assignQuality(Genome genome, File reportFile) {
        def results = parseOutput(reportFile)
        // Try matching by genome ID
        def info = results[genome.id] ?: results.values().first()
        if (info) {
            genome.assemblyInfo = info
            genome.mag = true
            log.info("Assigned CheckM2 quality to ${genome.id}: " +
                "completeness=${info.completeness}%, contamination=${info.contamination}%")
        }
    }
}
