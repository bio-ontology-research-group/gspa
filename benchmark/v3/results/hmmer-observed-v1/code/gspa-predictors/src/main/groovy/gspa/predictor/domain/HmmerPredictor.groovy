package gspa.predictor.domain

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

/**
 * HMMER predictor for direct Pfam domain searches.
 * Lighter-weight alternative to InterProScan when only Pfam is needed.
 * Uses hmmsearch against Pfam-A.hmm.
 */
class HmmerPredictor extends AbstractToolPredictor {

    /** Path to Pfam-A HMM database */
    String hmmDatabase

    /** E-value threshold for domain hits */
    double domainEvalue = 1e-5

    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'hmmer' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.PFAM] as Set
    }

    @Override
    String getExecutable() { 'hmmsearch' }

    /** hmmsearch needs a Pfam-A HMM file; without one the tool cannot run. */
    @Override
    boolean isAvailable() {
        hmmDatabase && new File(hmmDatabase).isFile() && probeVersion() != null
    }

    /** hmmsearch has no --version option; -h reports both program and version. */
    private String probeVersion() {
        try {
            def builder = new ProcessBuilder([executablePath ?: executable, '-h'])
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

    @Override
    String getVersion() { probeVersion() ?: 'unknown' }

    @Override
    String unavailableReason() {
        if (!hmmDatabase) return 'no HMM database configured (hmmDatabase / GSPA_PFAM_HMM)'
        if (!new File(hmmDatabase).isFile()) return "HMM database not found: ${hmmDatabase}"
        super.unavailableReason()
    }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (!hmmDatabase) {
            throw new IllegalStateException("${name}: hmmDatabase is unset (path to Pfam-A.hmm)")
        }
        def outputFile = new File(outputDir, 'hmmer_domtbl.tsv')
        [
            executablePath ?: executable,
            '--domtblout', outputFile.absolutePath,
            '--noali',
            '-E', '1',
            '--domE', domainEvalue.toString(),
            '--cpu', threads.toString(),
            hmmDatabase,
            inputFasta.absolutePath,
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // hmmsearch --domtblout always writes '#' header/footer lines, even with no hit,
        // so a missing or empty file is a failed run, not "no domains".
        def outputFile = requireOutput(outputDir, 'hmmer_domtbl.tsv',
            'hmmsearch --domtblout writes comment header/footer on every run')

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        outputFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split(/\s+/)
            if (fields.length < 22) {
                throw new ToolOutputContractException(
                    "hmmer domtblout row has ${fields.length} columns, the format has at least 22: ${line.take(80)}")
            }

            // hmmsearch domtblout format:
            // [0]=target name (protein), [1]=target acc, [2]=tlen
            // [3]=query name (HMM), [4]=query acc (Pfam ID), [5]=qlen
            // [6..8]=full seq E/score/bias, [9..10]=domain #/of
            // [11..13]=domain cE/iE/score, [14]=bias
            String proteinId = fields[0]
            String pfamAcc = fields[4]       // PF00001.21 (query accession)
            String pfamName = fields[3]      // HMM name
            double domEvalue
            double domScore
            try {
                domEvalue = Double.parseDouble(fields[12])
                domScore = Double.parseDouble(fields[13])
            } catch (NumberFormatException ignored) {
                throw new ToolOutputContractException("hmmer: invalid domain E-value/bit score for ${proteinId}")
            }
            if (!Double.isFinite(domEvalue) || domEvalue < 0.0d || !Double.isFinite(domScore)) {
                throw new ToolOutputContractException("hmmer: non-finite statistic or negative domain E-value for ${proteinId}")
            }
            if (pfamAcc == '-') {
                throw new ToolOutputContractException("hmmer: model '${pfamName}' lacks a query accession for ${proteinId}")
            }

            // Normalize accession (remove version)
            String acc = pfamAcc.replaceAll(/\.\d+$/, '')

            double score = Math.max(0.0d, Math.min(1.0d, domScore / 100.0d))

            results[proteinId] << new Annotation(
                type: AnnotationType.PFAM,
                value: acc,
                source: name,
                score: score,
                metadata: [
                    name: pfamName,
                    domainEvalue: domEvalue,
                    domainScore: domScore,
                ]
            )
        }
        results
    }
}
