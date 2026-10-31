package gspa.predictor.domain

import gspa.model.Annotation
import gspa.integration.EvidenceType
import groovy.json.JsonOutput
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

/**
 * AntiFam scanner for spurious ORFs and local gene-model errors.
 * A hit can indicate a wrong boundary (including the prfB frameshift example),
 * not whole-protein inactivity. PSEUDOGENE is retained as the legacy output
 * namespace; GENE_MODEL_QUALITY makes the non-functional evidence semantics
 * explicit. No score emitted here is a calibrated probability of inactivity.
 * See Eberhardt et al. (2012), doi:10.1093/database/bas003.
 *
 * <p>AntiFam HMMs ship with built-in gathering thresholds (GA), so the
 * canonical invocation is {@code hmmsearch --cut_ga …}. We honour that
 * by default; users wanting a stricter or looser cut can override via
 * {@link #domainEvalue} which switches to E-value gating.</p>
 *
 * <p>Database: download from
 * {@code https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/AntiFam.tar.gz}
 * and point {@link #hmmDatabase} at the concatenated {@code AntiFam.hmm}.</p>
 */
class AntiFamPredictor extends AbstractToolPredictor {

    /** Path to the AntiFam HMM database (e.g. AntiFam.hmm) */
    String hmmDatabase

    /**
     * Use the HMM-built gathering threshold (--cut_ga). True is the
     * recommended setting for AntiFam; flip to false to gate by
     * {@link #domainEvalue} instead.
     */
    boolean useCutGa = true

    /** E-value threshold for domain hits when {@link #useCutGa} is false. */
    double domainEvalue = 1e-5

    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'antifam' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.PSEUDOGENE] as Set
    }

    @Override
    String getExecutable() { 'hmmsearch' }

    @Override
    boolean isAvailable() {
        hmmDatabase && new File(hmmDatabase).isFile() &&
            HmmerToolProbe.version(executablePath ?: executable, environment) != null
    }

    @Override
    String getVersion() { HmmerToolProbe.version(executablePath ?: executable, environment) ?: 'unknown' }

    @Override
    String unavailableReason() {
        if (!hmmDatabase) return 'no AntiFam HMM database configured'
        if (!new File(hmmDatabase).isFile()) return "AntiFam HMM database not found: ${hmmDatabase}"
        super.unavailableReason()
    }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (!hmmDatabase) throw new IllegalStateException('antifam: hmmDatabase is unset')
        def outputFile = new File(outputDir, 'antifam_domtbl.tsv')
        def cmd = [
            executablePath ?: executable,
            '--domtblout', outputFile.absolutePath,
            '--noali',
            '--cpu', threads.toString(),
        ]
        if (useCutGa) {
            cmd << '--cut_ga'
        } else {
            cmd << '-E' << '1'
            cmd << '--domE' << domainEvalue.toString()
        }
        cmd << hmmDatabase
        cmd << inputFasta.absolutePath
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // hmmsearch --domtblout always writes '#' header/footer lines, even with no hit.
        def outputFile = requireOutput(outputDir, 'antifam_domtbl.tsv',
            'hmmsearch --domtblout writes comment header/footer on every run')

        Map<String, List<Annotation>> results = [:].withDefault { [] }
        Map<String, Annotation> seenPerProtein = [:]
        Map<String, List<Map>> domainsByKey = [:]

        outputFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split(/\s+/)
            if (fields.length < 22) {
                throw new ToolOutputContractException(
                    "antifam domtblout row has ${fields.length} columns, the format has at least 22: ${line.take(80)}")
            }

            // hmmsearch domtblout format:
            // [0]=target name (protein), [1]=target acc, [2]=tlen
            // [3]=query name (HMM), [4]=query acc (AntiFam ID), [5]=qlen
            // [6..8]=full seq E/score/bias, [11..13]=domain cE/iE/score
            String proteinId = fields[0]
            String antifamAcc  = fields[4]   // ANF00001
            String antifamName = fields[3]   // HMM name
            List<Double> statistics
            try {
                statistics = [6, 7, 12, 13].collect { Double.parseDouble(fields[it]) }
            } catch (NumberFormatException ignored) {
                throw new ToolOutputContractException("antifam: invalid E-value/bit score for ${proteinId}")
            }
            if (statistics.any { !Double.isFinite(it) } || statistics[0] < 0.0d || statistics[2] < 0.0d) {
                throw new ToolOutputContractException("antifam: non-finite statistic or negative E-value for ${proteinId}")
            }
            if (antifamAcc == '-') {
                throw new ToolOutputContractException("antifam: model '${antifamName}' lacks a query accession for ${proteinId}")
            }
            double seqEvalue = statistics[0]
            double seqScore = statistics[1]
            double domEvalue = statistics[2]
            double domScore = statistics[3]

            String acc = antifamAcc.replaceAll(/\.\d+$/, '')

            // Preserve every regional match while retaining one legacy flag per family.
            List<Integer> position
            try {
                position = [2, 5, 9, 10, 15, 16, 17, 18, 19, 20].collect { Integer.parseInt(fields[it]) }
            } catch (NumberFormatException ignored) {
                throw new ToolOutputContractException("antifam: invalid domain coordinates for ${proteinId}")
            }
            def (tlen, qlen, ordinal, total, hmmFrom, hmmTo, aliFrom, aliTo, envFrom, envTo) = position
            if (position.any { it < 1 } || ordinal > total || hmmFrom > hmmTo || hmmTo > qlen ||
                    aliFrom > aliTo || aliTo > tlen || envFrom > aliFrom || envTo < aliTo || envTo > tlen) {
                throw new ToolOutputContractException("antifam: domain coordinates out of bounds for ${proteinId}")
            }
            String key = "${proteinId}|${acc}"
            List<Map> domains = domainsByKey.computeIfAbsent(key) { [] }
            domains << [ordinal: ordinal, total: total, proteinLength: tlen, modelLength: qlen,
                hmmStart: hmmFrom, hmmEnd: hmmTo, alignmentStart: aliFrom, alignmentEnd: aliTo,
                envelopeStart: envFrom, envelopeEnd: envTo, evalue: domEvalue, bitScore: domScore]
            if (seenPerProtein.containsKey(key)) {
                seenPerProtein[key].metadata.domainsJson = JsonOutput.toJson(domains)
                seenPerProtein[key].metadata.domainCount = domains.size()
                return
            }

            // Legacy normalized score: threshold pass (GA), or capped -log10(E).
            // This is NOT a probability of whole-gene inactivity.
            double score = useCutGa ? 1.0 : Math.max(0.0, Math.min(1.0, -Math.log10(Math.max(domEvalue, 1e-300)) / 50.0))

            Annotation annotation = new Annotation(
                type: AnnotationType.PSEUDOGENE,
                value: acc,
                source: name,
                evidenceType: EvidenceType.GENE_MODEL_QUALITY,
                score: score,
                metadata: [
                    name        : antifamName,
                    modelAccession: antifamAcc,
                    qualityIssue: 'spurious_or_local_gene_model_error',
                    gating: useCutGa ? 'gathering_threshold' : 'domain_evalue',
                    domainEvalueThreshold: useCutGa ? null : domainEvalue,
                    scoreSemantics: useCutGa ? 'threshold_pass_indicator' : 'normalized_domain_evalue',
                    coordinateSystem: 'protein_1_based_inclusive',
                    domainsJson: JsonOutput.toJson(domains),
                    domainCount: domains.size(),
                    seqEvalue   : seqEvalue,
                    seqScore    : seqScore,
                    domainEvalue: domEvalue,
                    domainScore : domScore,
                ]
            )
            seenPerProtein[key] = annotation
            results[proteinId] << annotation
        }
        results
    }
}
