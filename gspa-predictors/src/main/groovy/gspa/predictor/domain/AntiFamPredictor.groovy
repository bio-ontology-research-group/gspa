package gspa.predictor.domain

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * AntiFam HMM scanner for spurious / pseudogenic ORF detection.
 *
 * AntiFam (EBI/Pfam) is a curated HMM library whose profiles represent
 * <i>known-spurious</i> protein families — ORFs that look real to gene
 * callers (Prodigal, Pyrodigal, Bakta, Prokka) but are pseudogenes,
 * misannotated coding regions, or sequencing artefacts. A hit against
 * any AntiFam HMM is strong evidence that the ORF should not be treated
 * as a functional gene.
 *
 * <p>The Bakta annotator uses this exact step as its pseudogene filter.
 * This predictor surfaces the same signal inside GSPA: it emits a
 * {@link AnnotationType#PSEUDOGENE} annotation per hit, leaving
 * downstream consumers (integrator, quality scorer, visualiser) free to
 * drop or down-weight claims attached to flagged proteins.</p>
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
    List<String> buildCommand(File inputFasta, File outputDir) {
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
        def outputFile = new File(outputDir, 'antifam_domtbl.tsv')
        if (!outputFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }
        Set<String> seenPerProtein = [] as Set

        outputFile.eachLine { line ->
            if (line.startsWith('#')) return
            def fields = line.split(/\s+/)
            if (fields.length < 22) return

            // hmmsearch domtblout format:
            // [0]=target name (protein), [1]=target acc, [2]=tlen
            // [3]=query name (HMM), [4]=query acc (AntiFam ID), [5]=qlen
            // [6..8]=full seq E/score/bias, [11..13]=domain cE/iE/score
            String proteinId = fields[0]
            String antifamAcc  = fields[4]   // ANF00001
            String antifamName = fields[3]   // HMM name
            double seqEvalue = fields[6] as double
            double seqScore  = fields[7] as double
            double domEvalue = fields[12] as double
            double domScore  = fields[13] as double

            String acc = antifamAcc.replaceAll(/\.\d+$/, '')

            // De-duplicate per (protein, AntiFam family) — a single pseudo flag
            // per family is enough, even when an HMM hits multiple domains.
            String key = "${proteinId}|${acc}"
            if (!seenPerProtein.add(key)) return

            // Score: AntiFam treats any --cut_ga pass as a confident pseudo
            // call, so we always emit score=1.0. When users switch to
            // E-value gating, scale by -log10(E) capped at 1.0.
            double score = useCutGa ? 1.0 : Math.min(1.0, -Math.log10(Math.max(domEvalue, 1e-300)) / 50.0)

            results[proteinId] << new Annotation(
                type: AnnotationType.PSEUDOGENE,
                value: acc,
                source: name,
                score: score,
                metadata: [
                    name        : antifamName,
                    seqEvalue   : seqEvalue,
                    seqScore    : seqScore,
                    domainEvalue: domEvalue,
                    domainScore : domScore,
                ]
            )
        }
        results
    }
}
