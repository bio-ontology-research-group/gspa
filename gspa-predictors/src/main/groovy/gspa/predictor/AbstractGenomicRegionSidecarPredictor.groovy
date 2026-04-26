package gspa.predictor

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType

/**
 * Base class for predictors that delegate to
 * {@code benchmark/neural/run_genomic_predictors.py}.
 *
 * <p>Output is the 6-column genomic-region TSV:
 * {@code contig_id<TAB>region_start<TAB>region_end<TAB>region_type<TAB>score<TAB>attributes}
 * with 1-based inclusive contig coordinates.
 *
 * <p>Each row becomes one {@link Annotation} with
 * {@link Annotation#contigId}, {@link Annotation#genomicStart},
 * {@link Annotation#genomicEnd}, and
 * {@link EvidenceType#GENOMIC_REGION_ML} (correlation group {@code viral}).
 *
 * <p>Subclasses provide:
 * <ul>
 *   <li>{@link #getPredictorName} — the {@code --predictor} flag value
 *   <li>{@link #getOutputAnnotationType} — JVM annotation type
 *       (typically {@code PROPHAGE}, {@code PLASMID}, or
 *       {@code VIRAL_CONTIG})
 *   <li>{@link #extraSidecarArgs} — per-tool flags
 * </ul>
 *
 * <p>Manifest contract differs from the protein sidecars: it carries a
 * {@code genome_fasta} column instead of {@code fasta_path} and an
 * optional {@code gff_path} for tools that need per-CDS features
 * (PhiSpy). The {@link #buildCommand} method writes the manifest
 * accordingly.
 */
abstract class AbstractGenomicRegionSidecarPredictor extends AbstractToolPredictor {

    /** Absolute path to {@code run_genomic_predictors.py}. */
    String sidecarScript

    String pythonExecutable = 'python3'

    /** Drop predictions below this score. */
    double minScore = 0.5

    /** Synthetic tag used in the manifest; identifies the output filename. */
    String tag = 'query'

    /** Optional GFF path for tools that need per-CDS features (PhiSpy). */
    String gffPath

    /** Evidence type reported on emitted annotations. */
    EvidenceType evidenceType = EvidenceType.GENOMIC_REGION_ML

    /** The {@code --predictor} flag value for the sidecar. */
    abstract String getPredictorName()

    /** AnnotationType this predictor primarily emits (e.g. PROPHAGE). */
    abstract AnnotationType getOutputAnnotationType()

    /** Extra sidecar CLI flags (per-tool). */
    protected List<String> extraSidecarArgs() { [] }

    @Override
    String getExecutable() { pythonExecutable }

    @Override
    boolean isAvailable() {
        if (!sidecarScript || !(new File(sidecarScript).exists())) return false
        try {
            def proc = [pythonExecutable, '--version'].execute()
            proc.waitForOrKill(10000)
            return proc.exitValue() == 0
        } catch (Exception ignored) {
            return false
        }
    }

    @Override
    Set<AnnotationType> getOutputTypes() { [outputAnnotationType] as Set }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (!sidecarScript) {
            throw new IllegalStateException(
                "${name}: sidecarScript is unset. Point it at run_genomic_predictors.py.")
        }
        def manifest = new File(outputDir, 'manifest.tsv')
        def gff = gffPath ?: '-'
        manifest.text = "tag\tgenome_fasta\tgff_path\toutput_dir\n" +
                "${tag}\t${inputFasta.absolutePath}\t${gff}\t${outputDir.absolutePath}\n"
        def cmd = [
            pythonExecutable,
            new File(sidecarScript).absolutePath,
            '--predictor', predictorName,
            '--manifest', manifest.absolutePath,
            '--min-score', minScore.toString(),
        ] as List<String>
        cmd.addAll(extraSidecarArgs())
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def tsv = new File(outputDir, "${tag}.${predictorName}.genomic.tsv")
        if (!tsv.exists()) return [:]
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        boolean firstLine = true
        tsv.eachLine { line ->
            if (firstLine) { firstLine = false; return }
            if (!line || line.startsWith('#')) return
            def parts = line.split('\t')
            if (parts.size() < 5) return
            try {
                String contig = parts[0]
                int start = parts[1] as int
                int end = parts[2] as int
                String regionType = parts[3]
                double score = parts[4] as double
                String attributesRaw = parts.size() >= 6 ? parts[5] : ''
                if (Double.isNaN(score) || score < minScore) return
                // Resolve annotation type from region_type when the predictor
                // emits multiple kinds (e.g. geNomad: prophage / plasmid /
                // viral_contig); fall back to the predictor's primary type.
                AnnotationType t = outputAnnotationType
                switch (regionType) {
                    case 'prophage':
                        t = AnnotationType.PROPHAGE
                        break
                    case 'plasmid':
                        t = AnnotationType.PLASMID
                        break
                    case 'viral_contig':
                        t = AnnotationType.VIRAL_CONTIG
                        break
                }
                Map<String, Object> meta = [:]
                if (attributesRaw) {
                    attributesRaw.tokenize('|').each { kv ->
                        def i = kv.indexOf('=')
                        if (i > 0) meta[kv.substring(0, i)] = kv.substring(i + 1)
                    }
                }
                // Genomic-region annotations are keyed by contig (not protein);
                // we use the contig id as the map key so downstream consumers
                // can group / flatten by chromosomal location.
                results[contig] << new Annotation(
                    type: t,
                    value: regionType,
                    source: name,
                    score: score,
                    evidence: 'IEA',
                    evidenceType: evidenceType,
                    contigId: contig,
                    genomicStart: start,
                    genomicEnd: end,
                    regionType: regionType,
                    metadata: meta,
                )
            } catch (NumberFormatException ignored) {
                // skip malformed row
            }
        }
        results
    }
}
