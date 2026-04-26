package gspa.predictor

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType

/**
 * Base class for predictors that delegate to
 * {@code benchmark/neural/run_site_predictors.py}.
 *
 * <p>Output format is the fixed 5-column site TSV:
 * {@code protein_id<TAB>position<TAB>site_type<TAB>score<TAB>annotation_type}.
 *
 * <p>Each row becomes one {@link Annotation} representing a single residue
 * (a region with {@code regionStart == regionEnd}). The
 * {@code annotation_type} column drives {@link Annotation#type} —
 * {@code PTM_SITE}, {@code PPI_INTERFACE}, etc.
 */
abstract class AbstractSiteSidecarPredictor extends AbstractToolPredictor {

    /** Absolute path to {@code run_site_predictors.py}. */
    String sidecarScript

    String pythonExecutable = 'python3'
    double minScore = 0.5
    String tag = 'query'

    EvidenceType evidenceType = EvidenceType.SEQUENCE_REGION_ML

    abstract String getPredictorName()

    /** AnnotationTypes the predictor may emit (used for getOutputTypes). */
    abstract Set<AnnotationType> getSiteAnnotationTypes()

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
    Set<AnnotationType> getOutputTypes() { siteAnnotationTypes }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (!sidecarScript) {
            throw new IllegalStateException(
                "${name}: sidecarScript is unset. Point it at run_site_predictors.py.")
        }
        def manifest = new File(outputDir, 'manifest.tsv')
        manifest.text = "tag\tfasta_path\toutput_dir\n" +
                "${tag}\t${inputFasta.absolutePath}\t${outputDir.absolutePath}\n"
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
        def tsv = new File(outputDir, "${tag}.${predictorName}.tsv")
        if (!tsv.exists()) return [:]
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        boolean firstLine = true
        tsv.eachLine { line ->
            if (firstLine) { firstLine = false; return }
            if (!line || line.startsWith('#')) return
            def parts = line.split('\t')
            if (parts.size() < 5) return
            AnnotationType t
            try {
                t = AnnotationType.valueOf(parts[4].toUpperCase(Locale.ROOT))
            } catch (IllegalArgumentException ignored) {
                return
            }
            try {
                String pid = parts[0]
                int pos = parts[1] as int
                String siteType = parts[2]
                double score = parts[3] as double
                if (Double.isNaN(score) || score < minScore) return
                results[pid] << new Annotation(
                    type: t,
                    value: siteType,
                    source: name,
                    score: score,
                    evidence: 'IEA',
                    evidenceType: evidenceType,
                    regionStart: pos,
                    regionEnd: pos,
                    regionType: siteType,
                )
            } catch (NumberFormatException ignored) {
                // Skip malformed row
            }
        }
        results
    }
}
