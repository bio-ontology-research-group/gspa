package gspa.predictor.viral

import gspa.model.AnnotationType
import gspa.predictor.AbstractGenomicRegionSidecarPredictor

/**
 * Wraps CheckV (Nayfach et al. 2021, BSD-3-Clause): viral contig
 * completeness + contamination assessment via marker-gene HMMs.
 *
 * <p>Emits {@link AnnotationType#VIRAL_CONTIG} per contig. The score
 * field carries the completeness percentage (0–1 normalised); the
 * {@code attributes} column carries CheckV's quality and contamination
 * estimates as {@code key=val|key=val}.
 */
class CheckVPredictor extends AbstractGenomicRegionSidecarPredictor {

    String dbPath
    String checkvSif
    int threads = 4

    @Override
    String getName() { 'checkv' }

    @Override
    String getPredictorName() { 'checkv' }

    @Override
    AnnotationType getOutputAnnotationType() { AnnotationType.VIRAL_CONTIG }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!dbPath) {
            throw new IllegalStateException("${name}: dbPath required (CheckV DB)")
        }
        def args = ['--db-path', dbPath, '--threads', threads.toString()]
        if (checkvSif) {
            args += ['--checkv-sif', checkvSif]
        }
        args
    }
}
