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

    /**
     * CheckV assesses contigs already believed to be viral. The sidecar reports one WHOLE-contig
     * row (1..contig_length, region_type viral_contig) per quality_summary.tsv line; CheckV's
     * trimmed provirus coordinates are not emitted. The score is completeness/100 (not a
     * probability that the contig is viral) and rows below minScore are dropped by the sidecar,
     * so a bacterial genome yields zero rows regardless of its prophage content.
     */
    @Override
    protected Map<String, Object> regionMetadata() {
        [scoreSemantics: 'checkv completeness/100', coordinateScope: 'whole_contig']
    }

    @Override
    protected String successNote() {
        'CheckV scores contig completeness, not viral origin: zero rows means no contig reached ' +
            'minScore completeness, not that no viral sequence is present'
    }

    @Override
    protected String configurationProblem() {
        dbPath ? null : 'dbPath (CheckV database) is not configured'
    }

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
