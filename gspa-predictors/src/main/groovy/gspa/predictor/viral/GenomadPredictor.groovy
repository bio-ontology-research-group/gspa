package gspa.predictor.viral

import gspa.model.AnnotationType
import gspa.predictor.AbstractGenomicRegionSidecarPredictor

/**
 * Wraps geNomad (Camargo et al. 2024, BSD-3-Clause): convolutional
 * neural network + marker-gene HMMs for virus and plasmid detection
 * on bacterial genomes / metagenomes.
 *
 * <p>Emits {@link AnnotationType#PROPHAGE} for integrated proviruses,
 * {@link AnnotationType#VIRAL_CONTIG} for whole-contig viral
 * classifications, and {@link AnnotationType#PLASMID} for plasmid
 * contigs. The actual sub-type per row is driven by the
 * {@code region_type} column emitted by the sidecar.
 */
class GenomadPredictor extends AbstractGenomicRegionSidecarPredictor {

    /** Path to the geNomad model database (download once with
     *  {@code genomad download-database}). */
    String dbPath

    /** Optional Singularity image (biocontainer). */
    String genomadSif

    @Override
    String getName() { 'genomad' }

    @Override
    String getPredictorName() { 'genomad' }

    @Override
    AnnotationType getOutputAnnotationType() { AnnotationType.PROPHAGE }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!dbPath) {
            throw new IllegalStateException("${name}: dbPath required (geNomad model DB)")
        }
        def args = ['--db-path', dbPath]
        if (genomadSif) {
            args += ['--genomad-sif', genomadSif]
        }
        args
    }
}
