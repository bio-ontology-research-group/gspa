package gspa.predictor.viral

import gspa.model.AnnotationType
import gspa.predictor.AbstractGenomicRegionSidecarPredictor

/**
 * Wraps PhiSpy (Akhter et al. 2012; McNair et al. 2018, MIT): random
 * forest on per-CDS features (codon usage, AA composition, transcription
 * strand, gene size) for prophage region detection in bacterial genomes.
 *
 * <p>Input contract: PhiSpy expects a GenBank file (not raw FASTA), so
 * callers must point {@code inputFasta} at a {@code .gbk} file. Emits
 * {@link AnnotationType#PROPHAGE} regions.
 */
class PhiSpyPredictor extends AbstractGenomicRegionSidecarPredictor {

    /** Optional Singularity image (biocontainer). */
    String phispySif

    /** Optional path to a custom trainset (default: genericAll). */
    String trainset

    @Override
    String getName() { 'phispy' }

    @Override
    String getPredictorName() { 'phispy' }

    @Override
    AnnotationType getOutputAnnotationType() { AnnotationType.PROPHAGE }

    @Override
    protected List<String> extraSidecarArgs() {
        def args = []
        if (phispySif) args += ['--phispy-sif', phispySif]
        if (trainset)  args += ['--phispy-trainset', trainset]
        args
    }
}
