package gspa.predictor.sites

import gspa.model.AnnotationType
import gspa.predictor.AbstractSiteSidecarPredictor

/**
 * Wraps ScanNet (jertubiana/ScanNet, Apache-2.0) for PPI interface-residue
 * prediction. Requires per-protein PDB / mmCIF structures under
 * {@link #structureDir}; structures can be provisioned by ESMFold or by
 * AFDB lookup (see Nextflow {@code STRUCTURE_PROVIDER} module).
 */
class ScanNetPredictor extends AbstractSiteSidecarPredictor {

    /** Path to the ScanNet repo clone. */
    String modelDir

    /** Directory of {@code <tag>/*.pdb} subdirectories. */
    String structureDir

    @Override
    String getName() { 'scannet' }

    @Override
    String getPredictorName() { 'scannet' }

    @Override
    Set<AnnotationType> getSiteAnnotationTypes() {
        [AnnotationType.PPI_INTERFACE] as Set
    }

    @Override
    boolean isAvailable() {
        super.isAvailable() && modelDir && structureDir &&
            new File(modelDir).exists() && new File(structureDir).exists()
    }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!modelDir)     throw new IllegalStateException("${name}: modelDir required")
        if (!structureDir) throw new IllegalStateException("${name}: structureDir required")
        ['--model-dir', modelDir, '--structure-dir', structureDir]
    }
}
