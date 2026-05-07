package gspa.predictor.structure

import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import gspa.predictor.neural.AbstractNeuralSidecarPredictor

/**
 * Wraps metagenomic-deepFRI (Bezshapkin et al. 2026, BSD-3-Clause) for
 * GO and EC prediction. Successor to the upstream {@link DeepFriPredictor}
 * with FoldComp-database structure retrieval + ONNX inference (2-12×
 * speedup over the original TensorFlow DeepFRI).
 *
 * <p>Two operating modes:
 * <ul>
 *   <li>sequence-only ({@code skipPdb = true}): runs the DeepFRI CNN
 *       branch only — no FoldComp DB needed. Fastest; matches the
 *       v1.5.0 mdF benchmark line in {@code benchmark/RESULTS.md}.</li>
 *   <li>structure-supplied ({@code foldcompDb} set): retrieves homolog
 *       structures via MMseqs2 against a FoldComp DB (AFDBv4 / ESM
 *       Atlas / PDB100), transfers contact map to query via PyOpal
 *       alignment, runs the DeepFRI GCN branch. Heaviest; needs the
 *       FoldComp DB on disk (~5 GB AFDBv4, ~14 GB ESM Atlas).</li>
 * </ul>
 *
 * <p>Emits {@link EvidenceType#STRUCTURE_DEEPLEARNING} on every
 * annotation so it correlation-groups with {@code DeepFriPredictor},
 * {@code FoldSeekPredictor}, and any future structure-side predictor —
 * the integrator's {@code structure} group keeps only the strongest
 * claim per (protein, term) before Noisy-OR.
 */
class MdFPredictor extends AbstractNeuralSidecarPredictor {

    /** Path to mDeepFRI weights dir (from {@code mDeepFRI get-models --version 1.0}). */
    String weightsDir

    /** Path to the {@code mDeepFRI} binary. {@code null} = first on PATH. */
    String mdfExecutable

    /** Sequence-only path; skip PDB100 + FoldComp DB lookup. */
    boolean skipPdb = true

    /** Optional FoldComp DB path (structure-supplied path). */
    String foldcompDb

    /** Threads forwarded to mDeepFRI. */
    int threads = 4

    MdFPredictor() {
        evidenceType = EvidenceType.STRUCTURE_DEEPLEARNING
        // mDeepFRI's own threshold defaults to 0.1; matches the sidecar.
        minScore = 0.1
    }

    @Override
    String getName() { 'mdf' }

    @Override
    String getPredictorName() { 'mdf' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.GO, AnnotationType.EC] as Set }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!weightsDir) {
            throw new IllegalStateException(
                "${name}: weightsDir must point at the mDeepFRI weights directory")
        }
        def args = ['--mdf-weights', weightsDir,
                    '--mdf-threads', threads.toString()] as List<String>
        if (skipPdb) {
            args << '--mdf-skip-pdb'
        }
        if (foldcompDb) {
            args.addAll(['--mdf-foldcomp-db', foldcompDb])
        }
        if (mdfExecutable) {
            args.addAll(['--mdf-executable', mdfExecutable])
        }
        args
    }
}
