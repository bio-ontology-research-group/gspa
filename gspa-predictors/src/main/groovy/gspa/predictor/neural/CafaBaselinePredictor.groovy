package gspa.predictor.neural

import gspa.model.AnnotationType

/**
 * GSPA's CAFA6 {@code cafa-baseline} predictor: a <b>learned stacker</b> over
 * the per-component GO scores already produced by the upstream GSPA
 * predictors (DIAMOND/BLAST-KNN, FoldSeek-KNN, CLEAN, InterPro, an ESM2 MLP
 * head and a structure-aware ProstT5 head).
 *
 * <p>It does <em>not</em> re-run those tools. Instead it consumes a directory
 * of precomputed per-component score TSVs ({@code <component>.tsv[.gz]} with
 * rows {@code protein\\tterm\\tscore}) and applies a <b>frozen</b> per-aspect
 * logistic-regression integrator (trained once via
 * {@code benchmark/neural/train_ltr_integrator.py --save-model} and saved as
 * JSON). This replaces the naive max-merge that, on our original CAFA6 entry,
 * destroyed the molecular-function signal; the learned integration recovered
 * novel-protein IA-weighted f_w from 0.359 to 0.483 on the CAFA6
 * reconstruction (vs the 0.524 winner).
 *
 * <p>The sidecar runner ({@code --predictor cafa-baseline}) propagates each
 * component to its GO ancestors (max), forms per-(protein, term) candidates
 * per aspect, and emits {@code sigmoid(w·x + b)} as the combined GO score.
 */
class CafaBaselinePredictor extends AbstractNeuralSidecarPredictor {

    /** Frozen integrator model JSON (train_ltr_integrator.py --save-model). Required. */
    String integrator

    /**
     * Directory of per-component score TSVs ({@code <component>.tsv[.gz]},
     * rows {@code protein\\tterm\\tscore}). The component file names must
     * match the {@code components} list in the integrator JSON. Required.
     */
    String componentsDir

    /** GO DAG file ({@code child\\tancestor}) for true-path propagation. Required. */
    String dag

    @Override
    String getName() { 'cafa-baseline' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.GO] as Set }

    @Override
    String getPredictorName() { 'cafa-baseline' }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!integrator) {
            throw new IllegalStateException(
                "${name}: integrator is unset. Provide the frozen integrator JSON path.")
        }
        if (!componentsDir) {
            throw new IllegalStateException(
                "${name}: componentsDir is unset. Provide the per-component scores directory.")
        }
        if (!dag) {
            throw new IllegalStateException(
                "${name}: dag is unset. Provide the go-dag.tsv (child\\tancestor) path.")
        }
        [
            '--integrator', new File(integrator).absolutePath,
            '--components-dir', new File(componentsDir).absolutePath,
            '--dag', new File(dag).absolutePath,
        ]
    }
}
