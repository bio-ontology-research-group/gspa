package gspa.predictor.neural

import gspa.model.AnnotationType

/**
 * GSPA's <b>DeepGO-PlusPlus-Light</b> predictor (sidecar id
 * {@code deepgo-plusplus-light}): the <b>self-contained, CPU-only</b> sibling of
 * {@link DeepGoPlusPlusPredictor}.
 *
 * <p>Where {@code deepgo-plusplus} is a learned stacker over <em>precomputed</em>
 * per-component score TSVs (it needs the whole upstream component pipeline run
 * first), {@code deepgo-plusplus-light} runs <b>DIAMOND BLAST-KNN</b> and the
 * <b>homology-bridged STRING Net-KNN</b> itself, directly from the query FASTA,
 * then applies the same frozen per-aspect logistic integrator. A single DIAMOND
 * search against the pre-t0 train DB powers both components; the STRING bridge is
 * a dict lookup over a precomputed index (no STRING files read per request), so it
 * predicts functions even for proteins that are <em>not</em> in STRING — the
 * realistic novel-protein case.
 *
 * <p>The default (fast) path needs only the {@code diamond} binary and the Python
 * standard library — no numpy, no torch. Two opt-in, heavier components exist:
 * {@link #cnn} (a CPU 1D-CNN giving a sequence signal for orphan proteins with no
 * homolog; needs PyTorch + a {@code cnn_model.pt}) and {@link #interpro}
 * (InterProScan domain GO terms).
 *
 * <p>Assets come from {@code deepgo-plusplus/service/make_assets.sh} (a bundle of
 * {@code train_db.dmnd}, {@code train_net_index.tsv}, {@code train_terms.tsv},
 * {@code go-dag.tsv}, optional {@code go.obo} + {@code cnn_model.pt}); the frozen
 * integrator JSONs ship in {@code deepgo-plusplus/models/} and are retrainable per
 * UniProt/STRING release via the {@code deepgo-plusplus/} pipeline. The inference
 * core is reused verbatim from {@code deepgo-plusplus/service/predict.py}, so the
 * CLI, the webservice, and the DeepGOWeb embed all share one model implementation.
 */
class DeepGoPlusPlusLightPredictor extends AbstractNeuralSidecarPredictor {

    /**
     * Runtime asset bundle directory ({@code make_assets.sh} output:
     * {@code train_db.dmnd}, {@code train_net_index.tsv}, {@code train_terms.tsv},
     * {@code go-dag.tsv}, optional {@code go.obo} + {@code cnn_model.pt}). Required.
     */
    String assets

    /** Directory of frozen integrator JSONs (default: {@code deepgo-plusplus/models}). */
    String modelsDir

    /** DIAMOND binary. */
    String diamond = 'diamond'

    /** Add the InterProScan domain component (needs {@link #interproscan}). */
    boolean interpro = false

    /** Add the CPU 1D-CNN component (orphan coverage; needs torch + a {@code cnn_model.pt}). */
    boolean cnn = false

    /** InterProScan launcher path (enables the {@link #interpro} component). */
    String interproscan

    /** Override CNN weights path (default {@code <assets>/cnn_model.pt}). */
    String cnnModel

    /** DIAMOND threads. */
    int threads = 8

    /** Homologs voted per query for the BLAST-KNN / Net-KNN components. */
    int topK = 5

    @Override
    String getName() { 'deepgo-plusplus-light' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.GO] as Set }

    @Override
    String getPredictorName() { 'deepgo-plusplus-light' }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!assets) {
            throw new IllegalStateException(
                "${name}: assets is unset. Provide the make_assets.sh bundle dir "
                + "(train_db.dmnd, train_net_index.tsv, train_terms.tsv, go-dag.tsv).")
        }
        List<String> out = [
            '--dgpp-light-assets', new File(assets).absolutePath,
            '--dgpp-light-diamond', diamond,
            '--dgpp-light-threads', threads.toString(),
            '--top-k', topK.toString(),
        ]
        if (modelsDir) {
            out += ['--dgpp-light-models', new File(modelsDir).absolutePath]
        }
        if (interpro) {
            out += ['--dgpp-light-interpro']
        }
        if (cnn) {
            out += ['--dgpp-light-cnn']
        }
        if (interproscan) {
            out += ['--dgpp-light-interproscan', new File(interproscan).absolutePath]
        }
        if (cnnModel) {
            out += ['--dgpp-light-cnn-model', new File(cnnModel).absolutePath]
        }
        out
    }
}
