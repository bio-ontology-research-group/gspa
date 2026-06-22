package gspa.cli

import gspa.config.ConfigLoader
import gspa.config.GspaConfig
import gspa.metrics.QualityPipeline
import gspa.metrics.QualityReportWriter
import gspa.predictor.AnnotationPipeline
import gspa.predictor.context.OperonPredictor
import picocli.CommandLine.Command
import picocli.CommandLine.Option

@Command(name = 'annotate', description = 'Annotate a genome or set of genomes')
class AnnotateCommand implements Runnable {

    @Option(names = ['-i', '--input'], required = true, description = 'Input FASTA/GenBank file')
    File input

    @Option(names = ['-o', '--output'], description = 'Output directory', defaultValue = 'gspa_output')
    File outputDir

    @Option(names = ['-c', '--config'], description = 'Configuration YAML file')
    File configFile

    @Option(names = ['-k', '--kingdom'], description = 'Organism kingdom: bacteria, archaea, eukaryote, virus, auto')
    String kingdom

    @Option(names = ['--mag'], description = 'Input is a MAG (adjust quality thresholds)')
    boolean mag

    @Option(names = ['--go-owl'], description = 'Path to GO OWL file for quality assessment')
    String goOwl

    @Option(names = ['-t', '--threads'], description = 'Number of threads', defaultValue = '0')
    int threads

    @Option(names = ['--db', '--database'], description = 'Path to sequence similarity database')
    String database

    // --- FoldSeek centroid mode ---

    @Option(names = ['--foldseek-centroid-db'],
            description = 'Path to a FoldSeek function-centroid DB (see benchmark/neural/build_foldseek_centroids.py)')
    String foldseekCentroidDb

    @Option(names = ['--foldseek-centroid-mode'],
            description = 'FoldSeek centroid mode: none | go | ec | both (default from config)')
    String foldseekCentroidMode

    // --- Disorder (Metapredict) ---

    @Option(names = ['--disorder'], description = 'Enable Metapredict disorder predictor')
    boolean disorderEnabled

    // --- Neural sidecar predictors ---

    @Option(names = ['--neural-sidecar'],
            description = 'Absolute path to benchmark/neural/run_neural_predictors.py (required for any --esm2-deepgoplus / --proteinfer / --clean flag)')
    String neuralSidecar

    @Option(names = ['--esm2-deepgoplus'],
            description = 'Enable the ESM2-backed DeepGO-Plus GO predictor (uses base ESM2Head only; no adapter/LP coupling)')
    boolean esm2DeepGoPlusEnabled

    @Option(names = ['--esm2-deepgoplus-checkpoint'],
            description = 'Trained ESM2Head checkpoint (.pt); overrides config')
    String esm2DeepGoPlusCheckpoint

    @Option(names = ['--esm2-deepgoplus-terms'],
            description = 'GO vocabulary file (one term per line); overrides config')
    String esm2DeepGoPlusTerms

    @Option(names = ['--esm2-deepgoplus-model'],
            description = 'ESM2 variant matching the checkpoint (e.g. esm2_t12_35M_UR50D)')
    String esm2DeepGoPlusModel

    @Option(names = ['--proteinfer'], description = 'Enable ProteInfer GO predictor')
    boolean proteinferEnabled

    @Option(names = ['--proteinfer-model-dir'],
            description = 'ProteInfer model directory; overrides config')
    String proteinferModelDir

    @Option(names = ['--clean'], description = 'Enable CLEAN EC predictor')
    boolean cleanEnabled

    @Option(names = ['--clean-model-dir'],
            description = 'CLEAN checkpoint directory; overrides config')
    String cleanModelDir

    @Option(names = ['--deepgo-plusplus', '--cafa-baseline'],
            description = 'Enable the DeepGO-PlusPlus learned stacker over precomputed component scores '
                    + '(legacy alias: --cafa-baseline)')
    boolean deepGoPlusPlusEnabled

    @Option(names = ['--deepgo-plusplus-integrator', '--cafa-baseline-integrator'],
            description = 'Frozen integrator JSON (deepgo-plusplus/pipeline/train_integrator.py --save-model); overrides config')
    String deepGoPlusPlusIntegrator

    @Option(names = ['--deepgo-plusplus-components-dir', '--cafa-baseline-components-dir'],
            description = 'Directory of per-component score TSVs (<component>.tsv[.gz]); overrides config')
    String deepGoPlusPlusComponentsDir

    @Option(names = ['--deepgo-plusplus-dag', '--cafa-baseline-dag'],
            description = 'GO DAG file (child\\tancestor) for true-path propagation; overrides config')
    String deepGoPlusPlusDag

    @Override
    void run() {
        println "GSPA annotate: ${input}"
        println "Output: ${outputDir}"
        println "Kingdom: ${kingdom ?: 'auto'}"
        println "MAG mode: ${mag}"
        println ""

        def overrides = [:]
        if (kingdom) overrides['kingdom'] = kingdom
        if (goOwl) overrides['go-owl'] = goOwl
        if (database) overrides['database'] = database
        def config = ConfigLoader.buildConfig(configFile, kingdom, overrides)
        if (mag) config.input.type = 'mag'

        applyPredictorFlags(config)

        // Build and configure annotation pipeline
        def pipeline = new AnnotationPipeline(config)
        pipeline.configure()

        // Run annotation
        def genome = pipeline.annotate(input)

        // Write output
        pipeline.writeOutput(genome, outputDir)

        // Always persist operon detection alongside the annotations so
        // downstream `gspa integrate --operons` and `gspa visualize` find it
        // without manual prep. Cheap to re-detect (O(n log n) over CDS).
        try {
            writeOperonFiles(genome, outputDir, config)
        } catch (Exception e) {
            System.err.println "  [warn] Failed to write operon TSVs: ${e.message}"
        }

        // Run quality evaluation if GO ontology provided
        if (goOwl && genome.allGoTerms().size() > 0) {
            println ""
            println "Running quality evaluation..."
            def qualityPipeline = new QualityPipeline(config)
                .goOwlFile(goOwl)
                .essentialFunctionsForDomain(config.resolveOrganismDomain())
                .initializeLite()

            def report = qualityPipeline.evaluate(genome)
            def reportFile = new File(outputDir, "${genome.id}_quality.json")
            QualityReportWriter.writeJson(report, reportFile)
            println "Quality report: ${reportFile}"

            qualityPipeline.dispose()
        }

        println ""
        println "Done. Output in: ${outputDir}"
    }

    /**
     * Detect operons on the annotated genome and persist three TSVs alongside
     * the other annotation outputs:
     *
     *   operons.tsv               human / viz-friendly: id, contig, start, end,
     *                             strand, n_members, members
     *   protein_to_operon.tsv     reverse index: protein_id, operon_id, position,
     *                             operon_size
     *   operons_for_integrate.tsv flat format that {@code gspa integrate
     *                             --operons} expects (one operon per line, tab-
     *                             separated protein IDs)
     *
     * Operon parameters come from the same config block that the inline
     * OperonPredictor uses inside AnnotationPipeline (so the persisted file
     * always matches whatever predictor calls happened during annotate).
     */
    private void writeOperonFiles(gspa.model.Genome genome, File outputDir, GspaConfig config) {
        if (genome == null || genome.proteins == null || genome.proteins.isEmpty()) return
        // Skip if no contig coordinates — operon detection is positional.
        if (!genome.proteins.any { it.start != null && it.end != null }) return
        def detector = new OperonPredictor(
            maxIntergenicDistance: config.predictors.operons.maxIntergenicDistance,
            requireSameStrand:     config.predictors.operons.requireSameStrand,
        )
        def operons = detector.detectOperons(genome)
        if (operons.isEmpty()) return
        outputDir.mkdirs()
        // 1. operons.tsv (verbose)
        new File(outputDir, "operons.tsv").withWriter { w ->
            w.writeLine(['operon_id','contig','start','end','strand','n_members','members'].join('\t'))
            operons.eachWithIndex { op, i ->
                String opId = String.format(Locale.ROOT, 'op_%05d', i + 1)
                w.writeLine([
                    opId, op.contigId, op.start.toString(), op.end.toString(), op.strand.symbol,
                    op.size.toString(), op.genes*.id.join(','),
                ].join('\t'))
            }
        }
        // 2. protein_to_operon.tsv
        new File(outputDir, "protein_to_operon.tsv").withWriter { w ->
            w.writeLine(['protein_id','operon_id','position','operon_size'].join('\t'))
            operons.eachWithIndex { op, i ->
                String opId = String.format(Locale.ROOT, 'op_%05d', i + 1)
                op.genes.eachWithIndex { gene, j ->
                    w.writeLine([gene.id, opId, (j + 1).toString(), op.size.toString()].join('\t'))
                }
            }
        }
        // 3. operons_for_integrate.tsv (the format `gspa integrate --operons` expects)
        new File(outputDir, "operons_for_integrate.tsv").withWriter { w ->
            operons.each { op ->
                w.writeLine(op.genes*.id.join('\t'))
            }
        }
        println "Operons: ${operons.size()} (written to operons.tsv, protein_to_operon.tsv, operons_for_integrate.tsv)"
    }

    private void applyPredictorFlags(GspaConfig config) {
        // FoldSeek centroid overrides (enable structure predictor if they're used)
        if (foldseekCentroidDb) {
            config.predictors.structure.enabled = true
            config.predictors.structure.centroidDb = foldseekCentroidDb
        }
        if (foldseekCentroidMode) {
            config.predictors.structure.enabled = true
            config.predictors.structure.centroidMode = foldseekCentroidMode
        }

        // Disorder
        if (disorderEnabled) {
            config.predictors.disorder.enabled = true
        }

        // Neural sidecar (path is shared across all three neural predictors)
        if (neuralSidecar) {
            config.predictors.neural.sidecarScript = neuralSidecar
        }

        if (esm2DeepGoPlusEnabled) {
            config.predictors.neural.esm2DeepGoPlus.enabled = true
        }
        if (esm2DeepGoPlusCheckpoint) {
            config.predictors.neural.esm2DeepGoPlus.checkpoint = esm2DeepGoPlusCheckpoint
        }
        if (esm2DeepGoPlusTerms) {
            config.predictors.neural.esm2DeepGoPlus.terms = esm2DeepGoPlusTerms
        }
        if (esm2DeepGoPlusModel) {
            config.predictors.neural.esm2DeepGoPlus.model = esm2DeepGoPlusModel
        }

        if (proteinferEnabled) {
            config.predictors.neural.proteinfer.enabled = true
        }
        if (proteinferModelDir) {
            config.predictors.neural.proteinfer.modelDir = proteinferModelDir
        }

        if (cleanEnabled) {
            config.predictors.neural.clean.enabled = true
        }
        if (cleanModelDir) {
            config.predictors.neural.clean.modelDir = cleanModelDir
        }

        if (deepGoPlusPlusEnabled) {
            config.predictors.neural.deepGoPlusPlus.enabled = true
        }
        if (deepGoPlusPlusIntegrator) {
            config.predictors.neural.deepGoPlusPlus.integrator = deepGoPlusPlusIntegrator
        }
        if (deepGoPlusPlusComponentsDir) {
            config.predictors.neural.deepGoPlusPlus.componentsDir = deepGoPlusPlusComponentsDir
        }
        if (deepGoPlusPlusDag) {
            config.predictors.neural.deepGoPlusPlus.dag = deepGoPlusPlusDag
        }
    }
}
