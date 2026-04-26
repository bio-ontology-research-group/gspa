package gspa.cli

import gspa.config.ConfigLoader
import gspa.config.GspaConfig
import gspa.metrics.QualityPipeline
import gspa.metrics.QualityReportWriter
import gspa.predictor.AnnotationPipeline
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
    }
}
