package gspa.cli

import gspa.config.ConfigLoader
import gspa.config.GspaConfig
import gspa.metrics.QualityPipeline
import gspa.metrics.QualityReportWriter
import gspa.predictor.AnnotationPipeline
import gspa.predictor.context.OperonPredictor
import picocli.CommandLine.Command
import picocli.CommandLine.Option

@Command(name = 'annotate', description = 'Annotate a genome or set of genomes', mixinStandardHelpOptions = true)
class AnnotateCommand implements Runnable {

    @Option(names = ['-i', '--input'],
            description = 'Primary input: nucleotide FASTA (genes called), or a GFF3 '
                    + '(optionally with an embedded ##FASTA block). Alias of --genome for FASTA input.')
    File input

    @Option(names = ['--genome'],
            description = 'Genome / metagenome nucleotide FASTA (one contig per sequence). '
                    + 'Pair with --gff3 to translate supplied CDS instead of calling genes.')
    File genomeFasta

    @Option(names = ['--gff3', '--gff'],
            description = 'GFF3 annotation paired with --genome (or --proteins). CDS coordinates '
                    + 'are reused: features are translated against the genome sequences and mapped '
                    + 'to their contig, rather than re-calling genes.')
    File gff3

    @Option(names = ['--proteins'],
            description = 'Pre-called protein FASTA. With --gff3, proteins are joined to contigs '
                    + 'via the CDS attributes (protein_id/ID/locus_tag); alone, they are annotated '
                    + 'on a single synthetic contig.')
    File proteins

    @Option(names = ['-o', '--output'], description = 'Output directory', defaultValue = 'gspa_output')
    File outputDir

    @Option(names = ['-c', '--config'], description = 'Configuration YAML file')
    File configFile

    @Option(names = ['-k', '--kingdom'], description = 'Organism kingdom: bacteria, archaea, eukaryote, virus, auto')
    String kingdom

    @Option(names = ['--mag'], description = 'Input is a MAG (adjust quality thresholds)')
    boolean mag

    @Option(names = ['--go-owl'], description = 'Path to GO OWL file for quality assessment (enables genome-scale metrics)')
    String goOwl

    @Option(names = ['--no-metrics'], description = 'Skip genome-scale metrics even when --go-owl is given')
    boolean noMetrics

    @Option(names = ['--metrics-scope'],
            description = 'Granularity of genome-scale metrics: contig (default, per contig), genome (pooled), or both')
    String metricsScope

    @Option(names = ['--enforce-consistency'],
            description = 'Run the SAT taxon-constraint enforcement post-pass on predicted GO annotations')
    boolean enforceConsistency

    @Option(names = ['--enforce-consistency-mode'],
            description = 'Action for inconsistent annotations: remove (default), downrank, flag, or '
                    + 'minimal-flip (joint weighted-MaxSAT min-cost demotion, resolves co-annotation conflicts)')
    String enforceConsistencyMode

    @Option(names = ['--enforce-completeness'],
            description = 'Promote each missing essential function onto the best-evidenced candidate protein '
                    + '(imputed annotation, evidence ISC). Provenance records the basis.')
    boolean enforceCompleteness

    @Option(names = ['--enforce-coherence'],
            description = 'Fix obligate heteromeric-complex singletons (demote / promote partner) and promote '
                    + 'missing has_part partners. Needs ELK (full init) and, for complexes, --complex-terms.')
    boolean enforceCoherence

    @Option(names = ['--complex-terms'],
            description = 'TSV of complex GO terms (term<TAB>h|n|a) for complex-coherence enforcement')
    String complexTermsFile

    @Option(names = ['--no-provenance'],
            description = 'Do not record enforcement provenance (omit the provenance column + enforcement_actions.tsv)')
    boolean noProvenance

    @Option(names = ['--taxon-constraints'],
            description = 'Explicit taxon-constraint source for consistency: go-computed-taxon-constraints.obo '
                    + '(.obo) or a TSV (GO-term<TAB>only_in|never_in<TAB>taxon). Needed when --go-owl is go-basic.')
    String taxonConstraintsFile

    @Option(names = ['--taxonomy'],
            description = 'Taxon hierarchy + disjointness file (Term<TAB>is_a|disjoint_from|union_of<TAB>term); '
                    + 'defaults to the bundled NCBI disjointness backbone')
    String taxonomyFile

    @Option(names = ['--taxon'],
            description = 'Assert the organism NCBI taxon during consistency enforcement so terms that cannot occur '
                    + 'in its lineage are flagged (e.g. eukaryote-only terms on a bacterium). Accepts NCBITaxon_2, 2, '
                    + 'or a kingdom name (bacteria/archaea/eukaryote/virus). For single-organism inputs.')
    String organismTaxon

    @Option(names = ['--infer-taxon'],
            description = 'Infer the organism domain taxon from the predicted functions via the GO taxon constraints '
                    + '(Asaad-style) and use it for consistency enforcement, instead of assuming one. Writes '
                    + '<genome>_taxon_inference.tsv with the per-candidate evidence. Ignored if --taxon is given.')
    boolean inferTaxon

    @Option(names = ['--infer-taxon-min-score'],
            description = 'Only predictions scoring at least this (0..1) inform taxon inference; the low-score tail '
                    + 'is cross-domain noise that biases the call. Default ${DEFAULT-VALUE}.')
    Double inferTaxonMinScore

    @Option(names = ['--min-score'],
            description = 'Drop GO annotations scoring below this (0..1) right after prediction, before metrics / '
                    + 'output. Prunes the low-confidence GO-DAG tail (a learned predictor propagates ~thousands of '
                    + 'near-zero terms per protein up the hierarchy), keeping genome-scale runs and outputs '
                    + 'tractable. 0 = keep everything (default).')
    Double minScore

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

    @Option(names = ['--base-predictor'],
            description = 'Use DeepGO-PlusPlus as the base function predictor: '
                    + 'full (GPU, precomputed components) or light (CPU, self-contained from FASTA). '
                    + 'Mutually exclusive. Provide the matching asset flags '
                    + '(full: --deepgo-plusplus-integrator/-components-dir/-dag; light: --deepgo-plusplus-light-assets).')
    String basePredictor

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

    @Option(names = ['--deepgo-plusplus-light'],
            description = 'Enable DeepGO-PlusPlus-Light: self-contained CPU GO predictor '
                    + '(runs DIAMOND BLAST-KNN + homology-bridged STRING Net-KNN itself)')
    boolean deepGoPlusPlusLightEnabled

    @Option(names = ['--deepgo-plusplus-light-assets'],
            description = 'make_assets.sh bundle dir (train_db.dmnd, train_net_index.tsv, '
                    + 'train_terms.tsv, go-dag.tsv [, go.obo, cnn_model.pt]); overrides config')
    String deepGoPlusPlusLightAssets

    @Option(names = ['--deepgo-plusplus-light-models'],
            description = 'Dir of frozen integrator JSONs (default: deepgo-plusplus/models); overrides config')
    String deepGoPlusPlusLightModels

    @Option(names = ['--deepgo-plusplus-light-interpro'],
            description = 'Add the InterProScan domain component (needs --deepgo-plusplus-light-interproscan)')
    boolean deepGoPlusPlusLightInterpro

    @Option(names = ['--deepgo-plusplus-light-cnn'],
            description = 'Add the CPU 1D-CNN component (orphan coverage; needs torch + cnn_model.pt)')
    boolean deepGoPlusPlusLightCnn

    @Option(names = ['--deepgo-plusplus-light-interproscan'],
            description = 'InterProScan launcher path (enables the interpro component); overrides config')
    String deepGoPlusPlusLightInterproscan

    @Override
    void run() {
        // Resolve the primary input across the supported forms:
        //   --genome g.fna [--gff3 a.gff3] [--proteins p.faa]
        //   -i annotated.gff3            (GFF3, possibly with ##FASTA)
        //   --proteins p.faa --gff3 a.gff3
        //   --proteins p.faa             (protein-only)
        //   -i g.fna                     (nucleotide FASTA; genes called)
        File primary = input ?: genomeFasta
        boolean proteinsOnly = false
        if (!primary) {
            if (gff3) {
                primary = gff3                 // GFF3 drives; proteins (if any) joined to contigs
            } else if (proteins) {
                primary = proteins             // bare protein FASTA
                proteinsOnly = true
            }
        }
        if (!primary) {
            throw new IllegalArgumentException(
                "No input given. Use --genome/-i (FASTA), --gff3 (GFF3), or --proteins (protein FASTA).")
        }

        println "GSPA annotate: ${primary}"
        if (gff3) println "GFF3: ${gff3}"
        if (proteins) println "Proteins: ${proteins}"
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
        if (gff3) config.input.gff3 = gff3.absolutePath
        if (proteins) config.input.proteinFasta = proteins.absolutePath
        config.input.proteinsOnly = proteinsOnly

        // Genome-scale metrics + SAT enforcement toggles
        if (noMetrics) config.quality.enabled = false
        if (metricsScope) config.quality.scope = metricsScope
        if (enforceConsistency) config.quality.consistency.enforce = true
        if (enforceConsistencyMode) config.quality.consistency.enforceMode = enforceConsistencyMode
        if (taxonConstraintsFile) config.quality.consistency.constraintsFile = taxonConstraintsFile
        if (taxonomyFile) config.quality.consistency.taxonomyFile = taxonomyFile
        if (organismTaxon) config.quality.consistency.organismTaxon = organismTaxon
        // Inference needs the constraints loaded even when enforcement is off.
        if (inferTaxon && !organismTaxon) config.quality.consistency.taxonConstraints = true
        if (inferTaxonMinScore != null) config.quality.consistency.inferTaxonMinScore = inferTaxonMinScore
        if (enforceCompleteness) config.quality.completeness.enforce = true
        if (enforceCoherence) config.quality.coherence.enforce = true
        if (noProvenance) config.quality.provenance = false

        applyPredictorFlags(config)

        // Wall-clock timing per phase — written to <genome>_timing.tsv and echoed,
        // so genome-scale runs report where the time went (prediction vs ontology
        // load vs metrics vs enforcement).
        def timings = new LinkedHashMap<String, Double>()
        final long tStart = System.nanoTime()
        long[] lastMark = [tStart]
        def lap = { String name ->
            long now = System.nanoTime()
            timings[name] = (timings[name] ?: 0.0d) + (now - lastMark[0]) / 1_000_000_000.0d
            lastMark[0] = now
        }

        // Build and configure annotation pipeline
        def pipeline = new AnnotationPipeline(config)
        pipeline.configure()
        lap('setup')

        // Run annotation
        def genome = pipeline.annotate(primary)
        lap('prediction')

        // Optionally prune the low-confidence GO-DAG tail before anything else
        // touches it. DeepGO-style predictors propagate every prediction up the
        // ontology, emitting thousands of near-zero terms per protein; on a
        // complete genome that is ~1e6 annotations that swamp IO, metrics and the
        // JSON the web shows. A small floor keeps runs tractable without losing
        // the confident calls (taxon inference, which uses its own high threshold,
        // is unaffected).
        if (minScore != null && minScore > 0.0d) {
            long before = genome.proteins.sum { (long) it.annotations.annotations.size() } ?: 0L
            genome.proteins.each { p ->
                p.annotations.annotations.removeAll { it.isGO() && it.score < minScore }
            }
            long after = genome.proteins.sum { (long) it.annotations.annotations.size() } ?: 0L
            println "Pruned GO annotations below score ${minScore}: ${before} -> ${after}"
            lap('filter')
        }

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
        lap('output')

        // Genome-scale metrics (+ optional enforcement). Needs a GO OWL file;
        // metrics run per contig by default ("not across").
        boolean wantMetrics = config.quality.enabled && goOwl
        boolean wantEnforce = (config.quality.consistency.enforce || config.quality.completeness.enforce
                || config.quality.coherence.enforce) && goOwl
        if ((wantMetrics || wantEnforce) && genome.allGoTerms().size() > 0) {
            println ""
            def qp = new QualityPipeline(config)
                .goOwlFile(goOwl)
                .essentialFunctionsForDomain(config.resolveOrganismDomain())
            if (complexTermsFile) qp.complexTermsFile(new File(complexTermsFile))
            // Coherence enforcement (process has_part) needs the ELK reasoner.
            def qualityPipeline = config.quality.coherence.enforce ? qp.initialize() : qp.initializeLite()
            def goOnt = qualityPipeline.goOntology
            lap('ontology_load')

            // Infer the organism taxon from the predictions (before enforcement, so
            // enforcement uses the inferred taxon) unless one was asserted explicitly.
            if (inferTaxon && !config.quality.consistency.organismTaxon) {
                def inf = qualityPipeline.inferOrganismTaxon(genome)
                writeTaxonInference(inf, outputDir, genome.id)
                lap('taxon_inference')
            }

            if (wantEnforce) {
                println "Enforcing quality constraints (consistency=${config.quality.consistency.enforce}/" +
                    "${config.quality.consistency.enforceMode}, completeness=${config.quality.completeness.enforce}, " +
                    "coherence=${config.quality.coherence.enforce})..."
                def report = qualityPipeline.enforceAll(genome)
                if (report.provenance && report.count() > 0) {
                    // Label every acted-upon term so the log never shows a bare GO id.
                    if (goOnt != null) report.actions.each { it.termLabel = goOnt.getLabel(it.term) }
                    def actionsFile = new File(outputDir, "${genome.id}_enforcement_actions.tsv")
                    report.writeTsv(actionsFile)
                    println "  ${report.count()} enforcement action(s) ${report.countsByDimension()}: ${actionsFile}"
                }
                // Re-persist annotations so the written outputs reflect enforcement.
                pipeline.writeOutput(genome, outputDir)
                lap('enforcement')
            }

            // Enrich each GO annotation with its label + aspect (MF/BP/CC) and
            // re-persist, so the annotations TSV never shows a bare GO id and
            // carries an `aspect` column for standard-format (GAF) export /
            // per-aspect summaries. After enforcement, so imputed terms are covered.
            if (goOnt != null) {
                genome.proteins.each { p ->
                    p.annotations.goAnnotations().each { a ->
                        if (!a.goAspect) a.goAspect = goOnt.getAspect(a.value)
                        if (!a.goLabel) a.goLabel = goOnt.getLabel(a.value)
                    }
                }
                pipeline.writeOutput(genome, outputDir)
            }

            if (wantMetrics) {
                String scope = config.quality.scope ?: 'contig'
                println "Running genome-scale metrics (scope=${scope})..."
                if (scope in ['contig', 'both']) {
                    def reports = qualityPipeline.evaluatePerContig(genome)
                    def tsv = new File(outputDir, "${genome.id}_quality_per_contig.tsv")
                    QualityReportWriter.writeSummaryTsv(reports, tsv)
                    println "  Per-contig metrics (${reports.size()} contigs): ${tsv}"
                    reports.each { r ->
                        def jf = new File(outputDir, "${r.genomeId.replaceAll('[:/]', '_')}_quality.json")
                        QualityReportWriter.writeJson(r, jf)
                    }
                }
                if (scope in ['genome', 'both']) {
                    def report = qualityPipeline.evaluate(genome)
                    def reportFile = new File(outputDir, "${genome.id}_quality.json")
                    QualityReportWriter.writeJson(report, reportFile)
                    println "  Whole-genome metrics: ${reportFile}"
                }
                lap('metrics')
            }

            qualityPipeline.dispose()
        }

        writeTimingReport(timings, (System.nanoTime() - tStart) / 1_000_000_000.0d,
            outputDir, genome, proteinsOnly)

        println ""
        println "Done. Output in: ${outputDir}"
    }

    /**
     * Persist per-phase wall-clock timing to {@code <genome>_timing.tsv}
     * (phase, seconds, proteins, seconds_per_1k_proteins) and echo a summary, so
     * genome-scale runs are transparent about where the time was spent.
     */
    private void writeTimingReport(Map<String, Double> timings, double totalSec,
                                   File outputDir, gspa.model.Genome genome, boolean proteinsOnly) {
        if (timings.isEmpty()) return
        outputDir.mkdirs()
        int nProteins = genome?.proteins?.size() ?: 0
        def file = new File(outputDir, "${genome.id}_timing.tsv")
        file.withWriter { w ->
            w.writeLine(['phase', 'seconds', 'percent'].join('\t'))
            timings.each { name, sec ->
                String pct = totalSec > 0 ? String.format(Locale.ROOT, '%.1f', 100.0d * sec / totalSec) : '0.0'
                w.writeLine([name, String.format(Locale.ROOT, '%.3f', sec), pct].join('\t'))
            }
            w.writeLine(['total', String.format(Locale.ROOT, '%.3f', totalSec), '100.0'].join('\t'))
            w.writeLine(['proteins', nProteins.toString(), ''].join('\t'))
        }
        double per1k = nProteins > 0 ? totalSec / nProteins * 1000.0d : 0.0d
        println ""
        println "Timing (${nProteins} proteins, ${String.format(Locale.ROOT, '%.1f', totalSec)}s total" +
            (nProteins > 0 ? ", ${String.format(Locale.ROOT, '%.1f', per1k)}s / 1k proteins" : '') + "):"
        timings.each { name, sec ->
            println "  ${name.padRight(16)} ${String.format(Locale.ROOT, '%6.2f', sec)}s" +
                (totalSec > 0 ? " (${String.format(Locale.ROOT, '%.0f', 100.0d * sec / totalSec)}%)" : '')
        }
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

    /**
     * Write the taxon-inference evidence to {@code <genome>_taxon_inference.tsv}
     * (taxon, label, forbidden, support, inferred) and echo a one-line summary.
     */
    private void writeTaxonInference(gspa.metrics.TaxonInference.Result inf, File outputDir, String genomeId) {
        if (inf == null) return
        outputDir.mkdirs()
        def file = new File(outputDir, "${genomeId}_taxon_inference.tsv")
        file.withWriter { w ->
            w.writeLine(['taxon', 'label', 'depth', 'forbidden_terms', 'support_terms',
                         'inferred', 'lineage'].join('\t'))
            inf.candidates.each { c ->
                boolean isBest = (c.taxon == inf.taxon)
                w.writeLine([c.taxon, c.label, c.depth.toString(), c.forbidden.toString(),
                             c.support.toString(), isBest.toString(),
                             isBest ? inf.lineage.join(' > ') : ''].join('\t'))
            }
        }
        if (inf.taxon) {
            String detail = inf.candidates.take(6).collect { "${it.label}=${it.forbidden}" }.join(' ')
            println "Inferred organism taxon: ${inf.label} (${inf.taxon})" +
                "${inf.confident ? '' : ' [low confidence]'} from ${inf.constrainedPresent} " +
                "constraint-bearing terms"
            if (inf.lineage) println "  Lineage: ${inf.lineage.join(' > ')}"
            println "  Top candidates (violations): ${detail}"
        } else {
            println "Organism taxon could not be inferred (no constraint-bearing terms predicted)."
        }
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

        if (basePredictor) {
            config.predictors.neural.basePredictor = basePredictor
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

        if (deepGoPlusPlusLightEnabled) {
            config.predictors.neural.deepGoPlusPlusLight.enabled = true
        }
        if (deepGoPlusPlusLightAssets) {
            config.predictors.neural.deepGoPlusPlusLight.assets = deepGoPlusPlusLightAssets
        }
        if (deepGoPlusPlusLightModels) {
            config.predictors.neural.deepGoPlusPlusLight.modelsDir = deepGoPlusPlusLightModels
        }
        if (deepGoPlusPlusLightInterpro) {
            config.predictors.neural.deepGoPlusPlusLight.interpro = true
        }
        if (deepGoPlusPlusLightCnn) {
            config.predictors.neural.deepGoPlusPlusLight.cnn = true
        }
        if (deepGoPlusPlusLightInterproscan) {
            config.predictors.neural.deepGoPlusPlusLight.interproscan = deepGoPlusPlusLightInterproscan
        }

        // Push the genome-scale score floor (--min-score) down into the neural
        // sidecar so it never computes/emits the sub-floor, DAG-propagated tail
        // that gspa would otherwise parse, propagate and then discard in the
        // post-prediction filter below. On a bacterial proteome this shrinks the
        // emitted set ~6x (170k -> 26k rows at 0.5) and cuts the prediction phase
        // ~3x (151s -> 55s) with identical surviving annotations. The filter
        // (further down) stays as a safety net.
        if (minScore != null && minScore > 0.0d) {
            def lc = config.predictors.neural.deepGoPlusPlusLight
            lc.minScore = Math.max(lc.minScore, minScore)
        }
        // Wire the -t/--threads budget into the DIAMOND/sidecar thread count.
        // (The option was previously parsed but never consumed by any predictor.)
        if (threads > 0) {
            config.predictors.neural.deepGoPlusPlusLight.threads = threads
        }
    }
}
