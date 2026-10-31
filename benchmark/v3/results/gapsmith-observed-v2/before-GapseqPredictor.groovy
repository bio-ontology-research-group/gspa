package gspa.predictor.pathway

import gspa.model.*
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor
import gspa.predictor.PredictorOutcome
import gspa.predictor.PredictorRunStatus
import gspa.predictor.ToolOutputContractException

/**
 * Pathway/gap predictor wrapping gapsmith (default) or gapseq.
 *
 * <p>gapsmith is a Rust reimplementation of gapseq with the same
 * output-file format (see github.com/bio-ontology-research-group/gapsmith).
 * The two binaries differ only in CLI flag spellings — the emitted
 * {@code *-Reactions.tbl} and {@code *-Pathways.tbl} files are
 * drop-in compatible, so downstream parsers (including gspa's
 * {@code parse_gapseq_gaps.py}) work unchanged.</p>
 *
 * <p>Switch via {@link #tool}: {@code 'gapsmith'} (default) or
 * {@code 'gapseq'}. {@code 'gapsmith'} is preferred because it's
 * significantly faster on large proteomes and can share a precomputed
 * alignment across genomes.</p>
 *
 * Supports community mode for crossfeeding analysis.
 */
class GapseqPredictor extends AbstractToolPredictor implements GenomePredictor {

    /**
     * External tool to invoke. {@code 'gapsmith'} (default) uses
     * {@code gapsmith find -p all -t Bacteria -o <out-dir>}.
     * {@code 'gapseq'} uses {@code gapseq find -p all -b Bacteria}
     * (legacy R/bash upstream).
     */
    String tool = 'gapsmith'

    /**
     * Reference-data directory for gapsmith ({@code --data-dir}).
     * When null, gapsmith auto-detects (assumes a sibling {@code data/}
     * dir, or {@code GAPSMITH_DATA_DIR}). Ignored for {@code tool=='gapseq'}.
     */
    String dataDir = null

    /** gapseq/gapsmith medium file for gap-filling */
    String medium

    /** Export SBML model */
    boolean exportSbml = false

    /** Taxonomy for pathway prediction */
    String taxonomy = 'Bacteria'

    /**
     * Phase 10: search target for gapseq's reaction-sequence database.
     * <ul>
     *   <li>{@code genome}  — tblastn against genome nucleotide (default, v1.0.0 behaviour)</li>
     *   <li>{@code proteome} — blastp against the full predicted proteome
     *       (~5× faster but loses small-ORF / frameshift rescue)</li>
     *   <li>{@code reps}    — blastp against cluster representatives only
     *       (requires intra-genome clustering to have run). A reaction is
     *       "found" if any cluster rep has a qualifying hit — equivalent
     *       to querying the whole proteome under the 90% identity
     *       equivalence, and the form Phase 11 will cache across genomes.</li>
     * </ul>
     */
    String target = 'genome'

    /**
     * When {@code target=='reps'}, the list of representative protein IDs
     * (from the intra-genome clusterer). {@code null} with {@code target=='reps'}
     * is an error.
     */
    Set<String> representativeIds = null

    @Override
    String getName() { tool }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.EC, AnnotationType.KEGG] as Set
    }

    @Override
    String getExecutable() {
        tool == 'gapseq' ? 'gapseq' : 'gapsmith'
    }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def exe = executablePath ?: getExecutable()
        if (tool == 'gapseq') {
            return [
                exe,
                'find',
                '-p', 'all',
                '-b', taxonomy,
                inputFasta.absolutePath,
            ]
        }
        // gapsmith: --data-dir is at the top level, out-dir goes to the subcommand
        def cmd = [exe]
        if (dataDir) cmd += ['--data-dir', dataDir]
        cmd += [
            'find',
            '-p', 'all',
            '-t', taxonomy,
            '-o', outputDir.absolutePath,
            inputFasta.absolutePath,
        ]
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // gapseq outputs reaction lists and pathway completeness
        // A successful gapseq/gapsmith find run writes <name>-Reactions.tbl with a header row.
        // A missing file, or the zero-byte table gapseq leaves behind when an internal step fails
        // while the shell still exits 0 (documented in tasks/todo.md, Phase C1), is a failed run.
        def reactionsFile = findOutputFile(outputDir, '*-Reactions.tbl')
        if (reactionsFile == null) {
            throw new ToolOutputContractException("${name}: no *-Reactions.tbl after a clean exit")
        }
        if (reactionsFile.length() == 0) {
            throw new ToolOutputContractException(
                "${name}: ${reactionsFile.name} is zero bytes (gapseq internal failure with exit 0); not 'no reactions'")
        }

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        reactionsFile.eachLine { line ->
            if (line.startsWith('#') || line.startsWith('Reaction') || line.trim().isEmpty()) return
            def fields = line.split('\t')
            if (fields.length < 6) return

            String reactionId = fields[0]
            String ecNumber = fields[2]
            String status = fields[4]      // active/inactive
            String geneName = fields.length > 5 ? fields[5] : ''

            if (ecNumber && ecNumber != '-' && status == 'active') {
                // Map EC to protein IDs based on gene names
                // gapseq reports gene locus tags
                if (geneName && geneName != '-') {
                    results[geneName] << new Annotation(
                        type: AnnotationType.EC,
                        value: ecNumber.startsWith('EC:') ? ecNumber : "EC:${ecNumber}",
                        source: name,
                        score: 0.7,
                        metadata: [reaction: reactionId, status: status]
                    )
                }
            }
        }
        results
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        if (target == 'reps' && (representativeIds == null || representativeIds.isEmpty())) {
            throw new IllegalStateException(
                "GapseqPredictor(target='reps') requires representativeIds to be populated " +
                "by an earlier intra-genome clustering step")
        }

        lastRunStatus = null
        def tmpDir = File.createTempDir("gspa_gapseq_", '')
        long t0 = System.currentTimeMillis()
        try {
            File inputFasta
            ContigExport export = null
            if (target == 'genome') {
                inputFasta = new File(tmpDir, 'genome.fna')
                export = writeContigFasta(genome, inputFasta)
                if (export.written.isEmpty()) {
                    recordStatus(PredictorRunStatus.failed(name,
                        'target=genome but no contig carries nucleotide sequence'), t0)
                    return [:]
                }
            } else {
                inputFasta = writeTargetFasta(genome, tmpDir)
                if (inputFasta.length() == 0) {
                    recordStatus(PredictorRunStatus.failed(name,
                        "target=${target} selected no protein sequences to search"), t0)
                    return [:]
                }
            }

            def result = runOrFail(buildCommand(inputFasta, tmpDir), tmpDir, t0)
            if (result == null) return [:]

            def parsed = parseGuarded(t0) { parseOutput(tmpDir) }
            int n = countAnnotations(parsed)
            if (export?.skipped) {
                recordStatus(PredictorRunStatus.partial(name, n, describeSkipped(export, genome.contigs.size()),
                    export.skipped, result.exitCode), t0)
            } else {
                recordStatus(PredictorRunStatus.completed(name, n, result.exitCode), t0)
            }
            return parsed
        } catch (Exception e) {
            if (lastRunStatus == null || lastRunStatus.outcome != PredictorOutcome.FAILED) {
                recordStatus(PredictorRunStatus.failed(name, "${e.class.simpleName}: ${e.message}"), t0)
            }
            throw e
        } finally {
            tmpDir.deleteDir()
        }
    }

    /** Emit the FASTA gapseq will search, based on the configured target. */
    private File writeTargetFasta(Genome genome, File dir) {
        switch (target) {
            case 'genome':
                def f = new File(dir, 'genome.fna')
                f.withWriter { w ->
                    genome.contigs.each { contig ->
                        if (contig.sequence) {
                            w.writeLine(">${contig.id}")
                            w.writeLine(contig.sequence)
                        }
                    }
                }
                return f
            case 'proteome':
                def f = new File(dir, 'proteins.faa')
                f.withWriter { w ->
                    genome.proteins.each { p ->
                        if (p.sequence) { w.writeLine(">${p.id}"); w.writeLine(p.sequence) }
                    }
                }
                return f
            case 'reps':
                def f = new File(dir, 'reps.faa')
                f.withWriter { w ->
                    genome.proteins.each { p ->
                        if (p.sequence && representativeIds.contains(p.id)) {
                            w.writeLine(">${p.id}"); w.writeLine(p.sequence)
                        }
                    }
                }
                return f
            default:
                throw new IllegalArgumentException("Unknown gapseq target: ${target}; expected genome|proteome|reps")
        }
    }

    /**
     * Run gapseq in community mode for crossfeeding analysis.
     * Takes multiple genome FASTA files and finds metabolic exchanges.
     */
    Map<String, Object> analyzeCommmunity(List<File> genomeFastas, File outputDir) {
        // gapseq community analysis would be done in a separate step
        // after individual models are built
        log.info("Community analysis with ${genomeFastas.size()} genomes")
        [:]
    }

    private File findOutputFile(File dir, String pattern) {
        dir.listFiles()?.find { it.name.endsWith('-Reactions.tbl') }
    }
}
