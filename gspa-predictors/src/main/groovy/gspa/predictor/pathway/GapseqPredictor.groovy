package gspa.predictor.pathway

import gspa.model.*
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor

/**
 * gapseq predictor for metabolic pathway reconstruction.
 * Identifies metabolic reactions and pathways, performs gap-filling
 * to create draft genome-scale metabolic models.
 *
 * Supports community mode for crossfeeding analysis.
 */
class GapseqPredictor extends AbstractToolPredictor implements GenomePredictor {

    /** gapseq medium file for gap-filling */
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
    String getName() { 'gapseq' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.EC, AnnotationType.KEGG] as Set
    }

    @Override
    String getExecutable() { 'gapseq' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def outputPrefix = new File(outputDir, 'gapseq_model')
        def cmd = [
            executablePath ?: executable,
            'find',
            '-p', 'all',
            '-b', taxonomy,
            inputFasta.absolutePath,
        ]
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // gapseq outputs reaction lists and pathway completeness
        def reactionsFile = findOutputFile(outputDir, '*-Reactions.tbl')
        if (reactionsFile == null) return [:]

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

        def tmpDir = File.createTempDir("gspa_gapseq_", '')
        try {
            File inputFasta = writeTargetFasta(genome, tmpDir)

            def command = buildCommand(inputFasta, tmpDir)
            log.info("Running gapseq (target=${target}): ${command.join(' ')}")
            def result = execute(command, tmpDir)

            if (result.exitCode != 0) {
                log.error("gapseq failed: ${result.stderr}")
                return [:]
            }

            return parseOutput(tmpDir)
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
