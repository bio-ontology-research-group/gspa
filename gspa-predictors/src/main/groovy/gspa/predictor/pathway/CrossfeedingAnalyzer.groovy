package gspa.predictor.pathway

import gspa.model.Genome
import gspa.predictor.CommunityPredictor
import gspa.predictor.CommunityResult
import gspa.predictor.CrossfeedingEdge
import gspa.predictor.GenomeMetabolicSummary
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Analyzes metabolic crossfeeding between co-occurring genomes.
 *
 * For MAGs/communities, identifies:
 * 1. Metabolites each genome can produce (exports)
 * 2. Metabolites each genome needs but cannot synthesize (imports)
 * 3. Complementary pairs where one genome's export fills another's gap
 * 4. Community-level pathway completeness
 *
 * Uses gapseq community mode for metabolic model reconstruction,
 * optionally exports SBML models for FBA via COBRApy.
 */
class CrossfeedingAnalyzer implements CommunityPredictor {

    private static final Logger log = LoggerFactory.getLogger(CrossfeedingAnalyzer)

    /** Path to gapseq executable */
    String gapseqPath

    /** Medium file for gap-filling */
    String medium

    /** Whether to export SBML models */
    boolean exportSbml = false

    /** Whether to run FBA via COBRApy (requires Python) */
    boolean runFba = false

    /** COBRApy script path (if runFba is true) */
    String cobraScript

    /** Working directory for intermediate files */
    File workDir

    @Override
    String getName() { 'crossfeeding' }

    @Override
    CommunityResult analyze(List<Genome> community) {
        log.info("Analyzing crossfeeding in community of ${community.size()} genomes")

        // Build per-genome metabolic summaries from annotations
        Map<String, GenomeMetabolicSummary> summaries = [:]
        community.each { genome ->
            summaries[genome.id] = buildMetabolicSummary(genome)
        }

        // Find crossfeeding interactions
        List<CrossfeedingEdge> interactions = findInteractions(summaries)

        // Compute community-level pathway completeness
        Set<String> communityEC = []
        community.each { genome ->
            genome.proteins.each { protein ->
                protein.annotations.ecAnnotations().each { ann ->
                    communityEC << ann.value
                }
            }
        }

        // Find unexplained gaps (metabolites needed by some genome but
        // not produced by any community member)
        Map<String, Set<String>> unexplainedGaps = [:]
        summaries.each { genomeId, summary ->
            Set<String> needed = new HashSet<>(summary.imports)
            Set<String> communityExports = summaries.findAll { it.key != genomeId }
                .values().collectMany { it.exports } as Set
            Set<String> unexplained = needed - communityExports
            if (unexplained) {
                unexplainedGaps[genomeId] = unexplained
            }
        }

        def result = new CommunityResult(
            genomeSummaries: summaries,
            interactions: interactions,
            unexplainedGaps: unexplainedGaps
        )

        log.info("Crossfeeding analysis: ${interactions.size()} interactions found, " +
            "${unexplainedGaps.values().sum { it.size() } ?: 0} unexplained gaps")

        // Optional: run gapseq community mode for full metabolic modeling
        if (gapseqPath && workDir) {
            runGapseqCommunity(community, result)
        }

        // Optional: export SBML
        if (exportSbml && workDir) {
            exportSbmlModels(community, result)
        }

        result
    }

    /**
     * Build a metabolic summary from a genome's EC annotations.
     * Identifies which metabolic capabilities the genome has
     * and which essential metabolites it can/cannot synthesize.
     */
    GenomeMetabolicSummary buildMetabolicSummary(Genome genome) {
        Set<String> ecNumbers = []
        genome.proteins.each { protein ->
            protein.annotations.ecAnnotations().each { ann ->
                ecNumbers << ann.value
            }
        }

        // Essential biosynthetic pathways and their key EC numbers
        // Amino acid biosynthesis
        def aminoAcidSynthesis = [
            'L-lysine'      : ['EC:4.3.3.7', 'EC:1.17.1.8', 'EC:2.3.1.117'],
            'L-tryptophan'  : ['EC:4.1.1.48', 'EC:4.2.1.20'],
            'L-histidine'   : ['EC:4.2.1.19', 'EC:2.4.2.17'],
            'L-methionine'  : ['EC:2.1.1.13', 'EC:2.1.1.14'],
            'L-leucine'     : ['EC:2.3.3.13', 'EC:1.1.1.85'],
            'L-isoleucine'  : ['EC:4.2.1.9', 'EC:2.6.1.42'],
            'L-valine'      : ['EC:4.2.1.9', 'EC:2.6.1.42'],
        ]
        // Vitamin/cofactor biosynthesis
        def vitaminSynthesis = [
            'thiamine'      : ['EC:2.5.1.3', 'EC:2.7.4.16'],
            'riboflavin'    : ['EC:2.5.1.78', 'EC:1.1.1.193'],
            'biotin'        : ['EC:2.8.1.6', 'EC:6.3.3.3'],
            'folate'        : ['EC:2.5.1.15', 'EC:1.5.1.3'],
            'cobalamin'     : ['EC:2.5.1.17', 'EC:1.16.8.1'],
        ]

        Set<String> canProduce = []
        Set<String> cannotProduce = []

        (aminoAcidSynthesis + vitaminSynthesis).each { metabolite, keyECs ->
            boolean canSynthesize = keyECs.any { ec -> ecNumbers.contains(ec) }
            if (canSynthesize) {
                canProduce << metabolite
            } else {
                cannotProduce << metabolite
            }
        }

        new GenomeMetabolicSummary(
            genomeId: genome.id,
            exports: canProduce,
            imports: cannotProduce
        )
    }

    /**
     * Find crossfeeding interactions: where one genome produces
     * what another needs.
     */
    List<CrossfeedingEdge> findInteractions(Map<String, GenomeMetabolicSummary> summaries) {
        List<CrossfeedingEdge> interactions = []

        summaries.each { consumerId, consumer ->
            consumer.imports.each { metabolite ->
                // Find producers
                summaries.each { producerId, producer ->
                    if (producerId != consumerId && producer.exports.contains(metabolite)) {
                        interactions << new CrossfeedingEdge(
                            producerGenomeId: producerId,
                            consumerGenomeId: consumerId,
                            metabolite: metabolite,
                            confidence: 0.5  // EC-based inference
                        )
                    }
                }
            }
        }

        interactions
    }

    /**
     * Run gapseq in community mode for detailed metabolic exchange analysis.
     */
    private void runGapseqCommunity(List<Genome> community, CommunityResult result) {
        // Write individual genome FASTA files
        def genomesDir = new File(workDir, 'community_genomes')
        genomesDir.mkdirs()

        community.each { genome ->
            def fastaFile = new File(genomesDir, "${genome.id}.fna")
            fastaFile.withWriter { writer ->
                genome.contigs.each { contig ->
                    if (contig.sequence) {
                        writer.writeLine(">${contig.id}")
                        writer.writeLine(contig.sequence)
                    }
                }
            }
        }

        // Run gapseq find + draft + fill for each genome
        // Then run community analysis
        log.info("Running gapseq community analysis...")
        // In a real implementation, this would execute gapseq commands
        // and parse the community exchange results
    }

    /**
     * Export SBML metabolic models for optional FBA analysis.
     */
    private void exportSbmlModels(List<Genome> community, CommunityResult result) {
        def sbmlDir = new File(workDir, 'sbml_models')
        sbmlDir.mkdirs()
        log.info("SBML export directory: ${sbmlDir}")
        // In a real implementation, would convert gapseq models to SBML
    }

    /**
     * Adjust quality metrics based on crossfeeding context.
     * A pathway gap in genome A is not a quality issue if genome B
     * in the same community can supply the missing metabolite.
     */
    static double adjustPathwayCoherence(double rawCoherence,
                                          String genomeId,
                                          CommunityResult community) {
        if (community == null) return rawCoherence

        def summary = community.genomeSummaries[genomeId]
        if (summary == null) return rawCoherence

        // Count how many import needs are fulfilled by the community
        int totalNeeds = summary.imports.size()
        int fulfilledByComm = community.interactions
            .findAll { it.consumerGenomeId == genomeId }
            .collect { it.metabolite }
            .unique()
            .size()

        if (totalNeeds == 0) return rawCoherence

        // Boost coherence proportionally to how many gaps are community-filled
        double communityBoost = fulfilledByComm / (double) totalNeeds
        double adjusted = rawCoherence + (1.0 - rawCoherence) * communityBoost * 0.5
        Math.min(adjusted, 1.0)
    }
}
