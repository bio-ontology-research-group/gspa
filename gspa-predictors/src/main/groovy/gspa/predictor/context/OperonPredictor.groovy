package gspa.predictor.context

import gspa.model.*
import gspa.predictor.GenomePredictor
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Operon predictor based on genomic context analysis.
 *
 * Identifies operons by analyzing:
 * 1. Intergenic distance between adjacent genes (short = likely co-operonic)
 * 2. Strand orientation (same strand required for operon membership)
 * 3. Optionally: conservation of gene order across related genomes
 *
 * Once operons are identified, transfers functional annotations between
 * co-operonic genes under the assumption that genes in the same operon
 * participate in related biological processes.
 */
class OperonPredictor implements GenomePredictor {

    private static final Logger log = LoggerFactory.getLogger(OperonPredictor)

    /** Maximum intergenic distance (bp) to consider genes co-operonic */
    int maxIntergenicDistance = 300

    /** Require genes to be on the same strand */
    boolean requireSameStrand = true

    /** Minimum operon size (number of genes) to trigger function transfer */
    int minOperonSize = 2

    /** Confidence score for operon-transferred annotations */
    double transferScore = 0.4

    @Override
    String getName() { 'operon' }

    @Override
    String getVersion() { '1.0' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.GO] as Set
    }

    @Override
    boolean isAvailable() { true }  // No external tool needed

    @Override
    List<Annotation> predict(Protein protein) {
        []  // Requires genome context; use predictGenome instead
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        log.info("Predicting operons for ${genome.id}")

        // Detect operons
        List<Operon> operons = detectOperons(genome)
        log.info("Detected ${operons.size()} operons (${operons.sum { it.genes.size() } ?: 0} genes total)")

        // Transfer annotations between co-operonic genes
        Map<String, List<Annotation>> transfers = transferAnnotations(operons)
        log.info("Generated ${transfers.values().sum { it.size() } ?: 0} annotation transfers")

        transfers
    }

    /**
     * Detect operons across all contigs in a genome.
     * Returns a list of operons, each containing ordered gene lists.
     */
    List<Operon> detectOperons(Genome genome) {
        List<Operon> operons = []

        genome.contigs.each { contig ->
            def proteins = contig.proteins.sort { it.start }
            if (proteins.size() < 2) return

            List<Protein> currentOperon = [proteins[0]]

            for (int i = 1; i < proteins.size(); i++) {
                def prev = proteins[i - 1]
                def curr = proteins[i]

                boolean coOperonic = isCoOperonic(prev, curr)

                if (coOperonic) {
                    currentOperon << curr
                } else {
                    if (currentOperon.size() >= minOperonSize) {
                        operons << new Operon(
                            contigId: contig.id,
                            genes: new ArrayList<>(currentOperon)
                        )
                    }
                    currentOperon = [curr]
                }
            }
            // Don't forget the last operon
            if (currentOperon.size() >= minOperonSize) {
                operons << new Operon(
                    contigId: contig.id,
                    genes: new ArrayList<>(currentOperon)
                )
            }
        }

        operons
    }

    /**
     * Check if two adjacent genes are likely co-operonic.
     */
    boolean isCoOperonic(Protein a, Protein b) {
        // Same strand check
        if (requireSameStrand && !a.strand.isSameStrand(b.strand)) {
            return false
        }

        // Intergenic distance check
        int distance = intergenicDistance(a, b)
        if (distance > maxIntergenicDistance) {
            return false
        }

        // Overlapping genes on same strand are co-operonic
        if (distance < 0 && a.strand.isSameStrand(b.strand)) {
            return true
        }

        true
    }

    /**
     * Transfer GO annotations between co-operonic genes.
     * For each operon, collect all BP (biological process) annotations
     * and transfer them to genes that lack them.
     */
    Map<String, List<Annotation>> transferAnnotations(List<Operon> operons) {
        Map<String, List<Annotation>> transfers = [:].withDefault { [] }

        operons.each { operon ->
            // Collect all GO BP annotations in the operon
            Set<String> operonBPTerms = []
            Map<String, Set<String>> proteinTerms = [:]

            operon.genes.each { protein ->
                def bpTerms = protein.annotations.goAnnotations()
                    .findAll { it.goAspect == 'BP' || it.goAspect == null }
                    .collect { it.value } as Set
                proteinTerms[protein.id] = bpTerms
                operonBPTerms.addAll(bpTerms)
            }

            // Transfer: each gene gets BP terms it doesn't already have
            operon.genes.each { protein ->
                Set<String> existing = proteinTerms[protein.id] ?: [] as Set
                Set<String> toTransfer = operonBPTerms - existing

                toTransfer.each { goTerm ->
                    // Find which co-operonic gene had this term
                    String sourceGene = operon.genes.find { g ->
                        g.id != protein.id && (proteinTerms[g.id]?.contains(goTerm))
                    }?.id

                    transfers[protein.id] << new Annotation(
                        type: AnnotationType.GO,
                        value: goTerm,
                        source: name,
                        score: transferScore,
                        evidence: 'IGC',  // Inferred from Genomic Context
                        goAspect: 'BP',
                        metadata: [
                            operonContig: operon.contigId,
                            operonSize: operon.genes.size(),
                            transferredFrom: sourceGene,
                        ]
                    )
                }
            }
        }

        transfers
    }

    private static int intergenicDistance(Protein a, Protein b) {
        if (a.end < b.start) return b.start - a.end - 1
        if (b.end < a.start) return a.start - b.end - 1
        return -(Math.min(a.end, b.end) - Math.max(a.start, b.start) + 1)
    }
}

/**
 * Represents a predicted operon: a set of co-transcribed genes.
 */
class Operon {
    String contigId
    List<Protein> genes = []

    int getStart() { genes.first().start }
    int getEnd() { genes.last().end }
    int getSize() { genes.size() }
    Strand getStrand() { genes.first().strand }

    @Override
    String toString() {
        "Operon(${contigId}:${start}-${end}, ${size} genes, ${strand.symbol})"
    }
}
