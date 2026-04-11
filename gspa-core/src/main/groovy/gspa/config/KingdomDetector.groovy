package gspa.config

import gspa.model.Genome
import gspa.model.OrganismDomain
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Auto-detects organism kingdom from genome properties.
 *
 * Uses heuristics based on:
 * - Genome size
 * - GC content
 * - Gene density (if features are present)
 * - Coding density
 *
 * This is a rough classifier — for accurate taxonomy, use GTDB-Tk.
 */
class KingdomDetector {

    private static final Logger log = LoggerFactory.getLogger(KingdomDetector)

    /**
     * Detect the likely organism domain from genome properties.
     */
    static OrganismDomain detect(Genome genome) {
        long assemblySize = genome.assemblySize
        if (assemblySize == 0) {
            log.warn("Cannot auto-detect kingdom: no sequence data")
            return OrganismDomain.UNKNOWN
        }

        double gcContent = computeGcContent(genome)
        int contigCount = genome.contigs.size()
        double geneDensity = genome.proteinCount > 0 ?
            genome.proteinCount / (assemblySize / 1_000_000.0) : 0.0

        log.info("Auto-detecting kingdom: size=${assemblySize}bp, GC=${String.format('%.1f%%', gcContent * 100)}, " +
            "contigs=${contigCount}, gene_density=${String.format('%.0f', geneDensity)} genes/Mbp")

        // Virus: very small genomes
        if (assemblySize < 50_000) {
            log.info("Detected: VIRUS (genome < 50kb)")
            return OrganismDomain.VIRUS
        }

        // Eukaryote: large genomes (> 10 Mbp for unicellular, much larger for multicellular)
        if (assemblySize > 15_000_000) {
            log.info("Detected: EUKARYA (genome > 15Mbp)")
            return OrganismDomain.EUKARYA
        }

        // Between 50kb and 15Mbp: prokaryote
        // Distinguish bacteria vs archaea is harder without phylogenetic analysis
        // Use a default of bacteria (most common case)
        // Extreme GC or very specific sizes might hint at archaea but unreliable

        if (assemblySize < 500_000) {
            // Very small prokaryote — could be obligate symbiont
            log.info("Detected: BACTERIA (small genome, likely symbiont/parasite)")
            return OrganismDomain.BACTERIA
        }

        // Standard prokaryote range (0.5 - 15 Mbp)
        // Gene density can help: prokaryotes typically have 800-1200 genes/Mbp
        // Eukaryotes that slip into this range (e.g., Encephalitozoon at 2.3Mbp)
        // would have much lower gene density

        if (geneDensity > 0 && geneDensity < 400) {
            // Low gene density suggests eukaryote even at small genome size
            log.info("Detected: EUKARYA (low gene density: ${String.format('%.0f', geneDensity)} genes/Mbp)")
            return OrganismDomain.EUKARYA
        }

        log.info("Detected: BACTERIA (default for prokaryotic-sized genome)")
        return OrganismDomain.BACTERIA
    }

    /**
     * Compute GC content of the genome.
     */
    static double computeGcContent(Genome genome) {
        long gc = 0
        long total = 0
        genome.contigs.each { contig ->
            if (contig.sequence) {
                contig.sequence.toUpperCase().each { ch ->
                    if (ch == 'G' || ch == 'C') gc++
                    if (ch in ['A', 'T', 'G', 'C']) total++
                }
            }
        }
        total > 0 ? gc / (double) total : 0.0
    }
}
