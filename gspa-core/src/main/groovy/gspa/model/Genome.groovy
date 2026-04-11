package gspa.model

import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * Top-level representation of a genome or metagenome-assembled genome (MAG).
 * Contains contigs, proteins, taxonomy, and assembly metadata.
 */
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class Genome {

    /** Genome identifier */
    String id

    /** Contigs in this genome */
    List<Contig> contigs = []

    /** Taxonomy information */
    TaxonomyInfo taxonomy

    /** Assembly quality metadata */
    AssemblyInfo assemblyInfo

    /** Quality assessment report (populated after evaluation) */
    QualityReport qualityReport

    /** Free-form metadata */
    Map<String, Object> metadata = [:]

    /** Whether this is a MAG */
    boolean mag = false

    /** All proteins across all contigs */
    List<Protein> getProteins() {
        contigs.collectMany { it.proteins }
    }

    /** All features across all contigs */
    List<Feature> getFeatures() {
        contigs.collectMany { it.features }
    }

    /** Total assembly size in bp */
    long getAssemblySize() {
        if (contigs.isEmpty()) return 0L
        contigs.sum { Contig c -> (long) c.length } as long
    }

    /** Number of proteins */
    int getProteinCount() {
        if (contigs.isEmpty()) return 0
        contigs.sum { Contig c -> c.proteins.size() } as int
    }

    /** N50 statistic */
    int getN50() {
        if (contigs.isEmpty()) return 0
        def lengths = contigs.collect { it.length }.sort().reverse()
        long half = assemblySize / 2
        long cumulative = 0
        for (int len : lengths) {
            cumulative += len
            if (cumulative >= half) return len
        }
        return lengths.last()
    }

    /** Add a contig */
    void addContig(Contig contig) {
        contig.genome = this
        contigs << contig
    }

    /** Find a contig by ID */
    Contig findContig(String contigId) {
        contigs.find { it.id == contigId }
    }

    /** Find a protein by ID */
    Protein findProtein(String proteinId) {
        proteins.find { it.id == proteinId }
    }

    /**
     * Collect all GO annotations across all proteins.
     * Returns a map of protein ID -> set of GO term IDs.
     */
    Map<String, Set<String>> proteinGoAnnotations() {
        proteins.collectEntries { p ->
            [(p.id): p.annotations.goTermIds()]
        }
    }

    /**
     * Get all unique GO terms annotated in this genome.
     */
    Set<String> allGoTerms() {
        proteins.collectMany { it.annotations.goTermIds() } as Set
    }

    @Override
    String toString() {
        "Genome(${id}, ${contigs.size()} contigs, ${proteinCount} proteins, ${assemblySize}bp)"
    }
}
