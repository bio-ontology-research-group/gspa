package gspa.model

import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * Represents a contiguous sequence (contig, chromosome, scaffold)
 * within a genome assembly.
 */
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class Contig {

    /** Contig identifier (matches GFF3 seqid) */
    String id

    /** Nucleotide sequence (may be null for lazy loading) */
    String sequence

    /** All features on this contig, sorted by position */
    List<Feature> features = []

    /** All proteins encoded on this contig */
    List<Protein> proteins = []

    /** Back-pointer to the genome */
    Genome genome

    int getLength() {
        sequence?.length() ?: 0
    }

    /** Get features of a specific type */
    List<Feature> featuresOfType(FeatureType type) {
        features.findAll { it.type == type }
    }

    /** Get CDS features */
    List<Feature> cdsFeatures() {
        featuresOfType(FeatureType.CDS)
    }

    /** Add a feature and keep the list sorted */
    void addFeature(Feature feature) {
        feature.seqId = this.id
        features << feature
    }

    /** Add a protein */
    void addProtein(Protein protein) {
        protein.contig = this
        proteins << protein
    }

    /** Sort features by genomic position */
    void sortFeatures() {
        features.sort()
    }

    /**
     * Get pairs of adjacent proteins on the same strand,
     * useful for operon prediction.
     */
    List<List<Protein>> adjacentProteinPairs() {
        def sorted = proteins.sort { it.start }
        (0..<(sorted.size() - 1)).collect { i ->
            [sorted[i], sorted[i + 1]]
        }
    }

    /**
     * Intergenic distance between two adjacent features.
     * Returns negative if features overlap.
     */
    static int intergenicDistance(Feature a, Feature b) {
        if (a.end < b.start) return b.start - a.end - 1
        if (b.end < a.start) return a.start - b.end - 1
        return -(Math.min(a.end, b.end) - Math.max(a.start, b.start) + 1)
    }

    @Override
    String toString() {
        "Contig(${id}, ${length}bp, ${proteins.size()} proteins)"
    }
}
