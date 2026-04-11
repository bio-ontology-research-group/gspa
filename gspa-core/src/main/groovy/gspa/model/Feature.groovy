package gspa.model

import groovy.transform.Canonical
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * A genomic feature in GFF3-compatible representation.
 * Can represent CDS, tRNA, rRNA, CRISPR arrays, etc.
 */
@Canonical
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class Feature implements Comparable<Feature> {

    /** Sequence ID (contig/chromosome) */
    String seqId

    /** Source of the feature (e.g., "prodigal", "gspa") */
    String source = 'gspa'

    /** Feature type (CDS, tRNA, rRNA, etc.) */
    FeatureType type

    /** 1-based start position (inclusive) */
    int start

    /** 1-based end position (inclusive) */
    int end

    /** Score (. if not applicable) */
    Double featureScore

    /** Strand: +, -, or . */
    Strand strand = Strand.NONE

    /** Reading frame (0, 1, 2 for CDS; . otherwise) */
    int phase = -1

    /** GFF3 attributes (ID, Name, product, etc.) */
    Map<String, List<String>> attributes = [:]

    int getLength() { end - start + 1 }

    String getId() {
        attributes['ID']?.first()
    }

    void setId(String id) {
        attributes['ID'] = [id]
    }

    String getName() {
        attributes['Name']?.first()
    }

    String getProduct() {
        attributes['product']?.first()
    }

    void setProduct(String product) {
        attributes['product'] = [product]
    }

    String getLocusTag() {
        attributes['locus_tag']?.first()
    }

    @Override
    int compareTo(Feature other) {
        int cmp = this.seqId <=> other.seqId
        if (cmp != 0) return cmp
        cmp = this.start <=> other.start
        if (cmp != 0) return cmp
        return this.end <=> other.end
    }

    /** Format as a GFF3 line */
    String toGff3() {
        String scoreStr = featureScore != null ? featureScore.toString() : '.'
        String strandStr = strand.symbol
        String phaseStr = phase >= 0 ? phase.toString() : '.'
        String attrsStr = attributes.collect { k, v ->
            "${k}=${v.join(',')}"
        }.join(';')
        "${seqId}\t${source}\t${type.gff3Type}\t${start}\t${end}\t${scoreStr}\t${strandStr}\t${phaseStr}\t${attrsStr}"
    }
}
