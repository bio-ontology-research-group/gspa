package gspa.model

import groovy.transform.Canonical

/**
 * Reference to a protein by (genomeId, proteinId).
 *
 * Phase 10 Part 1 populates genomeId with the single genome being processed;
 * Phase 11's cross-genome work uses distinct genomeIds in one ProteinClusterSet.
 * The shape is identical so cross-genome becomes a data-only change.
 */
@Canonical
class ProteinRef {
    String genomeId
    String proteinId
}
