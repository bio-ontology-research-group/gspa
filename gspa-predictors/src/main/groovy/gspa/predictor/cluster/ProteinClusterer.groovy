package gspa.predictor.cluster

import gspa.model.Genome
import gspa.model.ProteinClusterSet

/**
 * Groups proteins into equivalence clusters by sequence similarity.
 *
 * Part 1 implementation ({@link IntragenomeClusterer}) accepts a
 * single-element list of genomes. Phase 11 introduces a cross-genome
 * implementation that accepts any number of genomes and clusters their
 * proteomes jointly; the interface is identical.
 */
interface ProteinClusterer {

    /**
     * Cluster proteins across the given genomes.
     *
     * @param genomes   one (Part 1) or more (Phase 11) genomes to cluster
     * @param identity  minimum sequence identity for cluster membership
     *                  (e.g., 0.9 = 90%)
     * @param coverage  minimum alignment coverage (e.g., 0.8 = 80%)
     * @return {@link ProteinClusterSet} covering every protein in every
     *         input genome (singleton clusters for orphans)
     */
    ProteinClusterSet cluster(List<Genome> genomes, double identity, double coverage)
}
