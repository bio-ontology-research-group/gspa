package gspa.predictor.genecalling

import gspa.model.Genome

/**
 * Interface for gene calling tools.
 * Gene callers predict CDS features and extract protein sequences from nucleotide input.
 */
interface GeneCaller {

    String getName()
    boolean isAvailable()

    /**
     * Call genes on a genome.
     * Populates the genome with CDS features and corresponding proteins.
     */
    void callGenes(Genome genome)
}
