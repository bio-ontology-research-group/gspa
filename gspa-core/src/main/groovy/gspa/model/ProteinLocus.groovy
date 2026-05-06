package gspa.model

import groovy.transform.Canonical
import groovy.transform.CompileStatic

/**
 * A protein's position on a genome.
 *
 * <p>Part of the Phase 12 Reaction-Local Genomic Context framework. The
 * {@code genomeId} field is optional in single-genome mode (may be null
 * or empty), but carries the cross-genome identifier in multi-genome
 * flows — same shape as {@link ProteinRef} so cross-panel operations
 * are data-only.</p>
 *
 * <p>Coordinates are 1-based inclusive (GFF convention). {@code start} &le; {@code end};
 * use {@link #length()} for the CDS span and {@link #midpoint()} for the
 * centre position used by the genomic density kernel.</p>
 */
@Canonical
@CompileStatic
class ProteinLocus {
    String genomeId
    String proteinId
    String contig
    int start
    int end
    Strand strand

    int length() {
        end - start + 1
    }

    /** Coordinate used by kernel density: midpoint of the CDS. */
    int midpoint() {
        (start + end) / 2 as int
    }
}
