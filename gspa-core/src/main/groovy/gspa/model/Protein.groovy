package gspa.model

import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * Represents a protein encoded by a CDS feature in a genome.
 * Central entity for function annotation.
 */
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class Protein {

    /** Unique protein identifier */
    String id

    /** Amino acid sequence */
    String sequence

    /** The CDS feature that encodes this protein */
    Feature sourceFeature

    /** Back-pointer to the contig this protein is on */
    Contig contig

    /** Functional annotations from all sources */
    AnnotationSet annotations = new AnnotationSet()

    /** Path to predicted 3D structure file (PDB/mmCIF), if available */
    String structurePath

    int getLength() {
        sequence?.length() ?: 0
    }

    /** Strand of the encoding CDS */
    Strand getStrand() {
        sourceFeature?.strand ?: Strand.NONE
    }

    /** Start position on the contig */
    int getStart() {
        sourceFeature?.start ?: 0
    }

    /** End position on the contig */
    int getEnd() {
        sourceFeature?.end ?: 0
    }

    /** Get the genome this protein belongs to (via contig) */
    Genome getGenome() {
        contig?.genome
    }

    /** FASTA-formatted representation */
    String toFasta() {
        ">${id}\n${sequence}"
    }

    @Override
    String toString() {
        "Protein(${id}, ${length}aa, ${annotations.size()} annotations)"
    }
}
