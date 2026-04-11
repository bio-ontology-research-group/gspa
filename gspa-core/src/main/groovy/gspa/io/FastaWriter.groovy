package gspa.io

import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein

/**
 * Writes FASTA format files.
 */
class FastaWriter {

    static final int DEFAULT_LINE_WIDTH = 70

    /** Write genome contigs as nucleotide FASTA */
    static void writeGenome(Genome genome, File output, int lineWidth = DEFAULT_LINE_WIDTH) {
        output.withWriter { writer ->
            genome.contigs.each { contig ->
                writeEntry(writer, contig.id, contig.sequence, lineWidth)
            }
        }
    }

    /** Write proteins as protein FASTA */
    static void writeProteins(List<Protein> proteins, File output, int lineWidth = DEFAULT_LINE_WIDTH) {
        output.withWriter { writer ->
            proteins.each { protein ->
                writeEntry(writer, protein.id, protein.sequence, lineWidth)
            }
        }
    }

    /** Write a single FASTA entry */
    static void writeEntry(Writer writer, String id, String sequence, int lineWidth = DEFAULT_LINE_WIDTH) {
        writer.writeLine(">${id}")
        for (int i = 0; i < sequence.length(); i += lineWidth) {
            writer.writeLine(sequence.substring(i, Math.min(i + lineWidth, sequence.length())))
        }
    }
}
