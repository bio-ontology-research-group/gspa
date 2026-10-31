package gspa.io

import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein

/**
 * Reads FASTA format files (nucleotide or protein sequences).
 * Supports standard FASTA with multi-line sequences.
 */
class FastaReader {

    /**
     * Read a nucleotide FASTA file into a Genome with contigs.
     */
    static Genome readGenome(File fastaFile, String genomeId = null) {
        def genome = new Genome(id: genomeId ?: fastaFile.name.replaceAll(/\.(fasta|fa|fna|fsa)(\.gz)?$/, ''))
        parseFasta(fastaFile).each { id, seq ->
            genome.addContig(new Contig(id: id, sequence: seq))
        }
        genome
    }

    /**
     * Read a protein FASTA file into a list of Proteins.
     */
    static List<Protein> readProteins(File fastaFile) {
        parseFasta(fastaFile).collect { id, seq ->
            new Protein(id: id, sequence: seq)
        }
    }

    /**
     * Parse a FASTA file into an ordered map of ID -> sequence.
     * Handles multi-line sequences and comment lines.
     */
    static LinkedHashMap<String, String> parseFasta(File fastaFile) {
        def result = new LinkedHashMap<String, String>()
        String currentId = null
        StringBuilder currentSeq = new StringBuilder()

        def reader = fastaFile.name.endsWith('.gz') ?
            new BufferedReader(new InputStreamReader(new java.util.zip.GZIPInputStream(new FileInputStream(fastaFile)))) :
            fastaFile.newReader()

        try {
            reader.eachLine { line ->
                line = line.trim()
                if (line.startsWith('>')) {
                    if (currentId != null) {
                        result[currentId] = currentSeq.toString()
                    }
                    // Take first whitespace-delimited token as ID
                    currentId = line.substring(1).trim().split(/\s+/)[0]
                    currentSeq = new StringBuilder()
                } else if (!line.startsWith(';') && !line.isEmpty()) {
                    currentSeq.append(line)
                }
            }
            if (currentId != null) {
                result[currentId] = currentSeq.toString()
            }
        } finally {
            reader.close()
        }
        result
    }

    /**
     * Read FASTA from an InputStream.
     */
    static LinkedHashMap<String, String> parseFasta(InputStream input) {
        def tmpFile = File.createTempFile('gspa_fasta_', '.fa')
        try {
            tmpFile.bytes = input.bytes
            return parseFasta(tmpFile)
        } finally {
            tmpFile.delete()
        }
    }
}
