package gspa.io

import gspa.model.*

/**
 * Writes GenBank flat file format for genome annotations.
 * Produces files suitable for NCBI submission.
 */
class GenbankWriter {

    static void writeGenbank(Genome genome, File output) {
        output.withWriter { writer ->
            genome.contigs.each { contig ->
                writeContigRecord(writer, contig, genome)
            }
        }
    }

    private static void writeContigRecord(Writer writer, Contig contig, Genome genome) {
        int seqLen = contig.length

        // LOCUS line
        String molecule = 'DNA'
        String topology = 'linear'
        String division = genome.taxonomy?.domain?.isProkaryote() ? 'BCT' : 'UNK'
        writer.writeLine(String.format("LOCUS       %-16s %11d bp    %s     %s %s",
            contig.id.take(16), seqLen, molecule, topology, division))

        // DEFINITION
        String definition = genome.taxonomy?.scientificName ?: genome.id
        writer.writeLine("DEFINITION  ${definition}, ${contig.id}.")

        // ACCESSION
        writer.writeLine("ACCESSION   ${contig.id}")

        // VERSION
        writer.writeLine("VERSION     ${contig.id}")

        // SOURCE / ORGANISM
        if (genome.taxonomy) {
            writer.writeLine("SOURCE      ${genome.taxonomy.scientificName ?: 'Unknown'}")
            writer.writeLine("  ORGANISM  ${genome.taxonomy.scientificName ?: 'Unknown'}")
        }

        // FEATURES
        writer.writeLine("FEATURES             Location/Qualifiers")
        writer.writeLine("     source          1..${seqLen}")
        if (genome.taxonomy?.scientificName) {
            writer.writeLine("                     /organism=\"${genome.taxonomy.scientificName}\"")
        }

        contig.features.each { feature ->
            // Find the protein associated with this feature (if CDS)
            Protein protein = null
            if (feature.type == FeatureType.CDS) {
                protein = contig.proteins.find { it.sourceFeature == feature || it.id == feature.id }
            }
            writeFeature(writer, feature, protein)
        }

        // ORIGIN (sequence)
        if (contig.sequence) {
            writer.writeLine("ORIGIN")
            writeSequence(writer, contig.sequence)
        }

        writer.writeLine("//")
    }

    /**
     * Write a feature with optional functional annotations from the associated protein.
     */
    static void writeFeature(Writer writer, Feature feature, Protein protein = null) {
        String location = formatLocation(feature)
        String type = feature.type == FeatureType.CDS ? 'CDS' : feature.type.gff3Type
        writer.writeLine(String.format("     %-15s %s", type, location))

        if (feature.locusTag) {
            writer.writeLine("                     /locus_tag=\"${feature.locusTag}\"")
        }
        if (feature.product) {
            writer.writeLine("                     /product=\"${feature.product}\"")
        }
        if (feature.type == FeatureType.CDS) {
            writer.writeLine("                     /transl_table=11")
        }

        // Write functional annotations from the protein
        if (protein && !protein.annotations.isEmpty()) {
            // EC numbers
            protein.annotations.byType(AnnotationType.EC).collect { it.value }.unique().each { ec ->
                writer.writeLine("                     /EC_number=\"${ec.replaceFirst(/^EC:/, '')}\"")
            }
            // GO terms as /note qualifiers (standard GenBank practice)
            def goTerms = protein.annotations.byType(AnnotationType.GO).collect { it.value }.unique()
            if (goTerms) {
                writer.writeLine("                     /note=\"GO: ${goTerms.join(', ')}\"")
            }
            // Pfam domains as /note
            def pfamDomains = protein.annotations.byType(AnnotationType.PFAM).collect { it.value }.unique()
            if (pfamDomains) {
                writer.writeLine("                     /note=\"Pfam: ${pfamDomains.join(', ')}\"")
            }
            // InterPro
            protein.annotations.byType(AnnotationType.INTERPRO).collect { it.value }.unique().each { ipr ->
                writer.writeLine("                     /db_xref=\"InterPro:${ipr}\"")
            }
            // Protein sequence
            if (protein.sequence) {
                writer.writeLine("                     /translation=\"${protein.sequence}\"")
            }
        }
    }

    private static String formatLocation(Feature feature) {
        if (feature.strand == Strand.MINUS) {
            return "complement(${feature.start}..${feature.end})"
        }
        return "${feature.start}..${feature.end}"
    }

    private static void writeSequence(Writer writer, String sequence) {
        String lower = sequence.toLowerCase()
        for (int i = 0; i < lower.length(); i += 60) {
            int end = Math.min(i + 60, lower.length())
            StringBuilder line = new StringBuilder()
            line.append(String.format("%9d", i + 1))
            for (int j = i; j < end; j += 10) {
                line.append(' ')
                line.append(lower.substring(j, Math.min(j + 10, end)))
            }
            writer.writeLine(line.toString())
        }
    }
}
