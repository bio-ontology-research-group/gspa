package gspa.io

import gspa.model.AnnotationType
import gspa.model.FeatureType
import gspa.model.Genome
import gspa.model.Protein

/**
 * Writes GFF3 format genome annotation files.
 */
class GffWriter {

    /**
     * Write a genome's features to GFF3 format.
     * Optionally includes FASTA sequences after ##FASTA directive.
     * When includeAnnotations is true, adds GO/EC/Pfam to CDS attributes.
     */
    static void writeGff3(Genome genome, File output,
                          boolean includeFasta = false,
                          boolean includeAnnotations = true) {
        output.withWriter { writer ->
            writer.writeLine('##gff-version 3')

            // Sequence-region pragmas
            genome.contigs.each { contig ->
                if (contig.length > 0) {
                    writer.writeLine("##sequence-region ${contig.id} 1 ${contig.length}")
                }
            }

            // Features
            genome.contigs.each { contig ->
                contig.features.each { feature ->
                    if (includeAnnotations && feature.type == FeatureType.CDS) {
                        Protein protein = contig.proteins.find {
                            it.sourceFeature == feature || it.id == feature.id
                        }
                        if (protein && !protein.annotations.isEmpty()) {
                            addAnnotationAttributes(feature, protein)
                        }
                    }
                    writer.writeLine(feature.toGff3())
                }
            }

            // Optional FASTA section
            if (includeFasta) {
                writer.writeLine('##FASTA')
                genome.contigs.each { contig ->
                    if (contig.sequence) {
                        writer.writeLine(">${contig.id}")
                        for (int i = 0; i < contig.sequence.length(); i += 70) {
                            writer.writeLine(contig.sequence.substring(i, Math.min(i + 70, contig.sequence.length())))
                        }
                    }
                }
            }
        }
    }

    /**
     * Add functional annotations to a feature's GFF3 attributes.
     */
    private static void addAnnotationAttributes(gspa.model.Feature feature, Protein protein) {
        def goTerms = protein.annotations.byType(AnnotationType.GO).collect { it.value }.unique()
        if (goTerms) {
            feature.attributes['Ontology_term'] = goTerms
        }

        def ecNumbers = protein.annotations.byType(AnnotationType.EC).collect { it.value }.unique()
        if (ecNumbers) {
            feature.attributes['ec_number'] = ecNumbers
        }

        def pfam = protein.annotations.byType(AnnotationType.PFAM).collect { it.value }.unique()
        if (pfam) {
            feature.attributes['Dbxref'] = (feature.attributes['Dbxref'] ?: []) +
                pfam.collect { "Pfam:${it}" }
        }

        def interpro = protein.annotations.byType(AnnotationType.INTERPRO).collect { it.value }.unique()
        if (interpro) {
            feature.attributes['Dbxref'] = (feature.attributes['Dbxref'] ?: []) +
                interpro.collect { "InterPro:${it}" }
        }

        def cog = protein.annotations.byType(AnnotationType.COG).collect { it.value }.unique()
        if (cog) {
            feature.attributes['cog'] = cog
        }

        def kegg = protein.annotations.byType(AnnotationType.KEGG).collect { it.value }.unique()
        if (kegg) {
            feature.attributes['kegg'] = kegg
        }
    }
}
