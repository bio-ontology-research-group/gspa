package gspa.io

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Feature
import gspa.model.FeatureType
import gspa.model.Genome
import gspa.model.Protein
import gspa.model.Strand

/**
 * Writes GFF3 format genome annotation files.
 */
class GffWriter {

    /**
     * Write a genome's features to GFF3 format.
     * Optionally includes FASTA sequences after ##FASTA directive.
     * When includeAnnotations is true, adds GO/EC/Pfam to CDS attributes
     * and emits region-level annotations (signal peptide cleavage,
     * transmembrane helices, intrinsic disorder) as GFF3 subfeatures
     * whose {@code Derives_from} points back to the parent CDS.
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
                    Protein protein = null
                    if (includeAnnotations && feature.type == FeatureType.CDS) {
                        protein = contig.proteins.find {
                            it.sourceFeature == feature || it.id == feature.id
                        }
                        if (protein && !protein.annotations.isEmpty()) {
                            addAnnotationAttributes(feature, protein)
                        }
                    }
                    writer.writeLine(feature.toGff3())
                    if (protein) {
                        writeRegionSubfeatures(writer, feature, protein)
                    }
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
     * Emit region-level annotations (those with {@link Annotation#hasRegion()})
     * as GFF3 subfeature lines derived from the parent CDS. Residue indices
     * are mapped to genomic coordinates via the CDS start / strand.
     */
    private static void writeRegionSubfeatures(Writer writer, Feature cds, Protein protein) {
        String parentId = cds.id
        if (!parentId) return
        int subIdx = 0
        protein.annotations.annotations.each { Annotation ann ->
            if (!ann.hasRegion()) return
            int[] coords = residueToGenomic(cds, ann.regionStart, ann.regionEnd)
            if (coords == null) return
            String gff3Type = regionTypeToGff3(ann)
            subIdx++
            String id = "${parentId}_${ann.regionType ?: 'region'}_${subIdx}"
            List<String> attrs = [
                "ID=${id}".toString(),
                "Derives_from=${parentId}".toString(),
            ]
            if (ann.regionType) {
                attrs << "region_type=${ann.regionType}".toString()
            }
            if (ann.value) {
                attrs << "Name=${ann.value}".toString()
            }
            if (ann.source) {
                attrs << "source_tool=${ann.source}".toString()
            }
            String scoreStr = ann.score > 0 ? String.format(Locale.ROOT, '%.3f', ann.score) : '.'
            String strandStr = cds.strand.symbol
            writer.writeLine(
                "${cds.seqId}\t${ann.source ?: 'gspa'}\t${gff3Type}\t${coords[0]}\t${coords[1]}\t${scoreStr}\t${strandStr}\t.\t${attrs.join(';')}"
            )
        }
    }

    /**
     * Map 1-based residue positions within a protein to 1-based genomic
     * positions on the CDS contig. Returns {@code null} if the CDS is too
     * short for the requested residue range.
     */
    private static int[] residueToGenomic(Feature cds, int resStart, int resEnd) {
        int ntStart = (resStart - 1) * 3
        int ntEnd = resEnd * 3 - 1
        int cdsLen = cds.end - cds.start + 1
        if (ntEnd >= cdsLen) return null
        if (cds.strand == Strand.MINUS) {
            int gStart = cds.end - ntEnd
            int gEnd = cds.end - ntStart
            return [gStart, gEnd] as int[]
        }
        int gStart = cds.start + ntStart
        int gEnd = cds.start + ntEnd
        [gStart, gEnd] as int[]
    }

    private static String regionTypeToGff3(Annotation ann) {
        if (ann.type == AnnotationType.SIGNAL_PEPTIDE) return 'signal_peptide'
        if (ann.type == AnnotationType.TRANSMEMBRANE) return 'transmembrane_region'
        if (ann.regionType == 'signal_peptide') return 'signal_peptide'
        if (ann.regionType == 'tm_helix') return 'transmembrane_region'
        'polypeptide_region'
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
