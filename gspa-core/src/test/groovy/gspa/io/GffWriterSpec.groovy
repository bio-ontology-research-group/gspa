package gspa.io

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Feature
import gspa.model.FeatureType
import gspa.model.Genome
import gspa.model.Protein
import gspa.model.Strand
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Files
import java.nio.file.Path

class GffWriterSpec extends Specification {

    @TempDir
    Path tmp

    private static Feature cds(String id, int start, int end, Strand strand = Strand.PLUS) {
        new Feature(
            seqId: 'contig1',
            type: FeatureType.CDS,
            start: start,
            end: end,
            strand: strand,
            attributes: ['ID': [id]],
        )
    }

    private static Genome buildGenome(Feature cds, Protein protein) {
        def contig = new Contig(id: 'contig1', sequence: 'A' * 3000)
        contig.addFeature(cds)
        contig.addProtein(protein)
        def genome = new Genome(id: 'g1', contigs: [contig])
        contig.genome = genome
        genome
    }

    def "whole-protein annotations go on the CDS attribute line, not as subfeatures"() {
        given:
        def feat = cds('cds1', 100, 399)
        def prot = new Protein(id: 'cds1', sourceFeature: feat, sequence: 'M' * 100)
        prot.annotations.add(new Annotation(
            type: AnnotationType.GO, value: 'GO:0006412',
            source: 'diamond', score: 0.9,
        ))
        def out = Files.createTempFile(tmp, 'gff', '.gff').toFile()

        when:
        GffWriter.writeGff3(buildGenome(feat, prot), out)
        List<String> lines = out.text.readLines()

        then:
        lines.any { it.contains('Ontology_term=GO:0006412') }
        lines.count { it.startsWith('contig1') && !it.contains('##') } == 1
    }

    def "region annotation emits a subfeature line with Derives_from the parent CDS"() {
        given:
        def feat = cds('cds1', 100, 399) // 300 nt CDS = 100 residues
        def prot = new Protein(id: 'cds1', sourceFeature: feat, sequence: 'M' * 100)
        prot.annotations.add(new Annotation(
            type: AnnotationType.DISORDER,
            value: 'disorder',
            source: 'metapredict',
            score: 0.82,
            evidenceType: EvidenceType.SEQUENCE_REGION_ML,
            regionStart: 11,      // residue 11..20 → nt 30..59 → genomic 130..159
            regionEnd: 20,
            regionType: 'disorder',
        ))
        def out = Files.createTempFile(tmp, 'gff', '.gff').toFile()

        when:
        GffWriter.writeGff3(buildGenome(feat, prot), out)
        List<String> lines = out.text.readLines()
        def subLine = lines.find { it.contains('polypeptide_region') }

        then:
        subLine != null
        subLine.contains('Derives_from=cds1')
        subLine.contains('region_type=disorder')
        // genomic coordinates: (res-1)*3 offset → 100 + 30 = 130, res 20 end → 100 + 60 - 1 = 159
        def cols = subLine.split('\t')
        cols[0] == 'contig1'
        cols[2] == 'polypeptide_region'
        cols[3] == '130'
        cols[4] == '159'
    }

    def "signal peptide region maps to signal_peptide GFF3 type"() {
        given:
        def feat = cds('cds1', 1, 300)
        def prot = new Protein(id: 'cds1', sourceFeature: feat, sequence: 'M' * 100)
        prot.annotations.add(new Annotation(
            type: AnnotationType.SIGNAL_PEPTIDE,
            value: 'Sec/SPI',
            source: 'signalp',
            score: 0.95,
            evidenceType: EvidenceType.SEQUENCE_REGION_ML,
            regionStart: 1,
            regionEnd: 22,
            regionType: 'signal_peptide',
        ))
        def out = Files.createTempFile(tmp, 'gff', '.gff').toFile()

        when:
        GffWriter.writeGff3(buildGenome(feat, prot), out)
        def subLine = out.text.readLines().find {
            it.contains('signal_peptide') && it.contains('Derives_from=cds1')
        }

        then:
        subLine != null
        def cols = subLine.split('\t')
        cols[2] == 'signal_peptide'
        cols[3] == '1'
        cols[4] == '66'
    }

    def "reverse-strand CDS maps residues to inverted genomic coordinates"() {
        given:
        def feat = cds('cds1', 100, 399, Strand.MINUS)
        def prot = new Protein(id: 'cds1', sourceFeature: feat, sequence: 'M' * 100)
        prot.annotations.add(new Annotation(
            type: AnnotationType.TRANSMEMBRANE,
            value: 'tm1',
            source: 'deeptmhmm',
            score: 0.9,
            evidenceType: EvidenceType.SEQUENCE_REGION_ML,
            regionStart: 11,
            regionEnd: 20,
            regionType: 'tm_helix',
        ))
        def out = Files.createTempFile(tmp, 'gff', '.gff').toFile()

        when:
        GffWriter.writeGff3(buildGenome(feat, prot), out)
        def subLine = out.text.readLines().find {
            it.contains('transmembrane_region') && it.contains('Derives_from=cds1')
        }

        then:
        subLine != null
        def cols = subLine.split('\t')
        // residue 11..20 → nt 30..59 offset; reverse: end - nt_end .. end - nt_start
        // = 399 - 59 = 340, 399 - 30 = 369
        cols[3] == '340'
        cols[4] == '369'
        cols[6] == '-'
    }
}
