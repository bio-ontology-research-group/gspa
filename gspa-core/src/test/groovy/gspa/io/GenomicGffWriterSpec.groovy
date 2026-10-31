package gspa.io

import groovy.json.JsonSlurper
import gspa.integration.EvidenceType
import gspa.model.*
import spock.lang.Specification
import spock.lang.TempDir

class GenomicGffWriterSpec extends Specification {
    @TempDir File tmp

    private static Annotation region() {
        new Annotation(type: AnnotationType.VIRAL_CONTIG, value: 'viral; candidate=one',
            source: 'genomad', score: 0.0, contigId: 'ctg', genomicStart: 2,
            genomicEnd: 90, regionType: 'virus', evidence: 'observed,region',
            evidenceType: EvidenceType.GENOMIC_REGION_ML,
            metadata: [description: 'a;b=c,d% +\t\n', count: 0],
            provenance: ['tool:genomad', 'raw score=0'])
    }

    private static Genome genome(Annotation ann) {
        new Genome(id: 'test', contigs: [new Contig(id: 'ctg', sequence: 'A' * 100)],
            metadata: ['genomic_regions:genomad': [ann]])
    }

    private static List<List<String>> rows(File file) {
        file.readLines().findAll { !it.startsWith('#') && it.contains('\t') }
            .collect { it.split('\t', -1).toList() }
    }

    private static Map attributes(List<String> row) {
        row[8].split(';').collectEntries {
            def pair = it.split('=', 2)
            [(pair[0]): URLDecoder.decode(pair[1], 'UTF-8')]
        }
    }

    def 'full genomic region survives without CDS and preserves evidence including zero score'() {
        given:
        def ann = region()
        def g = genome(ann)
        def output = new File(tmp, 'no-cds.gff3')

        when:
        GffWriter.writeGff3(g, output)
        def records = rows(output)

        then:
        records.size() == 1
        records[0].take(8) == ['ctg', 'genomad', 'region', '2', '90', '0.0', '.', '.']
        def attrs = attributes(records[0])
        attrs.Name == ann.value
        attrs.region_type == 'virus'
        attrs.annotation_type == 'VIRAL_CONTIG'
        attrs.evidence == ann.evidence
        attrs.evidence_type == 'GENOMIC_REGION_ML'
        new JsonSlurper().parseText(attrs.metadata) == ann.metadata
        new JsonSlurper().parseText(attrs.provenance) == ann.provenance
        !attrs.containsKey('Parent')
        !attrs.containsKey('Derives_from')
        g.proteins.empty
    }

    def 'unscored genomic observations use a missing GFF score while genuine zero remains numeric'() {
        given:
        def ann = region()
        ann.metadata.score_available = false
        def output = new File(tmp, 'unscored.gff3')
        when:
        GffWriter.writeGff3(genome(ann), output)
        then:
        rows(output)[0][5] == '.'
        new JsonSlurper().parseText(attributes(rows(output)[0]).metadata).score_available == false
    }

    def 'overlapping CDS copies do not duplicate regions or collide with feature and residue IDs'() {
        given:
        def ann = region()
        def g = genome(ann)
        ['gspa_genomic_region_1', 'other'].each { id ->
            def cds = new Feature(seqId: 'ctg', type: FeatureType.CDS, start: 1,
                end: 90, strand: Strand.PLUS, attributes: [ID: [id]])
            def p = new Protein(id: id, sequence: 'M' * 30, sourceFeature: cds)
            p.annotations.add(ann)
            g.contigs[0].addFeature(cds)
            g.contigs[0].addProtein(p)
        }
        // This residue subfeature's existing naming convention occupies ID 2.
        def cds = new Feature(seqId: 'ctg', type: FeatureType.CDS, start: 1,
            end: 90, strand: Strand.PLUS, attributes: [ID: ['gspa_genomic']])
        def p = new Protein(id: 'gspa_genomic', sequence: 'M' * 30, sourceFeature: cds)
        [1, 2].each { p.annotations.add(new Annotation(type: AnnotationType.DISORDER,
            value: "part${it}", source: 'fixture', regionStart: it, regionEnd: it,
            regionType: it == 1 ? 'disorder' : 'region')) }
        g.contigs[0].addFeature(cds)
        g.contigs[0].addProtein(p)
        def output = new File(tmp, 'overlap.gff3')

        when:
        GffWriter.writeGff3(g, output)
        def records = rows(output)
        def originals = records.findAll { it[2] == 'region' }

        then:
        originals.size() == 1
        attributes(originals[0]).ID == 'gspa_genomic_region_3'
        originals[0][3..4] == ['2', '90']
        !attributes(originals[0]).containsKey('Derives_from')
        // The new genomic row must avoid every occupied ID.
        records.count { attributes(it).ID == 'gspa_genomic_region_3' } == 1
    }

    def 'annotation opt-out matches ordinary output byte for byte'() {
        given:
        def g = genome(region())
        def plain = new File(tmp, 'plain.gff3')
        def suppressed = new File(tmp, 'suppressed.gff3')

        when:
        GffWriter.writeGff3(g, suppressed, false, false)
        g.metadata.clear()
        GffWriter.writeGff3(g, plain)

        then:
        suppressed.bytes == plain.bytes
        plain.text == '##gff-version 3\n##sequence-region ctg 1 100\n'
    }

    def 'colocated CDS records retain their own protein annotations'() {
        given:
        def g = genome(region())
        g.metadata.clear()
        ['one': 'GO:0003674', 'two': 'GO:0008150'].each { id, term ->
            def cds = new Feature(seqId: 'ctg', type: FeatureType.CDS, start: 1,
                end: 90, strand: Strand.PLUS, attributes: [ID: [id]])
            def p = new Protein(id: id, sequence: 'M' * 30, sourceFeature: cds)
            p.annotations.add(new Annotation(type: AnnotationType.GO, value: term, source: 'fixture'))
            g.contigs[0].addFeature(cds)
            g.contigs[0].addProtein(p)
        }
        def output = new File(tmp, 'colocated.gff3')

        when:
        GffWriter.writeGff3(g, output)
        def byId = rows(output).collectEntries { row ->
            def attrs = attributes(row)
            [(attrs.ID): attrs]
        }

        then:
        byId.one.protein_id == 'one'
        byId.one.Ontology_term == 'GO:0003674'
        byId.two.protein_id == 'two'
        byId.two.Ontology_term == 'GO:0008150'
    }

    def 'genomic rows precede the optional FASTA section'() {
        given:
        def output = new File(tmp, 'with-fasta.gff3')

        when:
        GffWriter.writeGff3(genome(region()), output, true)
        def lines = output.readLines()

        then:
        rows(output).size() == 1
        lines.findIndexOf { it.contains('\tregion\t') } < lines.indexOf('##FASTA')
        lines[lines.indexOf('##FASTA') + 1] == '>ctg'
    }
}
