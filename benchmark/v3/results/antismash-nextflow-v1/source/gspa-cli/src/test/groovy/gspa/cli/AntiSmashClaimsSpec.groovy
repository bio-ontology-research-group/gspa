package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimExtractor
import gspa.predictor.specialized.AntiSmashPredictor
import picocli.CommandLine
import spock.lang.Specification
import spock.lang.TempDir

class AntiSmashClaimsSpec extends Specification {
    @TempDir File dir
    File fasta, proteins, gff, nativeJson, output

    void setup() {
        fasta = new File(dir, 'genome.fna'); fasta.text = '>chr\n' + 'ATG' * 100 + '\n'
        proteins = new File(dir, 'proteins.faa')
        // Native FASTA IDs deliberately differ from GFF IDs, and file orders differ.
        proteins.text = '>chr_2 # 61 # 90 # -1 # ID=1_2;partial=00\nMMMMMMMMMM\n' +
            '>chr_1 # 1 # 30 # 1 # ID=1_1;partial=00\nMMMMMMMMMM\n'
        gff = new File(dir, 'genes.gff')
        gff.text = '##gff-version 3\nchr\tpyrodigal\tCDS\t1\t30\t.\t+\t0\tID=1_1\n' +
            'chr\tpyrodigal\tCDS\t61\t90\t.\t-\t0\tID=1_2\n'
        nativeJson = new File(dir, 'result.json')
        nativeJson.text = new ObjectMapper().writeValueAsString([schema: 4, version: '8.0.4', records: [[
            id: 'chr', seq: [data: 'ATG' * 100], annotations: [:],
            modules: ['antismash.detection.hmm_detection': [record_id: 'chr', rule_results: [:], enabled_types: []]],
            areas: [[start: 29, end: 60, products: ['T3PKS']], [start: 100, end: 120, products: ['terpene']]]]]])
        output = new File(dir, 'claims.jsonl')
    }

    int execute() {
        new CommandLine(new GspaMain()).setErr(new PrintWriter(new StringWriter())).execute(
            'antismash-claims', '--results', nativeJson.path, '--genome', fasta.path,
            '--gff3', gff.path, '--proteins', proteins.path, '--output', output.path)
    }

    def 'native caller IDs join by identity and coordinates with exact endpoint behavior'() {
        expect:
        execute() == 0
        def claims = new ClaimExtractor().readClaimsJsonl(output)
        claims.size() == 1
        claims[0].proteinId == 'chr_1'
        claims[0].metadata.cdsOverlap == 'partial'
        !claims[0].hasFunctionalScore()
        claims[0].calibratedProb == null
        def receipt = new ObjectMapper().readValue(new File(output.path + '.import.json'), Map)
        receipt.coverage == 'complete'
        receipt.region_count == 2 // The intergenic region must survive even without a protein overlap.
        receipt.regions.size() == 2
        receipt.records.size() == 1
        receipt.protein_count == 2
        receipt.annotation_count == 1
        receipt.action == 'import_existing_output'
        [receipt.results_sha256, receipt.genome_sha256, receipt.gff3_sha256,
         receipt.proteins_sha256, receipt.output_sha256].every { it ==~ /[0-9a-f]{64}/ }
    }

    def 'explicit protein IDs preserve exactly the shared parser and projection claims'() {
        given:
        proteins.text = '>p1\nMMMMMMMMMM\n>p2\nMMMMMMMMMM\n'
        gff.text = gff.text.replace('ID=1_1', 'ID=p1').replace('ID=1_2', 'ID=p2')
        def genome = AntiSmashClaimsCommand.loadCoordinates(fasta, gff, proteins)
        def projected = new AntiSmashPredictor().importGenomeResults(genome, nativeJson)
        def expected = new File(dir, 'expected.jsonl')
        def extractor = new ClaimExtractor()
        ClaimExtractor.writeClaimsJsonl(projected.projection.byProtein.collectMany { id, calls ->
            calls.collect { extractor.fromAnnotation(id, it) }
        }, expected)
        expect:
        execute() == 0
        output.bytes == expected.bytes
    }

    def 'invalid sequence identity and CDS joins fail without publishing partial output'() {
        given:
        switch (fault) {
            case 'empty genome': fasta.text = ''; break
            case 'duplicate protein': proteins << proteins.text; break
            case 'duplicate contig': fasta << fasta.text; break
            case 'wrong contig': fasta.text = fasta.text.replace('chr', 'other'); break
            case 'wrong sequence': fasta.text = fasta.text.replace('ATG', 'AAA'); break
            case 'wrong CDS header': proteins.text = proteins.text.replace('# 61 #', '# 62 #'); break
            case 'wrong CDS identity': proteins.text = proteins.text.replace('ID=1_1', 'ID=unknown'); break
            case 'wrong strand': proteins.text = proteins.text.replace('# -1 #', '# 1 #'); break
            case 'missing protein': proteins.text = proteins.text.split('>chr_1')[0]; break
            case 'out of bounds': gff.text = gff.text.replace('\t90\t', '\t301\t'); break
            case 'ambiguous CDS': gff << gff.text; break
            case 'malformed GFF': gff << 'broken\n'; break
            case 'embedded FASTA': gff << '##FASTA\n'; break
            case 'invalid JSON': nativeJson.text = '[]'; break
        }
        expect:
        execute() != 0
        !output.exists()
        !new File(output.path + '.import.json').exists()
        where:
        fault << ['empty genome', 'duplicate protein', 'duplicate contig', 'wrong contig', 'wrong sequence',
            'wrong CDS header', 'wrong CDS identity', 'wrong strand', 'missing protein',
            'out of bounds', 'ambiguous CDS', 'malformed GFF', 'embedded FASTA', 'invalid JSON']
    }

    def 'a missing result record is partial coverage rather than evidence of absence'() {
        given:
        fasta << '>unassessed\nATGATGATG\n'
        expect:
        execute() == 0
        def receipt = new ObjectMapper().readValue(new File(output.path + '.import.json'), Map)
        receipt.coverage == 'partial'
        receipt.unassessed_inputs == ['unassessed']
    }

    def 'complete zero-region output is preserved as zero observations'() {
        given:
        def data = new ObjectMapper().readValue(nativeJson, Map)
        data.records[0].areas = []
        nativeJson.text = new ObjectMapper().writeValueAsString(data)
        expect:
        execute() == 0
        output.length() == 0
        def receipt = new ObjectMapper().readValue(new File(output.path + '.import.json'), Map)
        receipt.coverage == 'complete'
        receipt.region_count == 0
        receipt.protein_count == 2
    }

    def 'existing outputs are never overwritten'() {
        given:
        output.text = 'existing\n'
        expect:
        execute() != 0
        output.text == 'existing\n'
    }
}
