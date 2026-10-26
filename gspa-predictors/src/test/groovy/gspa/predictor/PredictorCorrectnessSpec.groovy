package gspa.predictor

import gspa.model.*
import gspa.predictor.context.OperonEnsemble
import gspa.predictor.context.OperonPredictor
import gspa.predictor.domain.HmmerPredictor
import gspa.predictor.domain.InterProScanPredictor
import gspa.predictor.function.EggNogMapperPredictor
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * Regression tests for v3 predictor correctness fixes, using the real output
 * formats of emapper 2.1.x, InterProScan 5 (TSV) and hmmsearch.
 */
class PredictorCorrectnessSpec extends Specification {

    @TempDir
    Path tmp

    // --- eggNOG-mapper: E-value vs bit score ---

    private static final String EMAPPER_HEADER =
        '#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\t' +
        'Preferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\t' +
        'BRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs'

    private File emapperFile(String body) {
        def f = tmp.resolve('eggnog_results.emapper.annotations').toFile()
        f.text = body
        f
    }

    private static String row(String id, String evalue, String score, String go = 'GO:0004672') {
        [id, '562.P0A7S3', evalue, score, 'COG0515@1|root', 'Bacteria', 'T', 'kinase', '-', go, '2.7.11.1',
         'ko:K08884', '-', '-', '-', '-', '-', '-', '-', '-', 'Pkinase'].join('\t')
    }

    def "eggNOG: confidence comes from the evalue column, not the bit score"() {
        given: "bit score 480 would give 0.089 if misread as an E-value; E-value 1e-15 gives 0.5"
        emapperFile("## emapper-2.1.12\n${EMAPPER_HEADER}\n${row('p1', '1e-15', '480')}\n" +
            "${row('p2', '0.0', '900')}\n## 2 queries scanned\n")

        when:
        def res = new EggNogMapperPredictor().parseOutput(tmp.toFile())

        then:
        res['p1'].find { it.type == AnnotationType.GO }.score == 0.5d
        res['p1'].find { it.type == AnnotationType.EC }.value == 'EC:2.7.11.1'
        and: "an E-value of exactly 0 is the strongest hit, not the 0.7 'unknown' default"
        res['p2'].every { it.score == 1.0d }
    }

    def "eggNOG: columns are located by header, even when the layout is shifted"() {
        given: "a header whose columns are in a different order than emapper 2.1"
        def shifted = '#query\tevalue\tscore\tGOs\tEC\tKEGG_ko\teggNOG_OGs\tCOG_category\n'
        emapperFile(shifted + 'p1\t1e-30\t300\tGO:0005524\t3.6.4.12\tko:K03657\tCOG0210@1|root\tL\n')

        when:
        def res = new EggNogMapperPredictor().parseOutput(tmp.toFile())

        then:
        res['p1'].find { it.type == AnnotationType.GO }.value == 'GO:0005524'
        res['p1'].find { it.type == AnnotationType.GO }.score == 1.0d
        res['p1'].find { it.type == AnnotationType.EC }.value == 'EC:3.6.4.12'
        res['p1'].find { it.type == AnnotationType.COG }.value == 'L'
    }

    def "eggNOG: headerless file uses the emapper 2.1 layout and still validates the evalue"() {
        given:
        emapperFile(row('p1', '1e-15', '480') + '\n')

        when:
        def res = new EggNogMapperPredictor().parseOutput(tmp.toFile())

        then:
        res['p1'].find { it.type == AnnotationType.GO }.score == 0.5d

        when: "a shifted layout puts text in the assumed evalue column"
        emapperFile(row('p1', '1e-15', '480') + '\n' + row('bad', 'COG0515@1|root', 'x') + '\n')
        new EggNogMapperPredictor().parseOutput(tmp.toFile())

        then:
        def e = thrown(ToolOutputContractException)
        e.message.contains('bad')
        e.message.contains('emapper version')
    }

    def "eggNOG: invalid E-values are a failed contract, never an invented confidence (#why)"() {
        given:
        emapperFile("${EMAPPER_HEADER}\n${row('p1', evalue, '480')}\n")

        when:
        new EggNogMapperPredictor().parseOutput(tmp.toFile())

        then:
        def e = thrown(ToolOutputContractException)
        e.message.contains('p1')
        e.message.contains('line 2')

        where:
        why        | evalue
        'negative' | '-1e-5'
        'NaN'      | 'NaN'
        'infinite' | 'Infinity'
        'dash'     | '-'
        'empty'    | ''
        'text'     | 'COG0515'
    }

    def "eggNOG: E-value 0 is valid; large E-values are weak (0) not errors; short rows are rejected"() {
        given:
        emapperFile("${EMAPPER_HEADER}\n${row('p1', '0.0', '900')}\n${row('p2', '5.0', '30')}\n")

        when:
        def res = new EggNogMapperPredictor().parseOutput(tmp.toFile())

        then:
        res['p1'].every { it.score == 1.0d }
        res['p2'].every { it.score == 0.0d }

        when:
        emapperFile("${EMAPPER_HEADER}\np3\tseed\t1e-9\n")
        new EggNogMapperPredictor().parseOutput(tmp.toFile())

        then:
        thrown(ToolOutputContractException)
    }

    def "eggNOG: a header missing the evalue column is rejected instead of guessed"() {
        given:
        emapperFile('#query\tseed_ortholog\tGOs\tEC\tKEGG_ko\n')

        when:
        new EggNogMapperPredictor().parseOutput(tmp.toFile())

        then:
        thrown(IllegalStateException)
    }

    def "eggNOG: missing annotations file after a run is an error, header-only file is zero"() {
        when:
        new EggNogMapperPredictor().parseOutput(tmp.toFile())

        then:
        thrown(IllegalStateException)

        when:
        emapperFile(EMAPPER_HEADER + '\n## 0 queries scanned\n')

        then:
        new EggNogMapperPredictor().parseOutput(tmp.toFile()).isEmpty()
    }

    // --- InterProScan: no EC in the TSV schema ---

    def "InterProScan: does not advertise or emit EC; '-' placeholders and GO sources are handled"() {
        given: "15-column TSV as written with -goterms -pa; the KEGG pathway cell contains an EC-like token"
        def ips = tmp.resolve('interproscan_results.tsv').toFile()
        ips.text = [
            ['p1', 'md5', '300', 'Pfam', 'PF00069', 'Pkinase', '10', '280', '1.2E-45', 'T', '01-01-2025',
             'IPR000719', 'Protein kinase domain', 'GO:0004672(InterPro)|GO:0006468(InterPro)',
             'KEGG: 00010+2.7.1.1'],
            ['p1', 'md5', '300', 'SMART', 'SM00220', 'S_TKc', '10', '280', '2.0E-40', 'T', '01-01-2025',
             '-', '-', '-', '-'],
        ]*.join('\t').join('\n') + '\n'

        when:
        def predictor = new InterProScanPredictor()
        def res = predictor.parseOutput(tmp.toFile())

        then:
        !(AnnotationType.EC in predictor.outputTypes)
        res['p1'].every { it.type != AnnotationType.EC }
        res['p1'].findAll { it.type == AnnotationType.GO }*.value.sort() == ['GO:0004672', 'GO:0006468']
        res['p1'].findAll { it.type == AnnotationType.INTERPRO }*.value == ['IPR000719']
    }

    // --- HMMER database ---

    def "HMMER: without a database file the predictor is unavailable and says why"() {
        expect:
        !new HmmerPredictor().isAvailable()
        new HmmerPredictor().unavailableReason().contains('database')
        !new HmmerPredictor(hmmDatabase: tmp.resolve('missing.hmm').toString()).isAvailable()
    }

    def "HMMER: command passes the configured database"() {
        given:
        def hmm = tmp.resolve('Pfam-A.hmm').toFile(); hmm.text = 'HMMER3/f\n'
        def faa = tmp.resolve('in.faa').toFile(); faa.text = '>p\nMKV\n'

        when:
        def cmd = new HmmerPredictor(hmmDatabase: hmm.absolutePath).buildCommand(faa, tmp.toFile())

        then:
        cmd[-2] == hmm.absolutePath
        cmd[-1] == faa.absolutePath

        when:
        new HmmerPredictor().buildCommand(faa, tmp.toFile())

        then:
        thrown(IllegalStateException)
    }

    // --- Operon: aspect-less GO terms are not biological process ---

    private Genome operonGenome() {
        def genome = new Genome(id: 'g')
        def contig = new Contig(id: 'c1', sequence: 'A' * 2000)
        [[id: 'p1', s: 1, e: 300], [id: 'p2', s: 310, e: 600]].each { spec ->
            def f = new Feature(seqId: 'c1', type: FeatureType.CDS, start: spec.s, end: spec.e, strand: Strand.PLUS)
            f.id = spec.id
            contig.addFeature(f)
            contig.addProtein(new Protein(id: spec.id, sequence: 'M' * 99, sourceFeature: f))
        }
        genome.addContig(contig)
        genome
    }

    def "Operon: aspect-less cytoplasm / enzyme GO terms are never copied to neighbours"() {
        given:
        def genome = operonGenome()
        def p1 = genome.findProtein('p1')
        // GO:0005737 cytoplasm (CC), GO:0016301 kinase activity (MF), both without aspect;
        // GO:0006468 protein phosphorylation (BP) with an explicit aspect.
        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0005737', source: 'eggnog-mapper'))
        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0016301', source: 'interproscan'))
        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0006468', source: 'diamond', goAspect: 'BP'))

        when:
        def transfers = new OperonPredictor().predictGenome(genome)

        then:
        transfers['p2']*.value == ['GO:0006468']
    }

    def "Operon ensemble: aspect-less shared GO terms do not support a co-operonic call"() {
        given:
        def genome = operonGenome()
        ['p1', 'p2'].each {
            genome.findProtein(it).annotations.add(
                new Annotation(type: AnnotationType.GO, value: 'GO:0005737', source: 'eggnog-mapper'))
        }

        expect: "no BP-labelled term in common, so no shared-function vote"
        new OperonEnsemble().detect(genome).every { !it.supportSet.contains('functional') }
    }
}
