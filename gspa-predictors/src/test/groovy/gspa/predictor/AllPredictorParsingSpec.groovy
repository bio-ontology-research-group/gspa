package gspa.predictor

import gspa.model.AnnotationType
import gspa.predictor.similarity.MMseqs2Predictor
import gspa.predictor.structure.FoldSeekPredictor
import gspa.predictor.domain.HmmerPredictor
import gspa.predictor.function.EggNogMapperPredictor
import gspa.predictor.pathway.GapseqPredictor
import gspa.predictor.localization.SignalPPredictor
import gspa.predictor.localization.DeepTmhmmPredictor
import gspa.predictor.specialized.AntiSmashPredictor
import gspa.predictor.specialized.CrisprPredictor
import spock.lang.Specification

/**
 * Comprehensive parsing tests for all predictor output formats.
 * Uses realistic mock output files based on actual tool documentation.
 */
class AllPredictorParsingSpec extends Specification {

    // --- MMseqs2 ---

    def "MMseqs2: should parse convertalis output"() {
        given:
        def predictor = new MMseqs2Predictor()
        def outputDir = mockDir('mmseqs2')

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.size() == 3
        // protein_001: GO terms from ribosomal protein
        results['protein_001'].any { it.type == AnnotationType.GO && it.value == 'GO:0003735' }
        results['protein_001'].any { it.type == AnnotationType.GO && it.value == 'GO:0006412' }
        results['protein_001'].any { it.type == AnnotationType.GO && it.value == 'GO:0005840' }
        // protein_002: EC number
        results['protein_002'].any { it.type == AnnotationType.EC && it.value == 'EC:3.1.26.3' }
        // protein_003: EC from prothrombin
        results['protein_003'].any { it.type == AnnotationType.EC && it.value == 'EC:3.4.21.5' }
        // All from mmseqs2
        results.values().flatten().every { it.source == 'mmseqs2' }
    }

    // --- FoldSeek ---

    def "FoldSeek: should parse easy-search output with TM-scores"() {
        given:
        def predictor = new FoldSeekPredictor()
        def outputDir = mockDir('foldseek')

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.size() == 2  // protein_003 has no GO/EC in header
        results['protein_001'].any { it.type == AnnotationType.GO && it.value == 'GO:0003735' }
        // Score should be based on TM-score
        results['protein_001'][0].score == 0.92
        results['protein_002'][0].score == 0.65
        // protein_002 has EC
        results['protein_002'].any { it.type == AnnotationType.EC && it.value == 'EC:3.1.26.3' }
        // Metadata should include TM-score
        results['protein_001'][0].metadata.tmscore == 0.92
    }

    // --- HMMER ---

    def "HMMER: should parse domtblout format"() {
        given:
        def predictor = new HmmerPredictor()
        def outputDir = mockDir('hmmer')

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.size() == 3
        // protein_001: Ribosomal_S12
        results['protein_001'].size() == 1
        results['protein_001'][0].type == AnnotationType.PFAM
        results['protein_001'][0].value == 'PF00177'  // version stripped
        // protein_002: two domains
        results['protein_002'].size() == 2
        results['protein_002'].any { it.value == 'PF03368' }  // Dicer_dimer
        results['protein_002'].any { it.value == 'PF00636' }  // RIBOc
        // protein_003: Trypsin
        results['protein_003'][0].value == 'PF00089'
    }

    // --- eggNOG-mapper ---

    def "eggNOG: should parse emapper.annotations output"() {
        given:
        def predictor = new EggNogMapperPredictor()
        def outputDir = mockDir('eggnog')

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.size() == 3
        // protein_001: GO + COG + KEGG
        results['protein_001'].any { it.type == AnnotationType.GO && it.value == 'GO:0003735' }
        results['protein_001'].any { it.type == AnnotationType.GO && it.value == 'GO:0006412' }
        results['protein_001'].any { it.type == AnnotationType.COG && it.value == 'J' }
        results['protein_001'].any { it.type == AnnotationType.KEGG && it.value == 'ko:K02950' }
        // protein_002: EC number
        results['protein_002'].any { it.type == AnnotationType.EC && it.value == 'EC:3.1.26.3' }
        // protein_003: EC + KEGG
        results['protein_003'].any { it.type == AnnotationType.EC && it.value == 'EC:3.4.21.5' }
        results['protein_003'].any { it.type == AnnotationType.KEGG && it.value == 'ko:K01313' }
    }

    // --- gapseq ---

    def "gapseq: should parse Reactions.tbl"() {
        given:
        def predictor = new GapseqPredictor()
        def outputDir = mockDir('gapseq')

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.size() > 0
        // Active reactions should be parsed
        results['gene_001']?.any { it.type == AnnotationType.EC && it.value == 'EC:2.7.1.1' }
        results['gene_005']?.any { it.value == 'EC:4.2.1.11' }
        // Inactive reactions (gene = '-') should not be included
        !results.containsKey('-')
    }

    // --- SignalP ---

    def "SignalP: should parse prediction_results.txt"() {
        given:
        def predictor = new SignalPPredictor()
        def outputDir = mockDir('signalp')

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        // protein_001: SP signal peptide
        results['protein_001'].any { it.type == AnnotationType.SIGNAL_PEPTIDE && it.value == 'SP' }
        results['protein_001'].any { it.type == AnnotationType.SUBCELLULAR_LOCALIZATION }
        results['protein_001'].find { it.type == AnnotationType.SIGNAL_PEPTIDE }.score == 0.9823

        // protein_002: OTHER (no signal peptide) — should not appear
        !results.containsKey('protein_002')

        // protein_003: LIPO
        results['protein_003'].any { it.type == AnnotationType.SIGNAL_PEPTIDE && it.value == 'LIPO' }

        // protein_004: TAT
        results['protein_004'].any { it.type == AnnotationType.SIGNAL_PEPTIDE && it.value == 'TAT' }

        // protein_005: OTHER — should not appear
        !results.containsKey('protein_005')
    }

    // --- DeepTMHMM ---

    def "DeepTMHMM: should parse predicted_topologies.3line"() {
        given:
        def predictor = new DeepTmhmmPredictor()
        def outputDir = mockDir('deeptmhmm')

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        // protein_001: no TM helices (all 'o')
        !results.containsKey('protein_001')

        // protein_002: 2 TM helices
        results['protein_002'].any { it.type == AnnotationType.TRANSMEMBRANE }
        results['protein_002'].find { it.type == AnnotationType.TRANSMEMBRANE }.metadata.tmHelixCount == 2
        results['protein_002'].any { it.type == AnnotationType.SUBCELLULAR_LOCALIZATION && it.value == 'integral membrane' }

        // protein_003: no TM helices
        !results.containsKey('protein_003')

        // protein_004: 3 TM helices
        results['protein_004'].find { it.type == AnnotationType.TRANSMEMBRANE }.metadata.tmHelixCount == 3
    }

    // --- antiSMASH ---

    def "antiSMASH: should parse JSON output"() {
        given:
        def predictor = new AntiSmashPredictor()
        def outputDir = mockDir('antismash')

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        results.size() > 0
        // NRPS/T1PKS cluster
        results['gene_010']?.any { it.type == AnnotationType.SECONDARY_METABOLITE && it.value.contains('NRPS') }
        results['gene_011']?.any { it.type == AnnotationType.SECONDARY_METABOLITE }
        results['gene_012']?.any { it.type == AnnotationType.SECONDARY_METABOLITE }
        // Terpene cluster
        results['gene_025']?.any { it.type == AnnotationType.SECONDARY_METABOLITE && it.value == 'terpene' }
    }

    // --- CRISPR/minced ---

    def "minced: should parse GFF output for CRISPR arrays"() {
        given:
        def predictor = new CrisprPredictor()
        def outputDir = mockDir('minced')

        when:
        def results = predictor.parseOutput(outputDir)

        then:
        // Two CRISPR arrays on two different contigs
        results['contig_1']?.any { it.type == AnnotationType.CRISPR }
        results['contig_2']?.any { it.type == AnnotationType.CRISPR }
        results['contig_1'][0].metadata.start == 10245
        results['contig_1'][0].metadata.end == 11456
    }

    // --- buildCommand tests ---

    def "all predictors should generate valid command structures"() {
        given:
        def inputFasta = File.createTempFile('test_', '.faa')
        def outputDir = File.createTempDir('test_', '')
        inputFasta.text = '>test\nMKAIL\n'

        expect:
        predictor.buildCommand(inputFasta, outputDir).size() > 0
        predictor.buildCommand(inputFasta, outputDir)[0] != null

        cleanup:
        inputFasta.delete()
        outputDir.deleteDir()

        where:
        predictor << [
            new gspa.predictor.similarity.DiamondPredictor(database: '/tmp/db'),
            new MMseqs2Predictor(database: '/tmp/db'),
            new FoldSeekPredictor(database: '/tmp/db'),
            new gspa.predictor.domain.InterProScanPredictor(),
            new HmmerPredictor(hmmDatabase: '/tmp/Pfam-A.hmm'),
            new SignalPPredictor(),
            new DeepTmhmmPredictor(),
            new gspa.predictor.specialized.AmrFinderPredictor(),
            new gspa.predictor.specialized.DbCanPredictor(),
            new CrisprPredictor(),
            new gspa.predictor.specialized.VfdbPredictor(database: '/tmp/vfdb'),
        ]
    }

    def "DIAMOND command should include all required flags"() {
        given:
        def predictor = new gspa.predictor.similarity.DiamondPredictor(
            database: '/data/swissprot.dmnd',
            evalue: 1e-10,
            maxTargetSeqs: 5,
        )
        def input = File.createTempFile('test_', '.faa')
        def outDir = File.createTempDir('test_', '')

        when:
        def cmd = predictor.buildCommand(input, outDir)

        then:
        cmd.contains('blastp')
        cmd.contains('--db')
        cmd.contains('/data/swissprot.dmnd')
        cmd.contains('--evalue')
        cmd.contains('1.0E-10')
        cmd.contains('--max-target-seqs')
        cmd.contains('5')

        cleanup:
        input.delete()
        outDir.deleteDir()
    }

    def "InterProScan command should include goterms and applications"() {
        given:
        def predictor = new gspa.predictor.domain.InterProScanPredictor(
            applications: ['Pfam', 'TIGRFAM'],
            goterms: true,
        )
        def input = File.createTempFile('test_', '.faa')
        def outDir = File.createTempDir('test_', '')

        when:
        def cmd = predictor.buildCommand(input, outDir)

        then:
        cmd.contains('-appl')
        cmd.contains('Pfam,TIGRFAM')
        cmd.contains('-goterms')
        cmd.contains('-f')
        cmd.contains('TSV')

        cleanup:
        input.delete()
        outDir.deleteDir()
    }

    // --- Helper ---

    private File mockDir(String tool) {
        def resourcePath = "/mock-output/${tool}"
        def url = getClass().getResource(resourcePath)
        if (url) {
            return new File(url.toURI())
        }
        throw new RuntimeException("Mock output not found: ${resourcePath}")
    }
}
