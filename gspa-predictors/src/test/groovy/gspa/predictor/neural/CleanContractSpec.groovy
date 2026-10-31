package gspa.predictor.neural

import gspa.config.GspaConfig
import gspa.integration.ClaimExtractor
import gspa.model.AnnotationType
import gspa.predictor.AnnotationPipeline
import gspa.predictor.ToolOutputContractException
import groovy.json.JsonOutput
import spock.lang.Specification
import spock.lang.TempDir
import java.security.MessageDigest

class CleanContractSpec extends Specification {
    @TempDir File directory

    private static String sha(File file) {
        MessageDigest.getInstance('SHA-256').digest(file.bytes).encodeHex().toString()
    }

    private Map fixture(String kind='gmm_confidence', double nativeScore=0.8d) {
        File input = new File(directory, 'input.faa')
        input.text = '>p|one\nMALW\n'
        Map metadata = [score_kind:kind, score_available:kind=='gmm_confidence', score_calibrated:false,
            native_score:nativeScore, native_score_text:Double.toString(nativeScore),
            native_protein_id:'opaque',
            native_score_direction:kind=='distance'?'lower':'higher',
            predictor_contract:'clean:maxsep/1:'+kind+':'+('d'*64), driver_sha256:'d'*64,
            driver_score_kind:kind, native_output_sha256:'a'*64, input_fasta_sha256:sha(input),
            candidate_coverage:'native_maxsep_selected_only; not the complete EC distance matrix']
        File tsv = new File(directory,'query.clean.tsv')
        tsv.text = 'protein_id\tterm\tscore\tannotation_type\tmetadata_json\n' +
            'p|one\tEC:1.1.1.1\t'+(kind=='distance'?'0.0':Double.toString(nativeScore))+'\tEC\t'+JsonOutput.toJson(metadata)+'\n'
        Map receipt = metadata + [schema:'gspa.clean-output/1', input_protein_ids:['p|one'],
            input_state:'nonempty',
            native_to_input_ids:[opaque:'p|one'],
            input_proteins:1, native_proteins:1, native_selected_calls:1, emitted_calls:1,
            filtered_calls:0, min_score:0d, output_sha256:sha(tsv)]
        File receiptFile = new File(directory,'query.clean.receipt.json')
        receiptFile.text=JsonOutput.toJson(receipt)
        [tsv:tsv, receipt:receipt, receiptFile:receiptFile, metadata:metadata]
    }

    def 'native confidence keeps exact score and typed provenance through raw claims'() {
        given:
        fixture('gmm_confidence', 0.00000001d)
        def predictor = new CleanPredictor(scoreKind:'gmm_confidence')
        when:
        def annotation = predictor.parseOutput(directory)['p|one'][0]
        def claim = new ClaimExtractor().fromAnnotation('p|one',annotation)
        File raw = new File(directory,'claims.jsonl')
        ClaimExtractor.writeClaimsJsonl([claim],raw)
        def restored = new ClaimExtractor().readClaimsJsonl(raw)[0]
        then:
        annotation.score == 0.00000001d
        annotation.type == AnnotationType.EC
        restored.rawScore == annotation.score
        restored.metadata.native_score == annotation.score
        restored.metadata.score_kind == 'gmm_confidence'
        restored.hasFunctionalScore()
    }

    def 'provisional EC identifiers from the native label vocabulary are preserved'() {
        given:
        def files=fixture()
        files.tsv.text=files.tsv.text.replace('EC:1.1.1.1','EC:1.1.1.n11')
        files.receipt.output_sha256=sha(files.tsv)
        files.receiptFile.text=JsonOutput.toJson(files.receipt)
        expect:
        new CleanPredictor(scoreKind:'gmm_confidence').parseOutput(directory)['p|one'][0].value=='EC:1.1.1.n11'
    }

    def 'native distance remains available without invented functional confidence'() {
        given:
        fixture('distance', 8.0627d)
        when:
        def annotation = new CleanPredictor(scoreKind:'distance').parseOutput(directory)['p|one'][0]
        def claim = new ClaimExtractor().fromAnnotation('p|one',annotation)
        then:
        annotation.metadata.native_score == 8.0627d
        annotation.metadata.native_score_direction == 'lower'
        claim.calibratedProb == null
        !claim.hasFunctionalScore()
    }

    def 'malformed or mismatched native output is rejected (#change)'() {
        given:
        def files = fixture()
        switch(change) {
            case 'header': files.tsv.text=files.tsv.text.replace('protein_id','other'); break
            case 'unknown': files.tsv.text=files.tsv.text.replace('p|one\t','not-input\t'); break
            case 'type': files.tsv.text=files.tsv.text.replace('\tEC\t','\tGO\t'); break
            case 'range': files.tsv.text=files.tsv.text.replace('\t0.8\t','\t8.0\t'); break
            case 'duplicate': files.tsv << files.tsv.readLines()[1]+'\n'; break
            case 'missing-row': files.tsv.text=files.tsv.readLines()[0]+'\n'; break
            case 'changed-input': new File(directory,'input.faa').text='>p|one\nOTHER\n'; break
            case 'wrong-mode': files.receipt.score_kind='distance'; break
            case 'missing-receipt': files.receiptFile.delete(); break
        }
        if (change != 'missing-receipt') {
            files.receipt.output_sha256=sha(files.tsv)
            files.receiptFile.text=JsonOutput.toJson(files.receipt)
        }
        when:
        new CleanPredictor(scoreKind:'gmm_confidence').parseOutput(directory)
        then:
        thrown(ToolOutputContractException)
        where:
        change << ['header','unknown','type','range','duplicate','missing-row','changed-input','wrong-mode','missing-receipt']
    }

    def 'zero cutoff cannot disguise missing calls as filtering'() {
        given:
        def files=fixture()
        files.tsv.text=files.tsv.readLines()[0]+'\n'
        files.receipt.output_sha256=sha(files.tsv)
        files.receipt.emitted_calls=0;files.receipt.filtered_calls=1
        files.receiptFile.text=JsonOutput.toJson(files.receipt)
        when:
        new CleanPredictor(scoreKind:'gmm_confidence').parseOutput(directory)
        then:
        thrown(ToolOutputContractException)
    }

    def 'clean registration forwards an explicit score kind and retains native selected weak calls by default'() {
        given:
        def config = new GspaConfig()
        config.predictors.neural.clean.enabled=true
        config.predictors.neural.clean.modelDir='/models/CLEAN/app'
        config.predictors.neural.clean.scoreKind='gmm_confidence'
        config.predictors.neural.sidecarScript='/tools/run_neural_predictors.py'
        def pipeline = new AnnotationPipeline(config)
        when:
        pipeline.configure()
        def predictor = pipeline.registry.get('clean')
        def command = predictor.buildCommand(new File(directory,'input.faa'),directory)
        then:
        predictor.minScore == 0d
        command[command.indexOf('--clean-score-kind')+1]=='gmm_confidence'
    }

    def 'clean cannot silently infer score semantics'() {
        when:
        new CleanPredictor(modelDir:'/models/CLEAN/app',sidecarScript:'/tools/sidecar.py')
            .buildCommand(new File(directory,'input.faa'),directory)
        then:
        thrown(IllegalArgumentException)
    }
}
