package gspa.cli

import gspa.integration.ClaimExtractor
import gspa.predictor.neural.CleanPredictor
import com.fasterxml.jackson.databind.ObjectMapper
import picocli.CommandLine
import spock.lang.Specification
import spock.lang.TempDir
import java.security.MessageDigest

class CleanClaimsSpec extends Specification {
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
            'p|one\tEC:1.1.1.1\t'+(kind=='distance'?'0.0':Double.toString(nativeScore))+'\tEC\t'+new ObjectMapper().writeValueAsString(metadata)+'\n'
        Map receipt = metadata + [schema:'gspa.clean-output/1', input_protein_ids:['p|one'],
            input_state:'nonempty',
            native_to_input_ids:[opaque:'p|one'],
            input_proteins:1, native_proteins:1, native_selected_calls:1, emitted_calls:1,
            filtered_calls:0, min_score:0d, output_sha256:sha(tsv)]
        File receiptFile = new File(directory,'query.clean.receipt.json')
        receiptFile.text=new ObjectMapper().writeValueAsString(receipt)
        [tsv:tsv, receipt:receipt, receiptFile:receiptFile, metadata:metadata, input:input]
    }


    private int invoke(Map files, File output, String kind) {
        new CommandLine(new GspaMain()).setErr(new PrintWriter(new StringWriter())).execute(
            'clean-claims', '--tsv', files.tsv.path, '--receipt', files.receiptFile.path,
            '--proteins', files.input.path, '--score-kind', kind, '--output', output.path)
    }

    def 'typed CLEAN import matches annotation claims for #kind'() {
        given:
        def files=fixture(kind, value)
        File output=new File(directory,'claims.jsonl')
        File expected=new File(directory,'expected.jsonl')
        def extractor=new ClaimExtractor()
        def annotations=new CleanPredictor(scoreKind:kind).parseOutput(directory)
        ClaimExtractor.writeClaimsJsonl(annotations['p|one'].collect { extractor.fromAnnotation('p|one',it) },expected)
        when:
        int result=invoke(files,output,kind)
        then:
        result==0
        output.text==expected.text
        new File(output.path+'.import.json').exists()
        extractor.readClaimsJsonl(output)[0].hasFunctionalScore()==(kind=='gmm_confidence')
        where:
        kind             | value
        'gmm_confidence' | 0.00000001d
        'distance'       | 8.0627d
    }

    def 'invalid or replaced CLEAN inputs do not produce claims (#fault)'() {
        given:
        def files=fixture()
        File output=new File(directory,'claims.jsonl')
        switch(fault) {
            case 'changed-tsv': files.tsv << 'bad row\n'; break
            case 'changed-fasta': files.input.text='>p|one\nDIFFERENT\n'; break
            case 'missing-receipt': files.receiptFile.delete(); break
            case 'missing-fasta': files.input.delete(); break
        }
        expect:
        invoke(files,output,'gmm_confidence')!=0
        !output.exists()
        !new File(output.path+'.import.json').exists()
        where:
        fault << ['changed-tsv','changed-fasta','missing-receipt','missing-fasta']
    }

    def 'import never overwrites an existing claim artifact'() {
        given:
        def files=fixture()
        File output=new File(directory,'claims.jsonl');output.text='existing\n'
        expect:
        invoke(files,output,'gmm_confidence')!=0
        output.text=='existing\n'
        !new File(output.path+'.import.json').exists()
    }
}
