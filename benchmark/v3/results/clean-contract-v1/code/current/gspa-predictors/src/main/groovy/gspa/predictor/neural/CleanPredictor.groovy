package gspa.predictor.neural

import gspa.model.AnnotationType
import gspa.model.Annotation
import gspa.io.FastaReader
import gspa.predictor.ToolOutputContractException
import groovy.json.JsonSlurper
import java.security.MessageDigest

/**
 * Wraps CLEAN's ESM-1b/contrastive EC predictor. The native max-separation
 * file can contain lower-is-better distances or GMM confidence values.
 * Distance calls retain their native values as unscored evidence.
 *
 * <p>Emits only EC annotations.
 */
class CleanPredictor extends AbstractNeuralSidecarPredictor {

    /** Directory with a CLEAN checkpoint (contains CLEAN_infer_fasta.py). */
    String modelDir
    /** Required native contract; checked against the driver's explicit GMM argument. */
    String scoreKind

    CleanPredictor() { minScore = 0d }

    @Override
    String getName() { 'clean' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.EC] as Set }

    @Override
    String getPredictorName() { 'clean' }

    @Override
    protected List<String> extraSidecarArgs() {
        if (!modelDir) {
            throw new IllegalStateException(
                "${name}: modelDir is unset. Point at a CLEAN checkpoint directory.")
        }
        validateScoreKind()
        ['--model-dir', new File(modelDir).absolutePath, '--clean-score-kind', scoreKind]
    }

    private void validateScoreKind() {
        if (!(scoreKind in ['gmm_confidence', 'distance'])) {
            throw new IllegalArgumentException('CLEAN requires scoreKind=gmm_confidence or distance')
        }
        if (!Double.isFinite(minScore) || minScore < 0d || minScore > 1d ||
                (scoreKind == 'distance' && minScore != 0d)) {
            throw new IllegalArgumentException('CLEAN minScore is a confidence cutoff; distance mode requires zero')
        }
    }

    private static String sha(File file) {
        MessageDigest.getInstance('SHA-256').digest(file.bytes).encodeHex().toString()
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        File tsv = requireOutput(outputDir, "${tag}.clean.tsv", 'typed CLEAN output and receipt required')
        File receiptFile = requireOutput(outputDir, "${tag}.clean.receipt.json", 'complete native-input coverage required')
        parseFiles(tsv, receiptFile, new File(outputDir, 'input.faa'))
    }

    /** Shared import/execution parser; validates the bound FASTA, receipt and every row. */
    Map<String, List<Annotation>> parseFiles(File tsv, File receiptFile, File input) {
        validateScoreKind()
        if (!tsv.isFile() || !receiptFile.isFile() || !input.isFile()) {
            throw new ToolOutputContractException('CLEAN requires TSV, receipt and original protein FASTA')
        }
        def receipt = new JsonSlurper().parse(receiptFile)
        if (!(receipt instanceof Map) || receipt.schema != 'gspa.clean-output/1' ||
                receipt.score_kind != scoreKind || receipt.driver_score_kind != scoreKind ||
                receipt.output_sha256 != sha(tsv) || receipt.min_score != minScore ||
                !(receipt.input_protein_ids instanceof List) ||
                !(receipt.native_to_input_ids instanceof Map) ||
                receipt.native_to_input_ids.values().toSet() != receipt.input_protein_ids.toSet() ||
                receipt.native_to_input_ids.size() != receipt.input_protein_ids.size() ||
                receipt.input_protein_ids.any { !(it instanceof String) || !it } ||
                receipt.input_protein_ids.toSet().size() != receipt.input_protein_ids.size() ||
                receipt.input_proteins != receipt.input_protein_ids.size() ||
                receipt.native_proteins != receipt.input_proteins ||
                receipt.input_state != (receipt.input_proteins == 0 ? 'empty' : 'nonempty') ||
                receipt.native_selected_calls < receipt.input_proteins ||
                receipt.emitted_calls < 0 || receipt.filtered_calls < 0 ||
                (minScore == 0d && receipt.filtered_calls != 0) ||
                receipt.native_selected_calls != receipt.emitted_calls + receipt.filtered_calls ||
                receipt.score_available != (scoreKind == 'gmm_confidence') || receipt.score_calibrated != false) {
            throw new ToolOutputContractException('CLEAN output receipt is inconsistent')
        }
        Set<String> allowed = receipt.input_protein_ids as Set<String>
        def proteins = FastaReader.readProteins(input)
        if (receipt.input_fasta_sha256 != sha(input) || proteins.size() != allowed.size() ||
                (proteins*.id as Set) != allowed) {
            throw new ToolOutputContractException('CLEAN receipt does not match input proteins')
        }
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        Set<List<String>> seen = [] as Set
        int lineNumber = 0
        tsv.eachLine { line ->
            lineNumber++
            List<String> fields = line.split('\t', -1) as List<String>
            if (lineNumber == 1) {
                if (fields != ['protein_id', 'term', 'score', 'annotation_type', 'metadata_json']) {
                    throw new ToolOutputContractException('CLEAN typed TSV header is missing or changed')
                }
                return
            }
            if (fields.size() != 5 || !allowed.contains(fields[0]) || fields[3] != 'EC' ||
                    !(fields[1] ==~ /EC:[1-7](\.(\d+|-)){2}\.(\d+|n\d+|-)/) || !seen.add([fields[0], fields[1]])) {
                throw new ToolOutputContractException("CLEAN line ${lineNumber}: invalid identity, type or duplicate call")
            }
            double score = fields[2] as double
            def metadata = new JsonSlurper().parseText(fields[4])
            if (!(metadata instanceof Map) || !(metadata.native_score instanceof Number)) {
                throw new ToolOutputContractException("CLEAN line ${lineNumber}: missing native score metadata")
            }
            double nativeScore = (metadata.native_score as Number).doubleValue()
            List<String> shared = ['score_kind', 'score_available', 'score_calibrated', 'predictor_contract',
                'driver_sha256', 'driver_score_kind', 'native_output_sha256', 'input_fasta_sha256', 'candidate_coverage']
            if (!Double.isFinite(score) || !Double.isFinite(nativeScore) || nativeScore < 0d ||
                    receipt.native_to_input_ids[metadata.native_protein_id] != fields[0] ||
                    shared.any { metadata[it] != receipt[it] } ||
                    metadata.native_score_direction != (scoreKind == 'distance' ? 'lower' : 'higher') ||
                    !(metadata.native_score_text instanceof String) ||
                    Double.parseDouble(metadata.native_score_text as String) != nativeScore ||
                    (scoreKind == 'gmm_confidence' ? (nativeScore > 1d || score != nativeScore || score < minScore) : score != 0d)) {
                throw new ToolOutputContractException("CLEAN line ${lineNumber}: incompatible score semantics")
            }
            results[fields[0]] << new Annotation(type: AnnotationType.EC, value: fields[1],
                source: name, score: score, evidence: 'IEA', evidenceType: evidenceType, metadata: metadata)
        }
        if (lineNumber == 0 || seen.size() != receipt.emitted_calls ||
                (minScore == 0d && results.keySet() != allowed)) {
            throw new ToolOutputContractException('CLEAN output row count differs from receipt')
        }
        results
    }
}
