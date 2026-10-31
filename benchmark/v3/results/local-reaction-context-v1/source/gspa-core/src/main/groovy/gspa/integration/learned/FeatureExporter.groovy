package gspa.integration.learned

import groovy.json.JsonOutput
import gspa.integration.EvidenceClaim
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.model.AnnotationType

import java.security.DigestOutputStream
import java.security.MessageDigest

/**
 * Streams candidate feature rows from the SAME extractor the learned stage scores with, so
 * training rows and inference vectors cannot drift. One JSON object per line:
 *
 * <pre>
 * {"schema":"gspa.context-feature-row/1","contract":{"id":..,"fingerprint":..},
 *  "genome_id":..,"protein_id":..,"function_type":"EC","function_id":"EC:1.1.1.1",
 *  "function_key":..,"role":"CATALYTIC","candidate_key":..,
 *  "features":{<every contract column>: value},
 *  "available_contexts":{operon,layout,gap,neighborhood,cohort,assembly,assembly_completeness,
 *                        fba,fba_suppressed,predictor_status,sequence_evidence : bool},
 *  "notes":[..],
 *  "model":{"id","sha256","kind","margin","probability","display_threshold","above_display_threshold"}}
 * </pre>
 *
 * <p>The {@code model} object is present only when a frozen model is supplied; export needs no
 * trained artifact, which is how real training data is acquired in the first place. Feature
 * values are written with {@link Double#toString} (shortest text that round-trips the exact
 * double), never rounded.
 *
 * <p><b>No labels.</b> The exporter reads no truth set and writes no {@code label}, {@code split},
 * {@code group} or {@code example_id}: an absent function is never auto-labelled negative. A
 * dataset builder joins independently curated truth and declared family/genome splits to the rows
 * by {@code candidate_key}, then feeds {@code scripts/train_context_model.py}, whose row schema
 * (`features` object keyed exactly by the contract columns, `role`, `genome_id`) these rows follow.
 *
 * <p>Streaming: claims are indexed once and one function's candidates are held at a time;
 * memory does not grow with the number of rows written.
 */
class FeatureExporter {

    static final String ROW_SCHEMA = 'gspa.context-feature-row/1'
    static final String MANIFEST_SCHEMA = 'gspa.context-feature-export/1'
    static final String CONTRACT_FILE = 'feature_contract.json'
    static final String ROWS_FILE = 'candidate_features.jsonl'
    static final String MANIFEST_FILE = 'export_manifest.json'

    /** Result of one export. */
    static class Summary {
        int rows
        int functions
        int rowsWithSequenceEvidence
        String contractId
        String contractFingerprint
        String modelId        // null when exported without a model
        String modelSha256
        String rowsSha256     // set by exportToDirectory
    }

    private final ContextFeatureExtractor extractor
    private final LearnedContextModel model

    /** @param model optional frozen model; when given each row also carries the model's score. */
    FeatureExporter(LearnedContextModel model = null, ContextFeatureExtractor extractor = new ContextFeatureExtractor()) {
        this.model = model
        this.extractor = extractor
        if (model != null && model.contract.fingerprint() != extractor.contract.fingerprint()) {
            throw new IllegalArgumentException('model and extractor use different feature contracts')
        }
    }

    /**
     * Stream rows to {@code out}.
     *
     * @param retained optional annotation set whose retained keys are also candidates (the same
     *        universe {@link LearnedIntegrationStage} rescores); null = claims and context only
     * @param ctx optional structured context; its candidate keys are exported unless
     *        {@code includeContextCandidates} is false
     * @param genomeId provenance only; defaults to the state's genome id
     */
    Summary export(List<EvidenceClaim> claims, IntegrationState state, ExternalContext ctx, Writer out,
                   IntegratedAnnotationSet retained = null, boolean includeContextCandidates = true, String genomeId = null) {
        if (state == null) throw new IllegalArgumentException('state must not be null')
        List<String> extra = []
        if (retained != null) extra.addAll(retained.provenance.keySet())
        if (includeContextCandidates && ctx != null) extra.addAll(ctx.candidateKeys())
        ContextFeatureExtractor.Prepared prepared = extractor.prepare(claims, state, ctx, extra)
        String gid = genomeId ?: state.genome?.id
        String head = '{"schema":' + JsonOutput.toJson(ROW_SCHEMA) +
            ',"contract":{"id":' + JsonOutput.toJson(extractor.contract.id) + ',"fingerprint":' + JsonOutput.toJson(extractor.contract.fingerprint()) + '}'
        List<String> names = extractor.contract.names()
        Summary s = new Summary(contractId: extractor.contract.id, contractFingerprint: extractor.contract.fingerprint(),
            modelId: model?.modelId, modelSha256: model?.sha256)
        s.functions = prepared.functionKeys().size()
        for (String key : prepared.functionKeys()) {
            for (CandidateFeatures f : prepared.candidatesOf(key)) {
                f.requireFinite()
                StringBuilder sb = new StringBuilder(2048)
                sb.append(head)
                if (gid != null) sb.append(',"genome_id":').append(JsonOutput.toJson(gid))
                sb.append(',"protein_id":').append(JsonOutput.toJson(f.proteinId))
                sb.append(',"function_type":').append(JsonOutput.toJson(f.functionType.name()))
                sb.append(',"function_id":').append(JsonOutput.toJson(f.functionId))
                sb.append(',"function_key":').append(JsonOutput.toJson(f.functionKey()))
                sb.append(',"role":').append(JsonOutput.toJson(f.role.name()))
                sb.append(',"candidate_key":').append(JsonOutput.toJson(f.candidateKey()))
                sb.append(',"sequence_sources":').append(JsonOutput.toJson(f.sequenceSources))
                if (f.qualityEvidence) sb.append(',"quality_evidence":').append(JsonOutput.toJson(f.qualityEvidence))
                sb.append(',"features":{')
                for (int i = 0; i < names.size(); i++) {
                    if (i > 0) sb.append(',')
                    sb.append(JsonOutput.toJson(names[i])).append(':').append(Double.toString(f.values[i]))
                }
                sb.append('},"available_contexts":').append(contexts(f))
                sb.append(',"notes":').append(JsonOutput.toJson(f.provenance))
                if (model != null) {
                    CandidateScore sc = null
                    String unsupported
                    try {
                        unsupported = model.unsupportedReason(f)
                        if (!unsupported) sc = model.score(f)
                    } catch (IllegalArgumentException e) {
                        throw new LearnedStageException("cannot score ${f.candidateKey()} for export: ${e.message}", e)
                    }
                    sb.append(',"model":{"id":').append(JsonOutput.toJson(model.modelId))
                    sb.append(',"sha256":').append(JsonOutput.toJson(model.sha256))
                    sb.append(',"kind":').append(JsonOutput.toJson(model.kind))
                    sb.append(',"score_status":').append(JsonOutput.toJson(unsupported ? 'unsupported' : 'scored'))
                    sb.append(',"unsupported_reason":').append(JsonOutput.toJson(unsupported))
                    sb.append(',"margin":').append(JsonOutput.toJson(sc?.margin))
                    sb.append(',"probability":').append(JsonOutput.toJson(sc?.probability))
                    sb.append(',"display_threshold":').append(Double.toString(model.displayThreshold))
                    sb.append(',"above_display_threshold":').append(JsonOutput.toJson(sc?.aboveDisplayThreshold)).append('}')
                }
                sb.append("}\n")
                out.write(sb.toString())
                s.rows++
                if (f.get('seq.n_groups_reporting') > 0) s.rowsWithSequenceEvidence++
            }
        }
        out.flush()
        s
    }

    /** Which families of context were AVAILABLE (not merely zero) for this candidate. */
    private static String contexts(CandidateFeatures f) {
        Map<String, Boolean> m = new LinkedHashMap<>()
        m.sequence_evidence = f.get('seq.n_groups_reporting') > 0
        m.predictor_status = f.get('seq.cov.status_available') > 0
        m.operon = f.get('ctx.operon.available') > 0
        m.layout = f.get('ctx.layout.available') > 0
        m.gap = f.get('ctx.gap.available') > 0
        m.neighborhood = f.get('ctx.nbr.available') > 0
        m.cohort = f.get(CohortFeatureColumns.AVAILABLE) > 0
        m.local_reaction = f.get(CohortFeatureColumns.groupColumn('local')) > 0
        m.assembly = f.get('asm.available') > 0
        m.assembly_completeness = f.get('asm.completeness_available') > 0
        m.fba = f.get('fba.available') > 0
        m.fba_suppressed = f.get('fba.suppressed_objective_forced') > 0
        JsonOutput.toJson(m)
    }

    /**
     * Write {@code candidate_features.jsonl}, the matching {@code feature_contract.json} (the file
     * the Python trainer reads) and an {@code export_manifest.json} with counts and hashes into
     * {@code dir}. Existing files of these names are replaced.
     */
    Summary exportToDirectory(List<EvidenceClaim> claims, IntegrationState state, ExternalContext ctx, File dir,
                              IntegratedAnnotationSet retained = null, boolean includeContextCandidates = true, String genomeId = null) {
        dir.mkdirs()
        File rows = new File(dir, ROWS_FILE)
        MessageDigest md = MessageDigest.getInstance('SHA-256')
        Summary s
        rows.withOutputStream { os ->
            def dos = new DigestOutputStream(os, md)
            def w = new BufferedWriter(new OutputStreamWriter(dos, 'UTF-8'))
            s = export(claims, state, ctx, w, retained, includeContextCandidates, genomeId)
            w.flush()
        }
        s.rowsSha256 = md.digest().collect { String.format('%02x', it & 0xff) }.join()
        new File(dir, CONTRACT_FILE).setText(extractor.contract.toJson(), 'UTF-8')
        new File(dir, MANIFEST_FILE).setText(JsonOutput.prettyPrint(JsonOutput.toJson([
            schema: MANIFEST_SCHEMA, rows_file: ROWS_FILE, rows: s.rows, functions: s.functions,
            rows_with_sequence_evidence: s.rowsWithSequenceEvidence, rows_sha256: s.rowsSha256,
            contract_file: CONTRACT_FILE, contract_id: s.contractId, contract_fingerprint: s.contractFingerprint,
            model: model == null ? null : [id: model.modelId, sha256: model.sha256, kind: model.kind, intended_use: model.intendedUse],
            labels_included: false,
        ])) + '\n', 'UTF-8')
        s
    }
}
