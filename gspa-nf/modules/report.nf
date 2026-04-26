/*
 * Per-sample multi-format report.
 *
 * Consumes: every neural predictor TSV (with its name) for a sample, the
 * ensemble TSV (if present), and the per-predictor eval JSON (if EVAL_PGAP
 * ran). Channels are joined on sample_id; each predictor enters as a
 * (name, file) tuple so the script's repeatable --predictor flag scales
 * to N predictors with no module change.
 *
 * Emits: <sample>.{html,ttl,jsonld} under
 *   ${outdir}/${sample_id}/report/
 *
 * Vocabulary: SIO-aligned; see benchmark/neural/make_report.py docstring.
 */

process MAKE_REPORT {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/report", mode: 'copy'

    input:
    // pred_specs: list of "name:filename" strings; pred_files: matching paths
    // eval_specs: same shape for eval json (may be empty list)
    tuple val(sample_id),
          val(pred_specs),  path(pred_files,  stageAs: 'pred_*'),
          val(eval_specs),  path(eval_files,  stageAs: 'eval_*')

    output:
    tuple val(sample_id),
          path("${sample_id}.html"),
          path("${sample_id}.ttl"),
          path("${sample_id}.jsonld"), emit: results

    script:
    def pred_args = pred_specs.collect { spec -> "--predictor ${spec}" }.join(' ')
    def eval_args = eval_specs.collect { spec -> "--eval ${spec}" }.join(' ')
    """
    python3 /opt/gspa/make_report.py \\
        --sample-id ${sample_id} \\
        ${pred_args} \\
        ${eval_args} \\
        --out-dir . \\
        --min-score ${params.report_min_score}
    """
}
