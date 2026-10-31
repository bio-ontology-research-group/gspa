/*
 * Evaluation: F-max + Smin + EC F-max for any 4-column predictor TSV
 * against a truth TSV (protein_id, aspect, function_id).
 *
 * Wraps benchmark/neural/evaluate_panel.py, which calls into
 * benchmark_pgap_v2 for the metric primitives.
 *
 * Opt-in via params.run_eval && params.truth_dir. The truth file for a
 * sample is matched on simpleName: ${sample_id}_truth.tsv in truth_dir.
 *
 * Output: one JSON record per (predictor, sample) — emitted as a TSV-like
 * .jsonl alongside the predictor TSV under results/<sample>/eval/.
 */

process EVAL_PGAP {
    tag "${sample_id}/${predictor}"
    publishDir "${params.outdir}/${sample_id}/eval", mode: 'copy'

    input:
    tuple val(sample_id), val(predictor), path(pred_tsv), path(truth_tsv)
    path go_aspect_map

    output:
    tuple val(sample_id), val(predictor), path("${sample_id}.${predictor}.eval.json"), emit: results

    script:
    def aspect_arg = go_aspect_map.name != 'NO_ASPECT_MAP' ? "--go-aspect-map ${go_aspect_map}" : ''
    def ann_type   = predictor == 'clean' ? 'EC' : 'GO'
    """
    python3 /opt/gspa/evaluate_panel.py \\
        --predictor-tsv ${pred_tsv} \\
        --truth ${truth_tsv} \\
        --annotation-type ${ann_type} \\
        ${aspect_arg} \\
        --tag ${sample_id} \\
        --predictor ${predictor} \\
        --truth-name ${params.eval_truth_name} \\
        --n-bootstrap ${params.eval_n_bootstrap} \\
        > ${sample_id}.${predictor}.eval.json
    """
}
