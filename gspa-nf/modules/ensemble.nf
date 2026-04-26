/*
 * Ensemble fusion of neural predictor outputs.
 *
 * Input: tuple(sample_id, [pred1.tsv, pred2.tsv, ...]) — collected via
 *        groupTuple from main.nf after mixing all enabled neural channels.
 * Output: tuple(sample_id, ensemble.tsv) with the same 4-column schema.
 *
 * Modes (params.ensemble_mode): max | mean | rank.
 * On the validated 21-genome panel, "mean" wins (RESULTS.md).
 */

process ENSEMBLE_PREDS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/ensemble", mode: 'copy'

    input:
    tuple val(sample_id), path(pred_tsvs, stageAs: 'pred_*.tsv')

    output:
    tuple val(sample_id), path("${sample_id}.ensemble.tsv"), emit: results

    script:
    def pred_args = pred_tsvs.collect { "--pred ${it}" }.join(' ')
    """
    python3 /opt/gspa/build_ensemble_preds.py \\
        ${pred_args} \\
        --out ${sample_id}.ensemble.tsv \\
        --mode ${params.ensemble_mode} \\
        --min-score ${params.ensemble_min_score}
    """
}
