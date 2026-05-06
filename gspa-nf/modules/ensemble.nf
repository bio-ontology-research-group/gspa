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

/*
 * Genomic-region ensemble — interval-overlap fusion across the per-tool
 * 6-col genomic-region TSVs (geNomad / CheckV / PhiSpy / VirSorter2 /
 * VIBRANT). Each input must arrive as 'NAME:FILENAME' so the fused
 * output records which tools voted for each consensus region.
 *
 * Default fusion is reciprocal-overlap clustering at 50% (params
 * .genomic_ensemble_min_overlap), with score = max | mean | wcov
 * (params.genomic_ensemble_mode, default 'max').
 */
process ENSEMBLE_GENOMIC {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/ensemble", mode: 'copy'

    input:
    tuple val(sample_id), val(specs), path(genomic_tsvs, stageAs: 'genomic_*.tsv')

    output:
    tuple val(sample_id), path("${sample_id}.ensemble.genomic.tsv"), emit: results

    script:
    def pred_args = specs.collect { "--pred ${it}" }.join(' ')
    """
    python3 /opt/gspa/build_genomic_ensemble.py \\
        --tag ${sample_id} \\
        ${pred_args} \\
        --out ${sample_id}.ensemble.genomic.tsv \\
        --mode ${params.genomic_ensemble_mode} \\
        --min-overlap ${params.genomic_ensemble_min_overlap} \\
        --min-score ${params.genomic_ensemble_min_score}
    """
}
