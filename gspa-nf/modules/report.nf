/*
 * Per-sample multi-format report.
 *
 * Inputs (joined on sample_id by main.nf):
 *  - predictor + ensemble TSVs (4-col)  → --predictor name:filename
 *  - region TSVs                  (5-col) → --region    name:filename
 *  - site TSVs                    (5-col) → --site      name:filename
 *  - eval JSON files              (1 obj) → --eval      name:filename
 *
 * Each kind enters as a (specs, files) pair so the script's repeatable
 * flags scale to N predictors without a module change.
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
    tuple val(sample_id),
          val(pred_specs),     path(pred_files,    stageAs: 'pred_*'),
          val(region_specs),   path(region_files,  stageAs: 'region_*'),
          val(site_specs),     path(site_files,    stageAs: 'site_*'),
          val(genomic_specs),  path(genomic_files, stageAs: 'genomic_*'),
          val(eval_specs),     path(eval_files,    stageAs: 'eval_*')

    output:
    tuple val(sample_id),
          path("${sample_id}.html"),
          path("${sample_id}.ttl"),
          path("${sample_id}.jsonld"), emit: results

    script:
    def pred_args    = pred_specs.collect    { "--predictor ${it}" }.join(' ')
    def region_args  = region_specs.collect  { "--region ${it}" }.join(' ')
    def site_args    = site_specs.collect    { "--site ${it}" }.join(' ')
    def genomic_args = genomic_specs.collect { "--genomic-region ${it}" }.join(' ')
    def eval_args    = eval_specs.collect    { "--eval ${it}" }.join(' ')
    """
    python3 /opt/gspa/make_report.py \\
        --sample-id ${sample_id} \\
        ${pred_args} \\
        ${region_args} \\
        ${site_args} \\
        ${genomic_args} \\
        ${eval_args} \\
        --out-dir . \\
        --min-score ${params.report_min_score}
    """
}
