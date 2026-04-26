/*
 * Additional term-level predictors that auto-join ENSEMBLE_PREDS:
 *   - DEEPFRI  (BSD-3): seq-only GO; complements ESM2-DGP / ProteInfer
 *   - DEEPEC   (AGPL): EC; complements CLEAN / DIAMOND (note network-clause)
 *   - DEEPARG  (MIT) : antimicrobial-resistance gene calls
 *
 * All wrap benchmark/neural/run_term_predictors.py and emit the standard
 * 4-column TSV (protein_id, term, score, annotation_type).
 */

def write_manifest(sample_id, fasta) {
    """tag\tfasta_path\toutput_dir
${sample_id}\t${fasta}\t."""
}

process DEEPFRI {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/term_extras/deepfri", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path model_dir

    output:
    tuple val(sample_id), path("${sample_id}.deepfri.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_term_predictors.py \\
        --predictor deepfri \\
        --manifest manifest.tsv \\
        --model-dir ${model_dir} \\
        --min-score ${params.deepfri_min_score}
    """
}

process DEEPEC {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/term_extras/deepec", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path model_dir

    output:
    tuple val(sample_id), path("${sample_id}.deepec.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_term_predictors.py \\
        --predictor deepec \\
        --manifest manifest.tsv \\
        --model-dir ${model_dir} \\
        --min-score ${params.deepec_min_score}
    """
}

process DEEPARG {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/term_extras/deeparg", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path model_dir

    output:
    tuple val(sample_id), path("${sample_id}.deeparg.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_term_predictors.py \\
        --predictor deeparg \\
        --manifest manifest.tsv \\
        --model-dir ${model_dir} \\
        --deeparg-type ${params.deeparg_type} \\
        --min-score ${params.deeparg_min_score}
    """
}
