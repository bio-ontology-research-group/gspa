/*
 * Region-level FOSS predictors. All wrap benchmark/neural/run_region_predictors.py
 * with a fixed --predictor flag. Output: 5-column TSV
 *   protein_id  region_start  region_end  region_type  score
 *
 * Container (nextflow.config): leechuck/gspa-region-stack:0.1
 */

def write_manifest(sample_id, fasta) {
    """tag\tfasta_path\toutput_dir
${sample_id}\t${fasta}\t."""
}

process METAPREDICT {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/region/metapredict", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}.metapredict.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_region_predictors.py \\
        --predictor metapredict \\
        --manifest manifest.tsv \\
        --min-score ${params.metapredict_min_score} \\
        --min-region-len ${params.metapredict_min_region_len}
    """
}

process DEEPSIG {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/region/deepsig", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}.deepsig.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_region_predictors.py \\
        --predictor deepsig \\
        --manifest manifest.tsv \\
        --min-score ${params.deepsig_min_score} \\
        --kingdom ${params.deepsig_kingdom}
    """
}

process TMBED {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/region/tmbed", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path prott5_model

    output:
    tuple val(sample_id), path("${sample_id}.tmbed.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    def model_arg = prott5_model.name != 'NO_PROTT5' ? "PROTT5_DIR=${prott5_model}" : ''
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    ${model_arg} python3 /opt/gspa/run_region_predictors.py \\
        --predictor tmbed \\
        --manifest manifest.tsv \\
        --min-score ${params.tmbed_min_score}
    """
}

process TPPRED3 {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/region/tppred3", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}.tppred3.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_region_predictors.py \\
        --predictor tppred3 \\
        --manifest manifest.tsv \\
        --min-score ${params.tppred3_min_score} \\
        --kingdom ${params.tppred3_kingdom}
    """
}
