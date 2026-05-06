/*
 * Neural function predictors.
 *
 * Each process wraps benchmark/neural/run_neural_predictors.py with a fixed
 * --predictor flag and emits the standard 4-column TSV
 *   protein_id<TAB>term<TAB>score<TAB>annotation_type
 *
 * Channel I/O matches the rest of gspa-nf: tuple(sample_id, file).
 *
 * Containers (set in nextflow.config via withName): four of these processes
 * share the ESM2 PyTorch image; PROTEINFER uses the TensorFlow image.
 *
 * GPU is opt-in via the `gpu` profile (containerOptions = '--gpus all' on
 * docker, '--nv' on singularity). All processes also work CPU-only at
 * reduced throughput; ESM2_DEEPGOPLUS and CLEAN want a GPU in practice.
 */

def write_manifest(sample_id, fasta) {
    """tag\tfasta_path\toutput_dir
${sample_id}\t${fasta}\t."""
}

process ESM2_DEEPGOPLUS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/esm2_deepgoplus", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path checkpoint
    path terms

    output:
    tuple val(sample_id), path("${sample_id}.esm2-deepgoplus.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_neural_predictors.py \\
        --predictor esm2-deepgoplus \\
        --manifest manifest.tsv \\
        --checkpoint ${checkpoint} \\
        --terms ${terms} \\
        --esm2-model ${params.esm2_dgp_model} \\
        --batch-size ${params.esm2_dgp_batch_size} \\
        --min-score ${params.neural_min_score}
    """
}

process ESM2_CENTROID {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/esm2_centroid", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path centroid_db

    output:
    tuple val(sample_id), path("${sample_id}.esm2-centroid.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_neural_predictors.py \\
        --predictor esm2-centroid \\
        --manifest manifest.tsv \\
        --centroid-db ${centroid_db} \\
        --esm2-model ${params.esm2_centroid_model} \\
        --top-k ${params.esm2_centroid_top_k} \\
        --batch-size ${params.esm2_centroid_batch_size} \\
        --min-score ${params.neural_min_score}
    """
}

process CLEAN {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/clean", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path model_dir

    output:
    tuple val(sample_id), path("${sample_id}.clean.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_neural_predictors.py \\
        --predictor clean \\
        --manifest manifest.tsv \\
        --model-dir ${model_dir} \\
        --batch-size ${params.clean_batch_size} \\
        --min-score ${params.neural_min_score}
    """
}

process PROTEINFER {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/proteinfer", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path model_dir

    output:
    tuple val(sample_id), path("${sample_id}.proteinfer.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_neural_predictors.py \\
        --predictor proteinfer \\
        --manifest manifest.tsv \\
        --model-dir ${model_dir} \\
        --batch-size ${params.proteinfer_batch_size} \\
        --min-score ${params.neural_min_score}
    """
}
