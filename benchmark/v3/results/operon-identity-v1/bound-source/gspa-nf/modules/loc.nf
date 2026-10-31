/*
 * PSORTb 3.0 — bacterial subcellular localization (GPL-3).
 * Uses upstream brinkmanlab/psortb_commandline image (no rebuild).
 * Output: 4-column TSV (joins ENSEMBLE_PREDS via ch_neural).
 */

def write_manifest(sample_id, fasta) {
    """tag\tfasta_path\toutput_dir
${sample_id}\t${fasta}\t."""
}

process PSORTB {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/loc/psortb", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}.psortb.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_term_predictors.py \\
        --predictor psortb \\
        --manifest manifest.tsv \\
        --min-score ${params.psortb_min_score} \\
        --gram ${params.psortb_gram}
    """
}
