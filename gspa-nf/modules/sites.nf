/*
 * Site-level FOSS predictors. Output: 5-column TSV
 *   protein_id  position  site_type  score  annotation_type
 *
 * - MUSITEDEEP (MIT) — PTM sites (phospho-S/T/Y by default)
 * - SCANNET    (Apache-2): PPI interface residues; needs structures
 *
 * Site outputs do NOT join the ensemble (different shape) but DO appear
 * in the per-sample MAKE_REPORT output.
 */

def write_manifest(sample_id, fasta) {
    """tag\tfasta_path\toutput_dir
${sample_id}\t${fasta}\t."""
}

process MUSITEDEEP {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/sites/musitedeep", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path model_dir

    output:
    tuple val(sample_id), path("${sample_id}.musitedeep.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_site_predictors.py \\
        --predictor musitedeep \\
        --manifest manifest.tsv \\
        --model-dir ${model_dir} \\
        --residue-types ${params.musitedeep_residue_types} \\
        --min-score ${params.musitedeep_min_score}
    """
}

process SCANNET {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/sites/scannet", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path model_dir
    path structures      // dir with <sample_id>/*.pdb

    output:
    tuple val(sample_id), path("${sample_id}.scannet.tsv"), emit: results

    script:
    def manifest_text = write_manifest(sample_id, proteins)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_site_predictors.py \\
        --predictor scannet \\
        --manifest manifest.tsv \\
        --model-dir ${model_dir} \\
        --structure-dir ${structures} \\
        --min-score ${params.scannet_min_score}
    """
}
