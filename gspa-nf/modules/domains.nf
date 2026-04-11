process HMMSEARCH {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/hmmer", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path pfam_db
    path pfam_h3f
    path pfam_h3i
    path pfam_h3m
    path pfam_h3p

    output:
    tuple val(sample_id), path("${sample_id}_pfam.domtbl"), emit: results

    script:
    """
    hmmsearch \\
        --domtblout ${sample_id}_pfam.domtbl \\
        --noali \\
        -E 1e-5 \\
        --cpu ${task.cpus} \\
        ${pfam_db} \\
        ${proteins} \\
        > /dev/null
    """
}

process EGGNOG_MAPPER {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/eggnog", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path eggnog_db

    output:
    tuple val(sample_id), path("${sample_id}_eggnog.emapper.annotations"), emit: results

    script:
    """
    emapper.py \\
        -i ${proteins} \\
        -o ${sample_id}_eggnog \\
        --data_dir ${eggnog_db} \\
        -m diamond \\
        --cpu ${task.cpus} \\
        --itype proteins \\
        --override
    """
}

process DBCAN {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/dbcan", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path dbcan_db

    output:
    tuple val(sample_id), path("${sample_id}_dbcan_overview.tsv"), emit: results

    script:
    """
    run_dbcan CAZyme_annotation \\
        --input_raw_data ${proteins} \\
        --mode protein \\
        --output_dir dbcan_out \\
        --db_dir ${dbcan_db}
    cp dbcan_out/overview.tsv ${sample_id}_dbcan_overview.tsv || touch ${sample_id}_dbcan_overview.tsv
    """
}
