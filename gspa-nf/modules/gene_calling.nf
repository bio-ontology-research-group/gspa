process PYRODIGAL {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/gene_calling", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}.gff"),  emit: gff
    tuple val(sample_id), path("${sample_id}.faa"),  emit: proteins
    tuple val(sample_id), path(fasta),                emit: genome

    script:
    def mode = params.mag ? '-p meta' : '-p single'
    """
    pyrodigal \\
        -i ${fasta} \\
        -o ${sample_id}.gff \\
        -a ${sample_id}.faa \\
        -f gff \\
        ${mode}
    """
}
