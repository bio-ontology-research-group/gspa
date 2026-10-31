process DIAMOND_BLASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/diamond", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path db

    output:
    tuple val(sample_id), path("${sample_id}_diamond.tsv"), emit: results

    script:
    """
    diamond blastp \\
        --db ${db} \\
        --query ${proteins} \\
        --out ${sample_id}_diamond.tsv \\
        --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore stitle \\
        --evalue 1e-5 \\
        --max-target-seqs 10 \\
        --threads ${task.cpus} \\
        --query-cover 50 \\
        --subject-cover 50 \\
        --id 30
    """
}

process MMSEQS2_SEARCH {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/mmseqs2", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path db

    output:
    tuple val(sample_id), path("${sample_id}_mmseqs2.tsv"), emit: results

    script:
    """
    mkdir -p tmp
    mmseqs createdb ${proteins} queryDB
    mmseqs search queryDB ${db} resultDB tmp \\
        -e 1e-5 --max-seqs 10 --threads ${task.cpus}
    mmseqs convertalis queryDB ${db} resultDB ${sample_id}_mmseqs2.tsv \\
        --format-output "query,target,pident,alnlen,qlen,tlen,evalue,bits,theader"
    """
}

process FOLDSEEK {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/foldseek", mode: 'copy'

    input:
    tuple val(sample_id), path(query_input)
    path db
    path prostt5_model

    output:
    tuple val(sample_id), path("${sample_id}_foldseek.tsv"), emit: results

    script:
    // ProstT5 mode: use protein sequences directly (no structures needed)
    // Without ProstT5: needs PDB/mmCIF structure files as input
    def prostt5_flag = prostt5_model.name != 'NO_PROSTT5' ? "--prostt5-model ${prostt5_model}" : ''
    def format = prostt5_flag ? 'query,target,pident,evalue,bits,theader' :
                                'query,target,pident,evalue,bits,alntmscore,theader'
    """
    foldseek easy-search \\
        ${query_input} \\
        ${db} \\
        ${sample_id}_foldseek.tsv \\
        tmp \\
        ${prostt5_flag} \\
        --format-output "${format}" \\
        --threads ${task.cpus} \\
        -e 1e-3
    """
}
