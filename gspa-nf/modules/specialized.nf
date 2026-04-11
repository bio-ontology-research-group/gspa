process BARRNAP {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/barrnap", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_rrna.gff"), emit: results

    script:
    def kingdom_flag = params.kingdom == 'archaea' ? '--kingdom arc' :
                       params.kingdom == 'eukaryote' ? '--kingdom euk' : '--kingdom bac'
    """
    barrnap ${kingdom_flag} --quiet ${fasta} > ${sample_id}_rrna.gff || true
    """
}

process MINCED {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/minced", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_crispr.gff"), emit: results

    script:
    """
    minced -minNR 3 ${fasta} ${sample_id}_crispr.gff || true
    """
}

process AMRFINDER {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/amrfinder", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path amrfinder_db

    output:
    tuple val(sample_id), path("${sample_id}_amr.tsv"), emit: results

    script:
    """
    amrfinder \\
        -p ${proteins} \\
        -o ${sample_id}_amr.tsv \\
        -d ${amrfinder_db} \\
        --threads ${task.cpus} \\
        --plus
    """
}

process ANTISMASH {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/antismash", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)
    path antismash_db

    output:
    tuple val(sample_id), path("antismash_out/${sample_id}.json"), emit: results

    script:
    """
    antismash \\
        ${fasta} \\
        --output-dir antismash_out \\
        --databases ${antismash_db} \\
        --cpus ${task.cpus} \\
        --minimal \\
        --genefinding-tool prodigal
    # Rename output to include sample_id
    mv antismash_out/*.json antismash_out/${sample_id}.json 2>/dev/null || true
    """
}

process SIGNALP {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/signalp", mode: 'copy'

    // No container — uses locally installed SignalP (requires license)
    container null

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${sample_id}_signalp.txt"), emit: results

    script:
    def org = params.kingdom == 'archaea' ? 'arch' :
              params.kingdom == 'eukaryote' ? 'eukarya' : 'gram-'
    def signalp = params.signalp_path ?: 'signalp'
    """
    ${signalp} \\
        --fasta ${proteins} \\
        --org ${org} \\
        --output_dir signalp_out \\
        --format short
    cp signalp_out/prediction_results.txt ${sample_id}_signalp.txt 2>/dev/null || touch ${sample_id}_signalp.txt
    """
}

process CHECKM2 {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/checkm2", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)
    path checkm2_db

    output:
    tuple val(sample_id), path("checkm2_out/quality_report.tsv"), emit: report

    script:
    """
    mkdir -p input_dir
    cp ${fasta} input_dir/
    checkm2 predict \\
        --input input_dir \\
        --output-directory checkm2_out \\
        --database_path ${checkm2_db} \\
        --threads ${task.cpus} \\
        --force \\
        -x fna
    """
}

process GTDBTK {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/gtdbtk", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)
    path gtdbtk_db

    output:
    tuple val(sample_id), path("gtdbtk_out/classify/*.summary.tsv"), emit: summary, optional: true

    script:
    """
    export GTDBTK_DATA_PATH=${gtdbtk_db}
    mkdir -p input_dir
    cp ${fasta} input_dir/
    gtdbtk classify_wf \\
        --genome_dir input_dir \\
        --out_dir gtdbtk_out \\
        --cpus ${task.cpus} \\
        -x fna \\
        --skip_ani_screen
    """
}
