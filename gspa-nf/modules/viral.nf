/*
 * v1.3 phage / prophage predictors. Each wraps
 * benchmark/neural/run_genomic_predictors.py with a fixed --predictor
 * flag. Output: 6-column genomic-region TSV
 *   contig_id  region_start  region_end  region_type  score  attributes
 *
 * Inputs come from PYRODIGAL.out.{genome, gff} so the pipeline does not
 * need a separate genome FASTA channel.
 *
 * Containers (set in nextflow.config via withName):
 *   GENOMAD  → quay.io/biocontainers/genomad
 *   CHECKV   → quay.io/biocontainers/checkv
 *   PHISPY   → quay.io/biocontainers/phispy
 * No GSPA-built image — all upstream biocontainers.
 */

def write_genomic_manifest(sample_id, genome, gff) {
    """tag\tgenome_fasta\tgff_path\toutput_dir
${sample_id}\t${genome}\t${gff ?: '-'}\t."""
}

process GENOMAD {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/viral/genomad", mode: 'copy'

    input:
    tuple val(sample_id), path(genome)
    path db_path

    output:
    tuple val(sample_id), path("${sample_id}.genomad.genomic.tsv"), emit: results

    script:
    def manifest_text = write_genomic_manifest(sample_id, genome, null)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_genomic_predictors.py \\
        --predictor genomad \\
        --manifest manifest.tsv \\
        --db-path ${db_path} \\
        --min-score ${params.genomad_min_score}
    """
}

process CHECKV {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/viral/checkv", mode: 'copy'

    input:
    tuple val(sample_id), path(genome)
    path db_path

    output:
    tuple val(sample_id), path("${sample_id}.checkv.genomic.tsv"), emit: results

    script:
    def manifest_text = write_genomic_manifest(sample_id, genome, null)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_genomic_predictors.py \\
        --predictor checkv \\
        --manifest manifest.tsv \\
        --db-path ${db_path} \\
        --threads ${task.cpus} \\
        --min-score ${params.checkv_min_score}
    """
}

process PHISPY {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/viral/phispy", mode: 'copy'

    input:
    tuple val(sample_id), path(genome_gbk)

    output:
    tuple val(sample_id), path("${sample_id}.phispy.genomic.tsv"), emit: results

    script:
    def manifest_text = write_genomic_manifest(sample_id, genome_gbk, null)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_genomic_predictors.py \\
        --predictor phispy \\
        --manifest manifest.tsv \\
        --min-score ${params.phispy_min_score}
    """
}

process VIRSORTER2 {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/viral/virsorter2", mode: 'copy'

    input:
    tuple val(sample_id), path(genome)
    path db_path

    output:
    tuple val(sample_id), path("${sample_id}.virsorter2.genomic.tsv"), emit: results

    script:
    def manifest_text = write_genomic_manifest(sample_id, genome, null)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_genomic_predictors.py \\
        --predictor virsorter2 \\
        --manifest manifest.tsv \\
        --db-path ${db_path} \\
        --threads ${task.cpus} \\
        --min-score ${params.virsorter2_min_score}
    """
}

process VIBRANT {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/viral/vibrant", mode: 'copy'

    input:
    tuple val(sample_id), path(genome)
    path db_path

    output:
    tuple val(sample_id), path("${sample_id}.vibrant.genomic.tsv"), emit: results

    script:
    def manifest_text = write_genomic_manifest(sample_id, genome, null)
    """
    cat > manifest.tsv <<'EOF'
${manifest_text}
EOF
    python3 /opt/gspa/run_genomic_predictors.py \\
        --predictor vibrant \\
        --manifest manifest.tsv \\
        --db-path ${db_path} \\
        --threads ${task.cpus} \\
        --min-score ${params.vibrant_min_score}
    """
}

/**
 * Flatten genomic-region calls (prophage/plasmid/viral_contig) to
 * per-CDS region annotations using the GFF gene coordinates. Emits a
 * 5-col region TSV that auto-feeds the existing per-protein region
 * pipeline (HTML "Regions" section, gspa:DisorderRegion-style RDF).
 *
 * For every CDS overlapping a predicted prophage region, one row with
 * region_type=prophage_cds spanning the full CDS protein length.
 */
process FLATTEN_PROPHAGE_CDS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/viral/flatten", mode: 'copy'

    input:
    tuple val(sample_id), path(genomic_tsv), path(gff)

    output:
    tuple val(sample_id), path("${sample_id}.prophage_cds.tsv"), emit: results

    script:
    """
    python3 - <<'PY'
import csv, re

# Parse the genomic-region TSV
regions = []
with open('${genomic_tsv}') as fh:
    next(fh, None)  # skip header
    for line in fh:
        f = line.rstrip('\\n').split('\\t')
        if len(f) < 5:
            continue
        contig, s, e, rtype, score = f[0], int(f[1]), int(f[2]), f[3], float(f[4])
        if rtype not in ('prophage', 'plasmid', 'viral_contig'):
            continue
        regions.append((contig, s, e, rtype, score))

# Walk the GFF and emit one row per CDS overlapping any region
with open('${sample_id}.prophage_cds.tsv', 'w') as out:
    out.write('protein_id\\tregion_start\\tregion_end\\tregion_type\\tscore\\n')
    if not regions:
        exit(0)
    with open('${gff}') as gh:
        for line in gh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\\n').split('\\t')
            if len(cols) < 9 or cols[2] != 'CDS':
                continue
            contig = cols[0]
            try:
                s, e = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            # Pull a protein_id from the attributes column
            m = re.search(r'(?:protein_id|ID)=([^;]+)', cols[8])
            if not m:
                continue
            pid = m.group(1)
            cds_len = e - s + 1
            for (rc, rs, re_, rt, sc) in regions:
                if rc == contig and not (e < rs or s > re_):
                    # Per-protein flat region spans full CDS length (1..cds_len)
                    out.write(f'{pid}\\t1\\t{cds_len // 3}\\t{rt}_cds\\t{sc:.4f}\\n')
                    break
PY
    """
}
