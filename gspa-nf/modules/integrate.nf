/*
 * Phase 7 evidence integrator parity for the Nextflow path.
 *
 * Closes the gap noted in gspa-nf/README.md: the JVM CLI's `gspa
 * integrate` produces per-(protein, function) posterior probabilities
 * with full prior stack and provenance, while the Nextflow pipeline
 * historically stopped at the merged TSV. Two processes:
 *
 *   BUILD_CLAIMS  — wraps benchmark/02b_parse_predictors_to_claims.py
 *                   to lift the per-tool TSVs into a single claims.jsonl
 *   INTEGRATE     — invokes `gspa-cli integrate` on that claims.jsonl
 *                   plus reference data; emits ${sample_id}_integrated.tsv
 *
 * Both opt-in via params.run_integrate. Reference-data paths
 * (params.goa, params.go_owl, params.ec2go, params.pathways,
 * params.gspa_jar) are required when the flag is set; the workflow
 * fails fast in main.nf if any are missing.
 */

process BUILD_CLAIMS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container 'python:3.12-slim'

    input:
    tuple val(sample_id), path(diamond_tsv), path(pfam_domtbl), path(interproscan_tsv), path(eggnog_tsv)
    path claims_parser
    path goa_file
    path pfam2go_file

    output:
    tuple val(sample_id), path("${sample_id}_claims.jsonl"), emit: claims

    script:
    def ips_arg = (interproscan_tsv.size() > 0) ? "--interproscan ${interproscan_tsv}" : ''
    def p2g_arg = (pfam2go_file.name != 'NO_PFAM2GO') ? "--pfam2go ${pfam2go_file}" : ''
    """
    mkdir -p preds
    # Parser keys on canonical filenames in --results-dir; stage with rename.
    # Tools that didn't run upstream are zero-byte placeholders here.
    cp ${diamond_tsv}  preds/diamond_results.tsv
    cp ${pfam_domtbl}  preds/pfam_results.domtbl
    cp ${eggnog_tsv}   preds/eggnog_results.emapper.annotations

    python3 ${claims_parser} \\
        --results-dir preds \\
        --goa ${goa_file} \\
        --output ${sample_id}_claims.jsonl \\
        ${ips_arg} \\
        ${p2g_arg}
    """
}

process INTEGRATE {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container 'eclipse-temurin:21-jre'

    input:
    tuple val(sample_id), path(claims)
    path gspa_jar
    path go_owl
    path ec2go
    path pathways

    output:
    tuple val(sample_id), path("${sample_id}_integrated.tsv"),  emit: integrated
    tuple val(sample_id), path("${sample_id}_provenance.json"), emit: provenance, optional: true

    script:
    def kingdom = params.essential_profile ?: 'bacteria'
    def priors  = params.enable_priors ?: 'essentiality,coherence,gap_filling,genomic_context'
    def theta   = params.theta_file ? "--theta ${params.theta_file}" : ''
    """
    java -jar ${gspa_jar} integrate \\
        --claims ${claims} \\
        --out ${sample_id}_integrated.tsv \\
        --go-owl ${go_owl} --lite \\
        --essential-profile ${kingdom} \\
        --pathways ${pathways} \\
        --ec2go ${ec2go} \\
        --enable-priors ${priors} \\
        ${theta}
    """
}
