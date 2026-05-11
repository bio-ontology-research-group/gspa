/*
 * Phase 7 evidence integrator + downstream artifacts for the Nextflow path.
 *
 * Closes the gap noted in gspa-nf/README.md: the JVM CLI's `gspa
 * integrate` produces per-(protein, function) posterior probabilities
 * with full prior stack and provenance, while the Nextflow pipeline
 * historically stopped at the merged TSV. Six processes:
 *
 *   BUILD_CLAIMS    — wraps benchmark/02b_parse_predictors_to_claims.py
 *                     to lift the per-tool TSVs into a single claims.jsonl
 *   SIDECAR_CLAIMS  — wraps the bundled sidecar_to_claims.py to emit a
 *                     parallel claims.jsonl from the sidecar tools
 *                     (mDeepFRI, ProteInfer, CLEAN). Together with
 *                     BUILD_CLAIMS this delivers all 8 source tools to
 *                     the integrator (vs the historical 4 that silently
 *                     dropped the rest — see ClaimExtractor.SOURCE_TO_TYPE
 *                     fix in v1.5).
 *   MERGE_CLAIMS    — concatenates builtin + sidecar claims.
 *   OPERONS         — runs the bundled predict_operons.py (3-predictor
 *                     ensemble with Noisy-OR per-pair posterior).
 *   INTEGRATE       — invokes `gspa-cli integrate` on the merged claims
 *                     plus reference data + operons; emits
 *                     ${sample_id}_integrated.tsv + provenance + GAEF.
 *   VISUALIZE       — runs `gspa visualize` to emit the self-contained
 *                     HTML browser (CDS, operons, BGCs, AMR, igv.js
 *                     genome view, GAEF detail with names, KEGG pathway
 *                     coverage).
 *
 * All opt-in via params.run_integrate. Reference-data paths
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

/*
 * Sidecar tool claims — mDeepFRI / ProteInfer / CLEAN. These tools are not
 * handled by 02b_parse_predictors_to_claims.py because they emit per-protein
 * (term, score) TSVs in the same shape, so a separate parser exists. The
 * v1.5 fix to ClaimExtractor lets these source names actually contribute to
 * the integrator (before, they were silently dropped at the EvidenceType
 * lookup).
 */
process SIDECAR_CLAIMS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container 'python:3.12-slim'

    input:
    tuple val(sample_id), path(mdf_tsv), path(proteinfer_tsv), path(clean_tsv)
    path sidecar_parser
    path go_obo

    output:
    tuple val(sample_id), path("${sample_id}_claims_sidecar.jsonl"), emit: claims

    script:
    def args = []
    if (mdf_tsv.size()       > 0) args << "--input mdf ${mdf_tsv}"
    if (proteinfer_tsv.size() > 0) args << "--input proteinfer ${proteinfer_tsv}"
    if (clean_tsv.size()     > 0) args << "--input clean ${clean_tsv}"
    def aspect = (go_obo.name != 'NO_GO_OBO') ? "--go-aspect-map ${go_obo}" : ''
    if (args.isEmpty()) {
        // Nothing to parse — emit empty file so the merge downstream is well-defined.
        """
        : > ${sample_id}_claims_sidecar.jsonl
        """
    } else {
        """
        python3 ${sidecar_parser} \\
            ${args.join(' ')} \\
            --output ${sample_id}_claims_sidecar.jsonl \\
            ${aspect}
        """
    }
}

/*
 * Combine builtin + sidecar claims so the integrator sees every tool. The
 * order is irrelevant — the integrator groups by (protein, function) and
 * applies Noisy-OR with per-EvidenceType correlation collapse.
 */
process MERGE_CLAIMS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container 'python:3.12-slim'

    input:
    tuple val(sample_id), path(builtin_claims), path(sidecar_claims)

    output:
    tuple val(sample_id), path("${sample_id}_claims.jsonl"), emit: claims

    script:
    """
    cat ${builtin_claims} ${sidecar_claims} > ${sample_id}_claims.jsonl
    """
}

/*
 * Operon prediction. Default: 3-predictor ensemble (distance + strict +
 * functional) with Noisy-OR per-pair posteriors; falls back to plain
 * distance when no GAF is available yet (first-pass run before integrate).
 *
 * Emits operons_for_integrate.tsv (the format `gspa integrate --operons`
 * consumes) plus operons.tsv + protein_to_operon.tsv for the visualizer.
 */
process OPERONS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container 'python:3.12-slim'

    input:
    tuple val(sample_id), path(gff)
    path operon_script

    output:
    tuple val(sample_id),
          path("${sample_id}_operons.tsv"),
          path("${sample_id}_protein_to_operon.tsv"),
          path("${sample_id}_operons_for_integrate.tsv"),
          emit: operons

    script:
    """
    mkdir -p out
    python3 ${operon_script} \\
        --gff ${gff} \\
        --out-dir out
    mv out/operons.tsv               ${sample_id}_operons.tsv
    mv out/protein_to_operon.tsv     ${sample_id}_protein_to_operon.tsv
    mv out/operons_for_integrate.tsv ${sample_id}_operons_for_integrate.tsv
    """
}

process INTEGRATE {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container 'eclipse-temurin:21-jre'

    input:
    tuple val(sample_id), path(claims), path(operons_tsv)
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
    def operons = (operons_tsv.size() > 0) ? "--operons ${operons_tsv}" : ''
    """
    java -jar ${gspa_jar} integrate \\
        --claims ${claims} \\
        --out ${sample_id}_integrated.tsv \\
        --provenance ${sample_id}_provenance.json \\
        --go-owl ${go_owl} --lite \\
        --essential-profile ${kingdom} \\
        --pathways ${pathways} \\
        --ec2go ${ec2go} \\
        --enable-priors ${priors} \\
        ${operons} \\
        ${theta}
    """
}

/*
 * Self-contained HTML browser for the per-sample workspace. Runs the
 * `gspa visualize` subcommand (which extracts a bundled Python templater
 * and invokes python3). Output is the same single HTML the tutorial
 * script builds, so anyone can scp it back and double-click it.
 */
process VISUALIZE {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container 'leechuck/gspa-cli:1.5.2'  // java 21 + python3 + bundled visualize templater

    input:
    tuple val(sample_id),
          path(prokka_dir, stageAs: 'prokka_out'),
          path(integrated_tsv),
          path(provenance_json),
          path(operons_tsv),
          path(quality_json),
          path(fasta)
    path gspa_jar
    path go_obo
    path ec2go

    output:
    tuple val(sample_id), path("${sample_id}_browser.html"), emit: html

    script:
    """
    # Lay out the canonical workspace shape that make_viz.py expects.
    mkdir -p workspace/prokka_out workspace/gspa_out workspace/input
    cp -r prokka_out/* workspace/prokka_out/  || true
    cp ${integrated_tsv}   workspace/gspa_out/integrated.tsv
    cp ${provenance_json}  workspace/gspa_out/provenance.json
    cp ${operons_tsv}      workspace/gspa_out/operons.tsv
    cp ${quality_json}     workspace/gspa_out/quality_gspa.json
    cp ${fasta}            workspace/input/${sample_id}_assembly.fa
    java -jar ${gspa_jar} visualize \\
        --workdir workspace \\
        --run-dir gspa_out \\
        --genome-id ${sample_id} \\
        --go-obo ${go_obo} \\
        --ec2go  ${ec2go} \\
        --fasta  workspace/input/${sample_id}_assembly.fa \\
        --out    ${sample_id}_browser.html
    """
}
