#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * GSPA — Genome-Scale Protein Annotation Pipeline
 *
 * Usage:
 *   nextflow run gspa-nf/main.nf \
 *     --input genome.fna \
 *     --diamond_db /path/to/uniprot_sprot.dmnd \
 *     --pfam_db /path/to/Pfam-A.hmm \
 *     -profile docker
 *
 * With database config:
 *   nextflow run gspa-nf/main.nf \
 *     -c databases.config \
 *     --input genome.fna \
 *     -profile docker
 *
 * Samplesheet CSV:
 *   nextflow run gspa-nf/main.nf --input samplesheet.csv -profile docker
 *   CSV format: sample_id,fasta
 */

include { PYRODIGAL }                                    from './modules/gene_calling'
include { DIAMOND_BLASTP; MMSEQS2_SEARCH; FOLDSEEK }     from './modules/similarity'
include { HMMSEARCH; INTERPROSCAN; EGGNOG_MAPPER; DBCAN } from './modules/domains'
include { BARRNAP; MINCED; AMRFINDER; ANTISMASH;
          SIGNALP; CHECKM2; GTDBTK }                    from './modules/specialized'
include { MERGE_ANNOTATIONS }                            from './modules/quality'
include { ESM2_DEEPGOPLUS; ESM2_CENTROID;
          CLEAN; PROTEINFER }                            from './modules/neural'
include { ENSEMBLE_PREDS }                               from './modules/ensemble'
include { EVAL_PGAP }                                    from './modules/eval'
include { MAKE_REPORT }                                  from './modules/report'

if (!params.input) {
    error "Please provide --input (FASTA file or samplesheet CSV)"
}

def create_input() {
    if (params.input.toString().endsWith('.csv')) {
        return Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample_id, file(row.fasta)) }
    } else {
        def f = file(params.input)
        return Channel.of(tuple(f.simpleName, f))
    }
}

def empty_file(name) {
    def f = file("${workDir}/${name}")
    f.text = ''
    return f
}

workflow {

    ch_input = create_input()

    // ===== Step 1: Gene calling =====
    PYRODIGAL(ch_input)

    // ===== Step 2: Sequence similarity =====
    ch_diamond = Channel.empty()
    if (params.run_diamond && params.diamond_db) {
        DIAMOND_BLASTP(PYRODIGAL.out.proteins, file(params.diamond_db))
        ch_diamond = DIAMOND_BLASTP.out.results
    }
    if (params.run_mmseqs2 && params.mmseqs2_db) {
        MMSEQS2_SEARCH(PYRODIGAL.out.proteins, file(params.mmseqs2_db))
    }

    // ===== Step 3: Structure similarity (FoldSeek) =====
    if (params.run_foldseek && params.foldseek_db) {
        // ProstT5 mode: use protein FASTA directly
        // Structure mode: would need a directory of PDB files
        def prostt5 = params.prostt5_model ? file(params.prostt5_model) : file('NO_PROSTT5')
        def query = params.prostt5_model ? PYRODIGAL.out.proteins : PYRODIGAL.out.proteins
        FOLDSEEK(query, file(params.foldseek_db), prostt5)
    }

    // ===== Step 4: Domain annotation =====
    ch_pfam = Channel.empty()
    if (params.run_hmmer && params.pfam_db) {
        def pfam = file(params.pfam_db)
        HMMSEARCH(PYRODIGAL.out.proteins, pfam,
                  file("${params.pfam_db}.h3f"), file("${params.pfam_db}.h3i"),
                  file("${params.pfam_db}.h3m"), file("${params.pfam_db}.h3p"))
        ch_pfam = HMMSEARCH.out.results
    }
    ch_interproscan = Channel.empty()
    if (params.run_interproscan && params.interproscan_dir) {
        INTERPROSCAN(PYRODIGAL.out.proteins, file(params.interproscan_dir))
        ch_interproscan = INTERPROSCAN.out.results
    }
    if (params.run_eggnog && params.eggnog_db) {
        EGGNOG_MAPPER(PYRODIGAL.out.proteins, file(params.eggnog_db))
    }
    if (params.run_dbcan && params.dbcan_db) {
        DBCAN(PYRODIGAL.out.proteins, file(params.dbcan_db))
    }

    // ===== Step 5: Genome-level annotation =====
    ch_rrna = Channel.empty()
    ch_crispr = Channel.empty()
    if (params.run_barrnap) {
        BARRNAP(PYRODIGAL.out.genome)
        ch_rrna = BARRNAP.out.results
    }
    if (params.run_minced) {
        MINCED(PYRODIGAL.out.genome)
        ch_crispr = MINCED.out.results
    }
    if (params.run_amrfinder && params.amrfinder_db) {
        AMRFINDER(PYRODIGAL.out.proteins, file(params.amrfinder_db))
    }
    if (params.run_antismash && params.antismash_db) {
        ANTISMASH(PYRODIGAL.out.genome, file(params.antismash_db))
    }
    if (params.run_signalp && params.signalp_path) {
        SIGNALP(PYRODIGAL.out.proteins)
    }

    // ===== Step 6: MAG quality (optional) =====
    if (params.run_checkm2 && params.checkm2_db) {
        CHECKM2(PYRODIGAL.out.genome, file(params.checkm2_db))
    }
    if (params.run_gtdbtk && params.gtdbtk_db) {
        GTDBTK(PYRODIGAL.out.genome, file(params.gtdbtk_db))
    }

    // ===== Step 6.5: Neural function predictors (opt-in) =====
    // Each predictor is enabled independently via params.run_<name>; their
    // outputs are NOT threaded into MERGE_ANNOTATIONS (which keeps its
    // 6-column schema). Neural TSVs land under
    //   ${outdir}/${sample_id}/{esm2_deepgoplus,esm2_centroid,clean,proteinfer}/
    // and the per-sample union is fed to ENSEMBLE_PREDS + (optionally) EVAL_PGAP.

    ch_neural = Channel.empty()
    if (params.run_esm2_deepgoplus && params.esm2_dgp_ckpt && params.esm2_dgp_terms) {
        ESM2_DEEPGOPLUS(PYRODIGAL.out.proteins,
                        file(params.esm2_dgp_ckpt),
                        file(params.esm2_dgp_terms))
        ch_neural = ch_neural.mix(ESM2_DEEPGOPLUS.out.results.map { id, f -> tuple(id, 'esm2-deepgoplus', f) })
    }
    if (params.run_esm2_centroid && params.esm2_centroid_db) {
        ESM2_CENTROID(PYRODIGAL.out.proteins, file(params.esm2_centroid_db))
        ch_neural = ch_neural.mix(ESM2_CENTROID.out.results.map { id, f -> tuple(id, 'esm2-centroid', f) })
    }
    if (params.run_clean && params.clean_model_dir) {
        CLEAN(PYRODIGAL.out.proteins, file(params.clean_model_dir))
        ch_neural = ch_neural.mix(CLEAN.out.results.map { id, f -> tuple(id, 'clean', f) })
    }
    if (params.run_proteinfer && params.proteinfer_model_dir) {
        PROTEINFER(PYRODIGAL.out.proteins, file(params.proteinfer_model_dir))
        ch_neural = ch_neural.mix(PROTEINFER.out.results.map { id, f -> tuple(id, 'proteinfer', f) })
    }

    // Ensemble: collect TSVs per sample
    ch_ensemble = Channel.empty()
    if (params.run_ensemble) {
        ch_for_ensemble = ch_neural
            .map { id, _name, f -> tuple(id, f) }
            .groupTuple()
        ENSEMBLE_PREDS(ch_for_ensemble)
        ch_ensemble = ENSEMBLE_PREDS.out.results.map { id, f -> tuple(id, 'ensemble-' + params.ensemble_mode, f) }
    }

    // Eval: each predictor TSV (and ensemble) joined to per-sample truth
    ch_eval_results = Channel.empty()
    if (params.run_eval && params.truth_dir) {
        ch_truth = Channel.fromPath("${params.truth_dir}/*_truth.tsv")
                          .map { f -> tuple(f.simpleName.replace('_truth',''), f) }
        ch_for_eval = ch_neural.mix(ch_ensemble)
            .combine(ch_truth, by: 0)
            .map { id, predictor, pred_tsv, truth_tsv -> tuple(id, predictor, pred_tsv, truth_tsv) }
        def aspect_map = params.go_obo ? file(params.go_obo) : file('NO_ASPECT_MAP')
        EVAL_PGAP(ch_for_eval, aspect_map)
        ch_eval_results = EVAL_PGAP.out.results
    }

    // Report: HTML + TTL (SIO-based RDF) + JSON-LD per sample
    // Collect (name, file) pairs across all neural + ensemble channels and
    // optional eval JSONs. Files are staged under unique names; --predictor
    // flag uses name:staged_filename so make_report.py reads each correctly.
    if (params.run_report) {
        ch_pred_pairs = ch_neural.mix(ch_ensemble)
            .map { id, name, f -> tuple(id, "${name}:${f.name}", f) }
            .groupTuple()
            .map { id, specs, files -> tuple(id, specs, files) }
        ch_eval_pairs = ch_eval_results
            .map { id, name, f -> tuple(id, "${name}:${f.name}", f) }
            .groupTuple()
            .map { id, specs, files -> tuple(id, specs, files) }
            .ifEmpty( Channel.empty() )
        // Outer-join: every sample with predictions reports, with or without eval
        ch_report = ch_pred_pairs
            .join(ch_eval_pairs, remainder: true)
            .map { id, pred_specs, pred_files, eval_specs, eval_files ->
                tuple(id, pred_specs, pred_files,
                      eval_specs ?: [], eval_files ?: file('NO_EVAL_FILE'))
            }
        MAKE_REPORT(ch_report)
    }

    // ===== Step 7: Merge annotations =====
    // Provide empty files for channels that weren't run
    ch_diamond_or_empty = params.run_diamond && params.diamond_db ?
        ch_diamond :
        PYRODIGAL.out.proteins.map { id, f -> tuple(id, empty_file("empty_diamond_${id}.tsv")) }

    ch_pfam_or_empty = params.run_hmmer && params.pfam_db ?
        ch_pfam :
        PYRODIGAL.out.proteins.map { id, f -> tuple(id, empty_file("empty_pfam_${id}.tsv")) }

    ch_rrna_or_empty = params.run_barrnap ?
        ch_rrna :
        PYRODIGAL.out.genome.map { id, f -> tuple(id, empty_file("empty_rrna_${id}.gff")) }

    ch_crispr_or_empty = params.run_minced ?
        ch_crispr :
        PYRODIGAL.out.genome.map { id, f -> tuple(id, empty_file("empty_crispr_${id}.gff")) }

    ch_interproscan_or_empty = params.run_interproscan && params.interproscan_dir ?
        ch_interproscan :
        PYRODIGAL.out.proteins.map { id, f -> tuple(id, empty_file("empty_interproscan_${id}.tsv")) }

    ch_merge = PYRODIGAL.out.gff
        .join(PYRODIGAL.out.proteins)
        .join(ch_diamond_or_empty)
        .join(ch_pfam_or_empty)
        .join(ch_interproscan_or_empty)
        .join(ch_rrna_or_empty)
        .join(ch_crispr_or_empty)

    MERGE_ANNOTATIONS(ch_merge)
}

workflow.onComplete {
    log.info """
    ========================================
    GSPA Pipeline Complete
    ========================================
    Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration : ${workflow.duration}
    Output   : ${params.outdir}
    ========================================
    """.stripIndent()
}
