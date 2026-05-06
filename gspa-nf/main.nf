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
include { BUILD_CLAIMS; INTEGRATE }                      from './modules/integrate'
include { ESM2_DEEPGOPLUS; ESM2_CENTROID;
          CLEAN; PROTEINFER }                            from './modules/neural'
include { ENSEMBLE_PREDS; ENSEMBLE_GENOMIC }              from './modules/ensemble'
include { EVAL_PGAP }                                    from './modules/eval'
include { MAKE_REPORT }                                  from './modules/report'
include { METAPREDICT; DEEPSIG; TMBED; TPPRED3 }         from './modules/region'
include { PSORTB }                                       from './modules/loc'
include { DEEPFRI; DEEPEC; DEEPARG }                     from './modules/term_extras'
include { MUSITEDEEP; SCANNET }                          from './modules/sites'
include { STRUCTURE_PROVIDER }                           from './modules/structure'
include { GENOMAD; CHECKV; PHISPY; VIRSORTER2; VIBRANT; FLATTEN_PROPHAGE_CDS } from './modules/viral'

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
    ch_eggnog = Channel.empty()
    if (params.run_eggnog && params.eggnog_db) {
        EGGNOG_MAPPER(PYRODIGAL.out.proteins, file(params.eggnog_db))
        ch_eggnog = EGGNOG_MAPPER.out.results
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

    // ===== Step 6.6: FOSS region predictors (opt-in) =====
    ch_regions = Channel.empty()
    if (params.run_metapredict) {
        METAPREDICT(PYRODIGAL.out.proteins)
        ch_regions = ch_regions.mix(METAPREDICT.out.results.map { id, f -> tuple(id, 'metapredict', f) })
    }
    if (params.run_deepsig) {
        DEEPSIG(PYRODIGAL.out.proteins)
        ch_regions = ch_regions.mix(DEEPSIG.out.results.map { id, f -> tuple(id, 'deepsig', f) })
    }
    if (params.run_tmbed) {
        def prott5 = params.tmbed_prott5_dir ? file(params.tmbed_prott5_dir) : file('NO_PROTT5')
        TMBED(PYRODIGAL.out.proteins, prott5)
        ch_regions = ch_regions.mix(TMBED.out.results.map { id, f -> tuple(id, 'tmbed', f) })
    }
    if (params.run_tppred3) {
        TPPRED3(PYRODIGAL.out.proteins)
        ch_regions = ch_regions.mix(TPPRED3.out.results.map { id, f -> tuple(id, 'tppred3', f) })
    }

    // ===== Step 6.7: FOSS term-extra predictors (auto-join ensemble) =====
    if (params.run_psortb) {
        PSORTB(PYRODIGAL.out.proteins)
        ch_neural = ch_neural.mix(PSORTB.out.results.map { id, f -> tuple(id, 'psortb', f) })
    }
    if (params.run_deepfri && params.deepfri_model_dir) {
        DEEPFRI(PYRODIGAL.out.proteins, file(params.deepfri_model_dir))
        ch_neural = ch_neural.mix(DEEPFRI.out.results.map { id, f -> tuple(id, 'deepfri', f) })
    }
    if (params.run_deepec && params.deepec_model_dir) {
        DEEPEC(PYRODIGAL.out.proteins, file(params.deepec_model_dir))
        ch_neural = ch_neural.mix(DEEPEC.out.results.map { id, f -> tuple(id, 'deepec', f) })
    }
    if (params.run_deeparg && params.deeparg_model_dir) {
        DEEPARG(PYRODIGAL.out.proteins, file(params.deeparg_model_dir))
        ch_neural = ch_neural.mix(DEEPARG.out.results.map { id, f -> tuple(id, 'deeparg', f) })
    }

    // ===== Step 6.8: Optional structure provisioning + site predictors =====
    ch_sites = Channel.empty()
    ch_structures = Channel.empty()
    if (params.structures_from && params.structures_from != 'none') {
        STRUCTURE_PROVIDER(PYRODIGAL.out.proteins)
        ch_structures = STRUCTURE_PROVIDER.out.structures
    }
    if (params.run_musitedeep && params.musitedeep_model_dir) {
        MUSITEDEEP(PYRODIGAL.out.proteins, file(params.musitedeep_model_dir))
        ch_sites = ch_sites.mix(MUSITEDEEP.out.results.map { id, f -> tuple(id, 'musitedeep', f) })
    }
    if (params.run_scannet && params.scannet_model_dir &&
        params.structures_from && params.structures_from != 'none') {
        ch_scannet_in = PYRODIGAL.out.proteins.join(ch_structures, by: 0)
                          .map { id, prot, struct -> tuple(id, prot) }
        SCANNET(ch_scannet_in, file(params.scannet_model_dir), ch_structures.map { id, d -> d })
        ch_sites = ch_sites.mix(SCANNET.out.results.map { id, f -> tuple(id, 'scannet', f) })
    }

    // ===== Step 6.9: v1.3 Phage / prophage track (genomic-level, opt-in) =====
    ch_genomic = Channel.empty()
    if (params.run_genomad && params.genomad_db) {
        GENOMAD(PYRODIGAL.out.genome, file(params.genomad_db))
        ch_genomic = ch_genomic.mix(GENOMAD.out.results.map { id, f -> tuple(id, 'genomad', f) })
    }
    if (params.run_checkv && params.checkv_db) {
        CHECKV(PYRODIGAL.out.genome, file(params.checkv_db))
        ch_genomic = ch_genomic.mix(CHECKV.out.results.map { id, f -> tuple(id, 'checkv', f) })
    }
    // PhiSpy needs a GenBank file (caller must supply it as params.gbk_input
    // or skip the predictor); kept as opt-in
    if (params.run_phispy && params.phispy_input) {
        ch_phispy_in = Channel.fromPath(params.phispy_input)
                              .map { f -> tuple(f.simpleName, f) }
        PHISPY(ch_phispy_in)
        ch_genomic = ch_genomic.mix(PHISPY.out.results.map { id, f -> tuple(id, 'phispy', f) })
    }
    if (params.run_virsorter2 && params.virsorter2_db) {
        VIRSORTER2(PYRODIGAL.out.genome, file(params.virsorter2_db))
        ch_genomic = ch_genomic.mix(VIRSORTER2.out.results.map { id, f -> tuple(id, 'virsorter2', f) })
    }
    if (params.run_vibrant && params.vibrant_db) {
        VIBRANT(PYRODIGAL.out.genome, file(params.vibrant_db))
        ch_genomic = ch_genomic.mix(VIBRANT.out.results.map { id, f -> tuple(id, 'vibrant', f) })
    }
    // v1.4 genomic-region ensemble (interval-overlap fusion)
    if (params.run_ensemble_genomic) {
        ch_for_genomic_ensemble = ch_genomic
            .map { id, name, f -> tuple(id, "${name}:${f.name}", f) }
            .groupTuple()
            .map { id, specs, files -> tuple(id, specs, files) }
        ENSEMBLE_GENOMIC(ch_for_genomic_ensemble)
        ch_genomic = ch_genomic.mix(
            ENSEMBLE_GENOMIC.out.results.map { id, f ->
                tuple(id, 'ensemble-genomic-' + params.genomic_ensemble_mode, f)
            }
        )
    }
    // Flatten genomic regions to per-CDS region annotations (joins ch_regions)
    if (params.run_flatten_prophage_cds) {
        ch_flat_in = ch_genomic
            .map { id, name, f -> tuple(id, f) }
            .combine(PYRODIGAL.out.gff, by: 0)
        FLATTEN_PROPHAGE_CDS(ch_flat_in)
        ch_regions = ch_regions.mix(FLATTEN_PROPHAGE_CDS.out.results.map { id, f -> tuple(id, 'prophage_cds', f) })
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
    // Collects (name, file) tuples across all neural + ensemble channels,
    // region channels, site channels, genomic-region channels, and
    // optional eval JSONs.
    if (params.run_report) {
        ch_pred_pairs = ch_neural.mix(ch_ensemble)
            .map { id, name, f -> tuple(id, "${name}:${f.name}", f) }
            .groupTuple()
            .map { id, specs, files -> tuple(id, specs, files) }
        ch_region_pairs = ch_regions
            .map { id, name, f -> tuple(id, "${name}:${f.name}", f) }
            .groupTuple()
            .map { id, specs, files -> tuple(id, specs, files) }
        ch_site_pairs = ch_sites
            .map { id, name, f -> tuple(id, "${name}:${f.name}", f) }
            .groupTuple()
            .map { id, specs, files -> tuple(id, specs, files) }
        ch_genomic_pairs = ch_genomic
            .map { id, name, f -> tuple(id, "${name}:${f.name}", f) }
            .groupTuple()
            .map { id, specs, files -> tuple(id, specs, files) }
        ch_eval_pairs = ch_eval_results
            .map { id, name, f -> tuple(id, "${name}:${f.name}", f) }
            .groupTuple()
            .map { id, specs, files -> tuple(id, specs, files) }
        // Outer-join everything so every sample with predictions reports
        ch_report = ch_pred_pairs
            .join(ch_region_pairs,  remainder: true)
            .join(ch_site_pairs,    remainder: true)
            .join(ch_genomic_pairs, remainder: true)
            .join(ch_eval_pairs,    remainder: true)
            .map { id, p_specs, p_files, r_specs, r_files, s_specs, s_files, g_specs, g_files, e_specs, e_files ->
                tuple(id,
                      p_specs, p_files,
                      r_specs ?: [], r_files ?: file('NO_REGION_FILE'),
                      s_specs ?: [], s_files ?: file('NO_SITE_FILE'),
                      g_specs ?: [], g_files ?: file('NO_GENOMIC_FILE'),
                      e_specs ?: [], e_files ?: file('NO_EVAL_FILE'))
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

    // ===== Step 8: Optional Phase 7 integration =====
    // Lifts per-tool TSVs into a single claims.jsonl, then runs the
    // JVM `gspa-cli integrate` to produce per-(protein, function)
    // posterior probabilities with the full prior stack.
    if (params.run_integrate) {
        if (!params.gspa_jar) error 'run_integrate=true requires --gspa_jar (path to gspa-cli shadowJar)'
        if (!params.goa)      error 'run_integrate=true requires --goa (Swiss-Prot/GOA annotation file, .gaf or .gaf.gz)'
        if (!params.go_owl)   error 'run_integrate=true requires --go_owl (GO ontology OWL file)'
        if (!params.ec2go)    error 'run_integrate=true requires --ec2go (EC -> GO mapping file)'
        if (!params.pathways) error 'run_integrate=true requires --pathways (KEGG pathway TSV)'

        ch_eggnog_or_empty = params.run_eggnog && params.eggnog_db ?
            ch_eggnog :
            PYRODIGAL.out.proteins.map { id, f -> tuple(id, empty_file("empty_eggnog_${id}.tsv")) }

        ch_claims_in = ch_diamond_or_empty
            .join(ch_pfam_or_empty)
            .join(ch_interproscan_or_empty)
            .join(ch_eggnog_or_empty)

        BUILD_CLAIMS(
            ch_claims_in,
            file("${projectDir}/../benchmark/02b_parse_predictors_to_claims.py"),
            file(params.goa),
            params.pfam2go ? file(params.pfam2go) : file('NO_PFAM2GO')
        )

        INTEGRATE(
            BUILD_CLAIMS.out.claims,
            file(params.gspa_jar),
            file(params.go_owl),
            file(params.ec2go),
            file(params.pathways)
        )
    }
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
