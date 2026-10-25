package gspa.integration.fba

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimProvenance
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.model.AnnotationType

/**
 * Stage 2B: StoichiometricLocusRefiner — Positive-Cost Flux-Certified Metabolic Repair,
 * Lexicographic pFBA, FVA/Knockout/Ensemble Certificates, and Backward Locus Responsibility r_j(h).
 */
class StoichiometricLocusRefiner {

    File fbaDataDir
    File sidecarScript = new File('scripts/stoichiometric_locus_refiner.py')
    String pythonBin = 'python3'
    String gram = 'neg'
    String medium = 'minimal'
    double tauCore = 0.50d
    double alphaBio = 0.10d
    int topM = 3
    double tauSub = 0.25d

    SparseRoleAggregator aggregator = new SparseRoleAggregator()
    LocusResponsibilityUpdater updater = new LocusResponsibilityUpdater()

    static class RefinementSummary {
        int nCoreRxns = 0
        int nRepairedRxns = 0
        int nUncertifiedAdditions = 0
        int nUpdatedLoci = 0
        double coreBio = 0.0d
        double maxUnivBio = 0.0d
        double solveTimeSec = 0.0d
        List<String> repairedRxns = []
        List<String> repairedEcs = []
    }

    /**
     * Extract current gene -> (EC -> posterior_prob) map from IntegratedAnnotationSet.
     */
    Map<String, Map<String, Double>> extractGeneEcScores(IntegratedAnnotationSet integrated) {
        Map<String, Map<String, Double>> geneEcScores = new LinkedHashMap<String, Map<String, Double>>()
        if (integrated == null) return geneEcScores
        for (ClaimProvenance prov : integrated.provenance.values()) {
            String[] parts = IntegrationState.splitFunctionKey(prov.functionKey)
            if (parts == null || parts.length < 3) continue
            if (parts[1] != 'EC') continue
            String pid = parts[0]
            String ec = parts[2].replace('EC:', '').trim()
            if (ec.isEmpty() || ec.contains('-')) continue
            Map<String, Double> ecMap = geneEcScores.get(pid)
            if (ecMap == null) {
                ecMap = new LinkedHashMap<String, Double>()
                geneEcScores.put(pid, ecMap)
            }
            ecMap.put(ec, prov.finalProbability)
        }
        return geneEcScores
    }

    /**
     * Extract contextual covariates (RLGC Bayes factors, operon/synteny boosts, and PLM similarities)
     * keyed by "geneId|ec" for backward locus responsibility attribution.
     */
    Map<String, Map<String, Double>> extractContextualCovariates(IntegratedAnnotationSet integrated) {
        Map<String, Double> rlgcBf = new LinkedHashMap<String, Double>()
        Map<String, Double> operonBoost = new LinkedHashMap<String, Double>()
        Map<String, Double> plmSim = new LinkedHashMap<String, Double>()
        if (integrated != null) {
            for (ClaimProvenance prov : integrated.provenance.values()) {
                String[] parts = IntegrationState.splitFunctionKey(prov.functionKey)
                if (parts == null || parts.length < 3 || parts[1] != 'EC') continue
                String pid = parts[0]
                String ec = parts[2].replace('EC:', '').trim()
                if (ec.isEmpty() || ec.contains('-')) continue
                String key = "${pid}|${ec}".toString()
                double gapPrior = (prov.priorContributions?.get('gap_filling') ?: 0.0d) as double
                if (gapPrior != 0.0d) {
                    rlgcBf.put(key, Math.exp(gapPrior / 1.5d))
                }
                double ctxPrior = ((prov.priorContributions?.get('genomic_context') ?: 0.0d) as double) +
                    ((prov.priorContributions?.get('coherence') ?: 0.0d) as double)
                if (ctxPrior > 0.0d) {
                    operonBoost.put(key, Math.min(1.0d, Math.max(0.0d, ctxPrior)))
                }
                double bestPlm = 0.0d
                for (def sc : (prov.supportingClaims ?: [])) {
                    String src = (sc?.source ?: '').toLowerCase(Locale.ROOT)
                    if (src.contains('esm') || src.contains('prott5') || src.contains('plm') || src.contains('clean')) {
                        bestPlm = Math.max(bestPlm, (sc.rawScore ?: 0.0d) as double)
                    }
                }
                if (bestPlm > 0.0d) {
                    plmSim.put(key, Math.min(1.0d, Math.max(0.0d, bestPlm)))
                }
            }
        }
        return [
            rlgc_bf_by_gene_ec             : rlgcBf,
            operon_partner_boost_by_gene_ec: operonBoost,
            plm_sim_by_gene_ec             : plmSim
        ]
    }

    /**
     * Run StoichiometricLocusRefiner on an IntegratedAnnotationSet and update posteriors via r_j(h).
     */
    RefinementSummary refine(IntegrationState state, IntegratedAnnotationSet integrated, File workDir = null) {
        RefinementSummary summary = new RefinementSummary()
        Map<String, Map<String, Double>> geneEcScores = extractGeneEcScores(integrated)
        if (geneEcScores.isEmpty() || fbaDataDir == null || !fbaDataDir.exists()) {
            return summary
        }

        aggregator.topM = topM
        aggregator.tauSub = tauSub
        Map<String, Double> roleProbs = aggregator.aggregateRoleProbabilities(geneEcScores)

        File tmpDir = workDir ?: File.createTempDir('gspa_fba_', '')
        tmpDir.mkdirs()
        File inJson = new File(tmpDir, 'fba_input.json')
        File outJson = new File(tmpDir, 'fba_output.json')

        Map<String, Map<String, Double>> covariates = extractContextualCovariates(integrated)
        ObjectMapper mapper = new ObjectMapper()
        Map<String, Object> payload = [
            gene_ec_scores                 : geneEcScores,
            role_probs                     : roleProbs,
            rlgc_bf_by_gene_ec             : covariates.rlgc_bf_by_gene_ec,
            operon_partner_boost_by_gene_ec: covariates.operon_partner_boost_by_gene_ec,
            plm_sim_by_gene_ec             : covariates.plm_sim_by_gene_ec,
            gram                           : gram,
            medium                         : medium,
            tau_core                       : tauCore,
            alpha_bio                      : alphaBio,
            top_m                          : topM,
            tau_sub                        : tauSub
        ]
        mapper.writeValue(inJson, payload)

        List<String> cmd = [
            pythonBin,
            sidecarScript.absolutePath,
            '--fba-data-dir', fbaDataDir.absolutePath,
            '--input-json', inJson.absolutePath,
            '--output-json', outJson.absolutePath,
            '--export-model', new File(tmpDir, 'metabolic_model.xml').absolutePath,
            '--gram', gram,
            '--medium', medium,
            '--tau-core', Double.toString(tauCore),
            '--alpha-bio', Double.toString(alphaBio)
        ]

        ProcessBuilder pb = new ProcessBuilder(cmd)
        pb.redirectErrorStream(true)
        File fbaLog = new File(tmpDir, 'fba.log')
        pb.redirectOutput(fbaLog)
        Process proc = pb.start()
        int exitCode = proc.waitFor()
        String procOut = fbaLog.text
        if (exitCode != 0 || !outJson.exists()) {
            throw new IllegalStateException("StoichiometricLocusRefiner failed (${exitCode}): ${procOut}")
        }

        Map result = mapper.readValue(outJson, Map)
        summary.nCoreRxns = (result.n_core_rxns as Number ?: 0).intValue()
        summary.nRepairedRxns = (result.n_repaired_rxns as Number ?: 0).intValue()
        summary.nUncertifiedAdditions = (result.n_uncertified_additions as Number ?: 0).intValue()
        summary.coreBio = (result.core_bio as Number ?: 0.0d).doubleValue()
        summary.maxUnivBio = (result.max_univ_bio as Number ?: 0.0d).doubleValue()
        summary.solveTimeSec = (result.solve_time_sec as Number ?: 0.0d).doubleValue()
        summary.repairedRxns = (result.repaired_rxns as List<String>) ?: []
        summary.repairedEcs = (result.repaired_ecs as List<String>) ?: []

        Map<String, Map<String, Double>> updatedScores = (Map<String, Map<String, Double>>) result.updated_gene_ec_scores
        if (updatedScores != null) {
            summary.nUpdatedLoci = updater.applyUpdatedScores(integrated, updatedScores)
        }
        return summary
    }
}
