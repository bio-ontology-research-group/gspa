package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.suggester.DisjunctiveSuggestion
import gspa.integration.suggester.SingletonSuggestion
import gspa.integration.suggester.Suggestion

/**
 * Export formats of the integration workflow. Every writer uses full-precision
 * posteriors; none re-derives a score. Per-genome paths are derived from the
 * genome's own output file so concurrent genomes never share an artifact.
 */
class IntegrationArtifacts {

    /** Evidence-gated KEGG module assessment before and after refinement, as {@code <out>.modules.json}. */
    static String writeModuleReport(IntegrationState state, IntegratedAnnotationSet integrated, File targetOutFile) {
        if (state?.moduleDefinitions == null) return null
        Map report = gspa.integration.ModuleReport.build(state, integrated)
        File f = new File(targetOutFile.parentFile ?: new File('.'), "${targetOutFile.name}.modules.json")
        f.text = new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString(report)
        "  Wrote: ${f} (${report.summary.modules_activated_after} activated modules, ${report.summary.modules_complete_after} complete)".toString()
    }

    /** Rescues made by the targeted priors, as {@code <out>.rescues.tsv} (only when any were made). */
    static String writeRescueLog(IntegrationState state, File targetOutFile) {
        if (state == null || state.rescueLog == null || state.rescueLog.isEmpty()) return null
        File f = new File(targetOutFile.parentFile ?: new File('.'), "${targetOutFile.name}.rescues.tsv")
        f.withWriter('UTF-8') { w ->
            w.writeLine("prior\tmode\tcategory\tprotein_id\tfunction_id\tlikelihood_prob\tscore\ttarget\tevidence")
            for (Map r : state.rescueLog) {
                String[] parts = IntegrationState.splitFunctionKey(r.functionKey as String)
                String pid = parts != null ? parts[0] : ''
                String fid = parts != null ? parts[2] : (r.functionKey as String)
                String sc = (r.score instanceof Number && !Double.isNaN((r.score as Number).doubleValue())) ? String.format(Locale.ROOT, '%.4f', (r.score as Number).doubleValue()) : ''
                String tg = (r.target instanceof Number && !Double.isNaN((r.target as Number).doubleValue())) ? String.format(Locale.ROOT, '%.2f', (r.target as Number).doubleValue()) : ''
                String ev = r.evidence instanceof Map ? (r.evidence as Map).collect { k, v -> "${k}=${v}" }.join(';') : ''
                w.writeLine("${r.prior}\t${r.mode}\t${r.category}\t${pid}\t${fid}\t${String.format(Locale.ROOT, '%.4f', (r.likelihood as Number).doubleValue())}\t${sc}\t${tg}\t${ev}")
            }
        }
        "  Wrote: ${f} (${state.rescueLog.size()} rescues)".toString()
    }

    static void writeIntegratedOutputTsv(IntegratedAnnotationSet integrated, File targetOutFile) {
        targetOutFile.parentFile?.mkdirs()
        targetOutFile.withWriter { w ->
            w.writeLine(['protein_id', 'type', 'function_id', 'go_aspect',
                         'posterior_prob', 'likelihood_logodds', 'final_logodds',
                         'n_supporting', 'priors_fired', 'convergence_iter'].join('\t'))
            integrated.provenance.values().each { prov ->
                def fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey == null) return
                def fired = prov.priorContributions.collect { k, v ->
                    "${k}:${String.format(Locale.ROOT, '%.3f', v)}"
                }.join(',')
                def firstClaim = prov.supportingClaims ? prov.supportingClaims.first() : null
                w.writeLine([
                    fkey[0],
                    fkey[1],
                    fkey[2],
                    firstClaim?.goAspect ?: '',
                    Double.toString(prov.finalProbability),
                    Double.toString(prov.likelihoodLogOdds),
                    Double.toString(prov.finalLogOdds),
                    (prov.supportingClaims?.size() ?: 0).toString(),
                    fired,
                    prov.convergenceIter.toString(),
                ].join('\t'))
            }
        }
    }

    /**
     * Write suggestions to TSV. Columns:
     *   kind, pathway_id, reaction_id, function_id, go_aspect, operon_id,
     *   bayes_factor, suggestion_score, protein_ids, q_values, provenance
     */
    static void writeSuggestionsTsv(List<Suggestion> suggestions, File out) {
        out.withWriter { w ->
            w.writeLine([
                'kind', 'pathway_id', 'reaction_id', 'function_id', 'go_aspect',
                'operon_id', 'bayes_factor', 'suggestion_score',
                'protein_ids', 'q_values', 'provenance', 'source', 'unassigned_probability', 'candidate_probabilities'
            ].join('\t'))
            for (Suggestion s : suggestions) {
                String proteinList
                String qList
                if (s instanceof SingletonSuggestion) {
                    proteinList = s.proteinId
                    qList = String.format(Locale.ROOT, '%.4f', s.q)
                } else if (s instanceof DisjunctiveSuggestion) {
                    proteinList = s.proteinIds.join(',')
                    qList = s.qValues.collect {
                        String.format(Locale.ROOT, '%.4f', it)
                    }.join(',')
                } else {
                    proteinList = ''
                    qList = ''
                }
                w.writeLine([
                    s.kind(),
                    s.pathwayId ?: '',
                    s.reactionId ?: '',
                    s.functionId ?: '',
                    s.goAspect ?: '',
                    s.operonId ?: '',
                    String.format(Locale.ROOT, '%.2f', s.bayesFactor),
                    String.format(Locale.ROOT, '%.4f', s.suggestionScore),
                    proteinList,
                    qList,
                    s.provenance ?: '',
                    s.source ?: '',
                    s.unassignedProbability != null ? String.format(Locale.ROOT, '%.6f', s.unassignedProbability) : '',
                    new ObjectMapper().writeValueAsString(s.candidateProbabilities ?: [:]),
                ].join('\t'))
            }
        }
    }

    /** Per-(protein, function) provenance, as JSON. */
    static void writeProvenanceJson(IntegratedAnnotationSet integrated, File target) {
        def provList = integrated.provenance.values().collect { prov ->
            [
                function_key: prov.functionKey,
                protein_id: prov.proteinId,
                likelihood_logodds: prov.likelihoodLogOdds,
                final_logodds: prov.finalLogOdds,
                final_probability: prov.finalProbability,
                prior_contributions: prov.priorContributions,
                n_supporting: prov.supportingClaims?.size() ?: 0,
                supporting_sources: prov.supportingClaims ? prov.supportingClaims*.source : [],
                convergence_iter: prov.convergenceIter,
            ]
        }
        target.parentFile?.mkdirs()
        new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(target, provList)
    }
}
