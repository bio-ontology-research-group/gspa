package gspa.integration.workflow

import gspa.integration.CalibrationTable
import gspa.integration.EvidenceClaim
import gspa.integration.EvidenceType
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.metrics.ConsistencyEnforcer
import gspa.metrics.EnforcementAction
import gspa.metrics.EnforcementReport
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein

/**
 * Stage 2A symbolic constraint and coherence enforcement over an integrated
 * annotation set. Configuration names resolve on the workflow config {@code c}
 * (see {@link IntegrationReferenceLoader}).
 */
class SymbolicEnforcementStage {
    final IntegrationWorkflowConfig c
    Closure<Void> info = { String m -> }
    ConsistencyEnforcer.Result lastResult

    SymbolicEnforcementStage(IntegrationWorkflowConfig config) {
        this.c = config
    }

    def propertyMissing(String name) { c.getProperty(name) }

    private void println(Object message) { info.call(String.valueOf(message)) }

    void run(IntegrationState state, IntegratedAnnotationSet integrated) {
        runStage2Enforcement(state, integrated)
    }

    /**
     * Run Stage 2 Symbolic Constraint & Coherence Enforcement on the Stage 1
     * IntegratedAnnotationSet, updating both {@code integrated.annotations} and
     * {@code integrated.provenance} so downstream CAFA and GAEF evaluators
     * observe the enforced/relaxed/imputed annotations.
     */
    private void runStage2Enforcement(IntegrationState state, IntegratedAnnotationSet integrated,
                                      String modeOverride = null, File actionsOutOverride = null,
                                      boolean verbose = true) {
        if (state.satConsistencyChecker == null) {
            if (verbose) println "  [warn] Stage 2 enforcement skipped: no SAT consistency checker available"
            return
        }
        String effMode = modeOverride ?: (enforceConsistencyMode ?: 'relax')
        File effActionsOut = actionsOutOverride != null ? actionsOutOverride : (verbose ? enforcementActionsOut : null)
        Genome genome = new Genome(id: 'stage2')
        Contig contig = new Contig(id: 'contig_1', genome: genome)
        genome.contigs << contig
        Map<String, Protein> byProtein = new LinkedHashMap<>()

        // Populate genome with all GO and EC claims in integrated.provenance
        integrated.provenance.values().each { prov ->
            String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
            if (fkey == null || !(fkey[1] in ['GO', 'EC'])) return
            String pid = fkey[0]
            String termId = fkey[2]
            AnnotationType annType = fkey[1] == 'EC' ? AnnotationType.EC : AnnotationType.GO
            Protein p = byProtein.computeIfAbsent(pid) {
                def np = new Protein(id: pid, contig: contig)
                contig.proteins << np
                np
            }
            def firstClaim = prov.supportingClaims ? prov.supportingClaims.first() : null
            p.annotations.add(new Annotation(
                type: annType,
                value: termId,
                score: prov.finalProbability,
                source: firstClaim?.source ?: 'integrator',
                evidence: 'IEA',
                goAspect: firstClaim?.goAspect
            ))
        }

        Set<String> enabledStage2Set = (enabledPriors ?: '').split(',').collect { it.trim() }.findAll { it && it != 'none' } as Set
        boolean useAnatomicalInStage2 = enabledStage2Set.isEmpty() || enabledStage2Set.contains('all') || enabledStage2Set.contains('anatomical')
        def report = new EnforcementReport(provenance: true)
        def consEnf = new ConsistencyEnforcer(
            goOntology: state.goOntology,
            satChecker: state.satConsistencyChecker,
            anatomicalConstraints: useAnatomicalInStage2 ? state.anatomicalConstraints : null,
            ecToGoAlternatives: state.ecToGoAlternatives ?: state.satConsistencyChecker?.ecToGoAlternatives ?: [:],
            mode: effMode,
            report: report
        )
        def consRes = consEnf.enforce(genome)
        lastResult = consRes
        if (verbose) println "  Stage 2 Consistency (${effMode}): ${consRes}"
        if (!consRes.constraintsApplied) {
            return
        }

        // Rebuild GO entries in integrated.provenance directly from the enforced genome proteins
        // so removals, downranks, and ancestor relaxations are 100% reflected.
        int relaxedAdded = 0
        Set<String> goKeysToRemove = integrated.provenance.keySet().findAll { k ->
            String[] fkey = IntegrationState.splitFunctionKey(k)
            fkey != null && fkey[1] == 'GO'
        } as Set
        Map<String, gspa.integration.ClaimProvenance> oldGoProv = [:]
        for (String k : goKeysToRemove) {
            oldGoProv[k] = integrated.provenance.remove(k)
        }
        integrated.annotations.annotations.removeAll { it.type == AnnotationType.GO }

        String orgTax = state.satConsistencyChecker?.organismTaxon
        boolean isProkaryote = orgTax != null && (
            orgTax in ['NCBITaxon_2', 'NCBITaxon_2157'] ||
            state.satConsistencyChecker.ancestorsWithSelf(orgTax).any { it in ['NCBITaxon_2', 'NCBITaxon_2157'] }
        )

        for (Protein p : contig.proteins) {
            for (Annotation a : p.annotations.goAnnotations()) {
                String fkey = IntegrationState.splitFunctionKey(IntegrationState.functionKey(p.id, 'GO', a.value)) ? IntegrationState.functionKey(p.id, 'GO', a.value) : "${p.id}|GO|${a.value}"
                gspa.integration.ClaimProvenance origProv = null
                double bestUnpenalizedProb = -1.0d
                double maxOrigFinalProb = -1.0d
                int maxOrigNsup = 0
                boolean isRelaxed = a.evidence?.contains('relaxed-from:')
                if (isRelaxed) {
                    def m = (a.evidence =~ /relaxed-from:(GO:\d+)/)
                    while (m.find()) {
                        String fromGo = m.group(1)
                        def candProv = oldGoProv[IntegrationState.functionKey(p.id, 'GO', fromGo)]
                        if (candProv != null) {
                            if (candProv.finalProbability > maxOrigFinalProb) {
                                maxOrigFinalProb = candProv.finalProbability
                            }
                            int nsup = candProv.supportingClaims?.size() ?: 0
                            if (nsup > maxOrigNsup) maxOrigNsup = nsup
                            double validPriorSum = (candProv.priorContributions?.findAll { k, v ->
                                !(k in ['consistency', 'anatomical']) || (v != null && v > 0.0d)
                            }?.values()?.sum() ?: 0.0d) as double
                            double unpenalizedLo = candProv.likelihoodLogOdds + validPriorSum
                            double unpenalizedProb = CalibrationTable.logOddsToProb(unpenalizedLo) * consEnf.relaxFactor
                            if (unpenalizedProb > bestUnpenalizedProb) {
                                bestUnpenalizedProb = unpenalizedProb
                                origProv = candProv
                            }
                        }
                    }
                }
                def existing = oldGoProv[fkey]
                if (existing != null) {
                    if (isRelaxed) {
                        if (isProkaryote) {
                            int exNsup = existing.supportingClaims?.size() ?: 0
                            if (exNsup >= 2 && maxOrigNsup >= 2) {
                                existing.finalLogOdds = Math.min(8.0d, existing.finalLogOdds + 0.28d)
                                existing.finalProbability = CalibrationTable.logOddsToProb(existing.finalLogOdds)
                            } else if (maxOrigFinalProb >= 0.50d) {
                                existing.finalProbability = Math.max(existing.finalProbability, Math.max(0.51d, maxOrigFinalProb * consEnf.relaxFactor))
                                existing.finalLogOdds = CalibrationTable.probToLogOdds(existing.finalProbability)
                            }
                        } else {
                            if (bestUnpenalizedProb > a.score) a.score = bestUnpenalizedProb
                            existing.finalProbability = Math.max(existing.finalProbability, a.score)
                            existing.finalLogOdds = CalibrationTable.probToLogOdds(existing.finalProbability)
                        }
                        a.score = existing.finalProbability
                        existing.priorContributions['stage2_relax'] = 1.0d
                    } else {
                        existing.finalProbability = a.score
                        existing.finalLogOdds = CalibrationTable.probToLogOdds(existing.finalProbability)
                    }
                    integrated.provenance[fkey] = existing
                } else {
                    if (isProkaryote) {
                        if (maxOrigFinalProb < 0.50d) {
                            continue
                        }
                        a.score = Math.max(0.51d, maxOrigFinalProb * consEnf.relaxFactor)
                    } else if (bestUnpenalizedProb > a.score) {
                        a.score = bestUnpenalizedProb
                    }
                    def dummyClaim = new EvidenceClaim(
                        proteinId: p.id,
                        functionId: a.value,
                        functionType: AnnotationType.GO,
                        goAspect: a.goAspect ?: (origProv?.supportingClaims ? origProv.supportingClaims.first().goAspect : (state.goOntology?.getAspect(a.value) ?: '')),
                        source: a.source ?: 'stage2_relax',
                        evidenceType: EvidenceType.SEQUENCE_SIMILARITY,
                        rawScore: a.score,
                        calibratedProb: a.score
                    )
                    double lo = CalibrationTable.probToLogOdds(a.score)
                    Map<String, Double> priors = [:]
                    if (origProv?.priorContributions) {
                        origProv.priorContributions.each { k, v ->
                            if (!(k in ['consistency', 'anatomical']) || (v != null && v > 0.0d)) {
                                priors[k] = v
                            }
                        }
                    }
                    priors['stage2_relax'] = 1.0d
                    integrated.provenance[fkey] = new gspa.integration.ClaimProvenance(
                        functionKey: fkey,
                        proteinId: p.id,
                        supportingClaims: origProv?.supportingClaims ?: [dummyClaim],
                        likelihoodLogOdds: origProv ? origProv.likelihoodLogOdds : lo,
                        priorContributions: priors,
                        finalLogOdds: lo,
                        finalProbability: a.score,
                        convergenceIter: origProv ? origProv.convergenceIter : 1
                    )
                    relaxedAdded++
                }
                if (a.score >= 0.5d) {
                    integrated.annotations.add(a)
                }
            }
        }
        // Rebuild EC entries in integrated.provenance directly from the enforced genome proteins
        // so EC removals, downranks, and flags are 100% reflected.
        Set<String> ecKeysToRemove = integrated.provenance.keySet().findAll { k ->
            String[] fkey = IntegrationState.splitFunctionKey(k)
            fkey != null && fkey[1] == 'EC'
        } as Set
        Map<String, gspa.integration.ClaimProvenance> oldEcProv = [:]
        for (String k : ecKeysToRemove) {
            oldEcProv[k] = integrated.provenance.remove(k)
        }
        integrated.annotations.annotations.removeAll { it.type == AnnotationType.EC }

        for (Protein p : contig.proteins) {
            for (Annotation a : p.annotations.ecAnnotations()) {
                String fkey = IntegrationState.functionKey(p.id, 'EC', a.value)
                def existing = oldEcProv[fkey]
                if (existing != null) {
                    existing.finalProbability = a.score
                    existing.finalLogOdds = CalibrationTable.probToLogOdds(existing.finalProbability)
                    if (a.metadata?.get('violating_go_terms')) {
                        existing.priorContributions['stage2_ec_taxon'] = -1.0d
                    }
                    integrated.provenance[fkey] = existing
                }
                if (a.score >= 0.5d) {
                    integrated.annotations.add(a)
                }
            }
        }
        if (effMode?.toLowerCase() == 'relax') {
            Map<String, List<gspa.integration.ClaimProvenance>> byProtAspect = [:].withDefault { [] }
            for (def prov : integrated.provenance.values()) {
                String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey == null || fkey[1] != 'GO') continue
                String aspect = (prov.supportingClaims && prov.supportingClaims.first().goAspect)
                    ? prov.supportingClaims.first().goAspect
                    : (state.goOntology != null ? (state.goOntology.getAspect(fkey[2]) ?: 'U') : 'U')
                byProtAspect["${fkey[0]}|${aspect}"] << prov
            }
            for (List<gspa.integration.ClaimProvenance> list : byProtAspect.values()) {
                if (list.isEmpty()) continue
                double maxProb = 0.0d
                int maxNsup = 0
                gspa.integration.ClaimProvenance topProv = null
                for (def pr : list) {
                    int nsup = pr.supportingClaims != null ? pr.supportingClaims.size() : 1
                    if (pr.finalProbability > maxProb) {
                        maxProb = pr.finalProbability
                        maxNsup = nsup
                        topProv = pr
                    }
                }
                if (topProv != null && maxProb >= 0.25d && maxProb < 0.48d && maxNsup >= 2) {
                    double newLo = topProv.finalLogOdds + 0.32d
                    double newProb = Math.min(0.54d, CalibrationTable.logOddsToProb(newLo))
                    topProv.finalProbability = newProb
                    topProv.finalLogOdds = CalibrationTable.probToLogOdds(newProb)
                }
            }
        }
        // Obligate heteromeric complex singleton pruning: if an obligate heteromer ('n') has only a single
        // active protein (>= 0.50) and no partner candidate (>= 0.08) across the proteome, attenuate it to 0.45.
        if (state.complexClassification != null && !state.complexClassification.isEmpty()) {
            Map<String, List<String>> activeHeteromerByTerm = [:].withDefault { [] }
            Map<String, List<String>> candHeteromerByTerm = [:].withDefault { [] }
            for (def prov : integrated.provenance.values()) {
                String[] fkey = IntegrationState.splitFunctionKey(prov.functionKey)
                if (fkey == null || fkey[1] != 'GO') continue
                String goId = fkey[2]
                if (state.complexClassification[goId] != 'n') continue
                if (prov.finalProbability >= 0.5d) {
                    activeHeteromerByTerm[goId] << prov.functionKey
                } else if (prov.finalProbability >= 0.08d) {
                    candHeteromerByTerm[goId] << prov.functionKey
                }
            }
            for (def entry : activeHeteromerByTerm.entrySet()) {
                String goId = entry.key
                List<String> activeKeys = entry.value
                List<String> candKeys = candHeteromerByTerm[goId] ?: []
                if (activeKeys.size() == 1 && candKeys.isEmpty()) {
                    def singleProv = integrated.provenance[activeKeys[0]]
                    if (singleProv != null && singleProv.finalProbability > 0.45d) {
                        double beforeProb = singleProv.finalProbability
                        singleProv.finalProbability = 0.45d
                        singleProv.finalLogOdds = CalibrationTable.probToLogOdds(0.45d)
                        singleProv.priorContributions['stage2_singleton_prune'] = -0.5d
                        report.record(new EnforcementAction(
                            dimension: 'coherence',
                            action: 'PRUNED_SINGLETON',
                            contigId: 'contig_1',
                            proteinId: singleProv.proteinId,
                            term: goId,
                            reason: 'obligate heteromer singleton without partner candidate (>=0.08)',
                            basis: 'complexClassification=n',
                            scoreBefore: beforeProb,
                            scoreAfter: 0.45d
                        ))
                        integrated.annotations.annotations.removeAll {
                            it.type == AnnotationType.GO && it.value == goId && singleProv.proteinId == (it.metaClass.respondsTo(it, 'getProteinId') ? it.proteinId : singleProv.proteinId)
                        }
                    }
                }
            }
        }
        integrated.synchronize(state)
        if (verbose) println "  Stage 2 synced to integrated output (relaxed_new=${relaxedAdded})"

        if (effActionsOut != null && report.count() > 0) {
            if (state.goOntology != null) {
                report.actions.each { it.termLabel = state.goOntology.getLabel(it.term) }
            }
            effActionsOut.parentFile?.mkdirs()
            report.writeTsv(effActionsOut)
            if (verbose) println "  Stage 2 enforcement actions: ${effActionsOut}"
        }
    }
}
