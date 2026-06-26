package gspa.metrics

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.model.Protein
import gspa.ontology.GoOntology
import gspa.ontology.SatConsistencyChecker
import groovy.util.logging.Slf4j

/**
 * Optional coherence enforcement (Asaad Stage-2 + process coherence):
 *
 * <ul>
 *   <li><b>Complex coherence.</b> An obligate heteromeric-complex GO term
 *       (classification {@code 'n'}) annotated to exactly ONE protein is
 *       incoherent (a heteromer needs &ge;2 distinct subunits). It is either
 *       <i>promoted</i> onto a plausible partner protein (best near-ancestor
 *       evidence) or, if {@link #promotePartner} is off or no partner is found,
 *       <i>demoted</i> from the lone protein.</li>
 *   <li><b>Process coherence.</b> For an unsatisfied {@code has_part} pair
 *       (C present, partner F missing), F is promoted onto the best candidate
 *       protein.</li>
 * </ul>
 *
 * Promotions reuse {@link PromotionCandidates}; terms the organism cannot carry
 * are never promoted. Every action is recorded for provenance.
 */
@Slf4j
class CoherenceEnforcer {

    GoOntology goOntology
    SatConsistencyChecker satChecker
    EnforcementReport report

    /** GO term -> 'h' (homomeric) | 'n' (heteromeric) | 'a' (ambiguous). */
    Map<String, String> complexClassification = [:]

    boolean promotePartner = true
    boolean enforceProcess = true
    String promoteEvidence = 'ISC'
    double promoteMinScore = 0.05
    int maxHops = 1
    double hopDecay = 0.5

    static class Result {
        int complexDemoted = 0
        int complexPromoted = 0
        int processPromoted = 0
        @Override
        String toString() {
            "CoherenceEnforcer.Result(complexDemoted=${complexDemoted}, " +
                "complexPromoted=${complexPromoted}, processPromoted=${processPromoted})"
        }
    }

    /**
     * @param unsatisfiedProcessPairs (C, F) has_part pairs with C present, F missing
     */
    Result enforce(Genome genome, List<Map.Entry<String, String>> unsatisfiedProcessPairs = []) {
        def result = new Result()
        enforceComplex(genome, result)
        if (enforceProcess) enforceProcessCoherence(genome, unsatisfiedProcessPairs, result)
        log.info("Coherence enforcement: ${result}")
        result
    }

    private void enforceComplex(Genome genome, Result result) {
        if (complexClassification.isEmpty()) return
        // Proteins annotated with each complex term.
        Map<String, List<Protein>> byComplex = [:].withDefault { [] }
        genome.proteins.each { Protein p ->
            p.annotations.goAnnotations()*.value.unique().each { String t ->
                if (complexClassification[t] == 'n') byComplex[t] << p
            }
        }
        byComplex.each { String complexTerm, List<Protein> proteins ->
            if (proteins.size() != 1) return   // coherent or absent
            Protein lone = proteins[0]
            def partner = promotePartner ?
                PromotionCandidates.best(genome, complexTerm, goOntology, maxHops, hopDecay,
                    promoteMinScore, [lone.id] as Set) : null
            if (partner != null && !organismForbids(complexTerm)) {
                promote(partner.protein as Protein, complexTerm, partner.score as double,
                    "complex partner; ${partner.basis}", 'coherence', result, { result.complexPromoted++ })
            } else {
                // Demote the singleton from the lone protein.
                def victim = lone.annotations.goAnnotations().find { it.value == complexTerm }
                if (victim) {
                    double before = victim.score
                    lone.annotations.annotations.remove(victim)
                    result.complexDemoted++
                    report?.record(new EnforcementAction(
                        dimension: 'coherence', action: 'remove',
                        contigId: lone.contig?.id, proteinId: lone.id, term: complexTerm,
                        reason: 'obligate heteromeric complex with a single subunit',
                        basis: "${victim.source ?: 'predictor'}@${String.format(Locale.ROOT, '%.3f', before)}",
                        scoreBefore: before, scoreAfter: 0.0d))
                }
            }
        }
    }

    private void enforceProcessCoherence(Genome genome, List<Map.Entry<String, String>> pairs, Result result) {
        if (!pairs) return
        Set<String> done = new HashSet<>()
        pairs.each { Map.Entry<String, String> pair ->
            String f = pair.value   // missing partner
            if (!f || !done.add(f) || organismForbids(f)) return
            def cand = PromotionCandidates.best(genome, f, goOntology, maxHops, hopDecay, promoteMinScore)
            if (cand == null) return
            promote(cand.protein as Protein, f, cand.score as double,
                "has_part partner of ${pair.key}; ${cand.basis}", 'coherence', result, { result.processPromoted++ })
        }
    }

    private void promote(Protein p, String term, double score, String basis, String dim,
                         Result result, Closure tally) {
        def ann = new Annotation(type: AnnotationType.GO, value: term, score: score,
            source: 'gspa-enforce:coherence', evidence: promoteEvidence)
        p.annotations.add(ann)
        tally()
        report?.record(new EnforcementAction(
            dimension: dim, action: 'promote',
            contigId: p.contig?.id, proteinId: p.id, term: term,
            reason: 'coherence repair', basis: basis,
            scoreBefore: 0.0d, scoreAfter: score), ann)
    }

    private boolean organismForbids(String f) {
        if (satChecker?.organismTaxon == null) return false
        !satChecker.check([f] as Set).consistent
    }
}
