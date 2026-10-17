package gspa.metrics

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.model.Protein
import gspa.ontology.AnatomicalConstraints
import gspa.ontology.GoOntology
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConsistencyMaxSat
import groovy.util.logging.Slf4j

/**
 * Post-annotation taxon and anatomical constraint enforcement on predicted GO annotations.
 *
 * Modes:
 * <ul>
 *   <li>{@code remove} / {@code downrank} / {@code flag} — act on each
 *       annotation found taxon- or anatomically-inconsistent. With an organism asserted on the
 *       {@link SatConsistencyChecker}, a term is inconsistent iff it (or an
 *       ancestor) cannot occur in the organism's lineage (precise per-term).
 *       Without one, per-protein co-annotation conflicts are used.</li>
 *   <li>{@code relax} — Semantic Ancestor Relaxation: replaces each violating GO term
 *       with its closest non-root ancestor in the GO hierarchy ({@code is_a}/{@code part_of})
 *       that satisfies all taxonomic and anatomical constraints for the organism,
 *       or removes the annotation if no valid non-root ancestor exists.</li>
 *   <li>{@code minimal-flip} — joint, minimum-cost demotion via
 *       {@link TaxonConsistencyMaxSat} (Asaad Stage-1): resolves co-annotation
 *       conflicts optimally (drops the lower-weight of two disjoint
 *       requirements) so no residual inconsistency remains.</li>
 * </ul>
 *
 * Every action is recorded in the optional {@link EnforcementReport} for
 * provenance (and removed annotations, which leave no trace in the set, are
 * captured there).
 */
@Slf4j
class ConsistencyEnforcer {

    GoOntology goOntology
    SatConsistencyChecker satChecker
    AnatomicalConstraints anatomicalConstraints
    String mode = 'remove'
    double downrankFactor = 0.5
    double relaxFactor = 1.0d

    /** Optional provenance sink; when set, every action is logged and trail lines appended. */
    EnforcementReport report

    static class Result {
        int proteinsAffected = 0
        int annotationsAffected = 0
        int annotationsRelaxed = 0
        int violations = 0

        @Override
        String toString() {
            "ConsistencyEnforcer.Result(violations=${violations}, " +
                "proteins=${proteinsAffected}, annotations=${annotationsAffected}, relaxed=${annotationsRelaxed})"
        }
    }

    Result enforce(Genome genome) {
        def result = new Result()
        if (mode?.toLowerCase() == 'minimal-flip') {
            enforceMinimalFlip(genome, result)
        } else if (satChecker?.organismTaxon) {
            enforceOrganism(genome, result)
        } else {
            genome.proteins.each { Protein p -> enforceProtein(p, result) }
        }
        log.info("Consistency enforcement (${mode}): acted on ${result.annotationsAffected} annotation(s) " +
            "(relaxed=${result.annotationsRelaxed}) across ${result.proteinsAffected} protein(s); " +
            "${result.violations} violation(s) found")
        result
    }

    // --- minimal-flip (joint MaxSAT) -----------------------------------------

    private void enforceMinimalFlip(Genome genome, Result result) {
        // Pass A: with an organism asserted, drop the bulk of impossible terms
        // per-term (fast) so the joint solver only faces the residual.
        if (satChecker?.organismTaxon) {
            enforceOrganism(genome, result)
            log.info("  minimal-flip Pass A done: ${result.annotationsAffected} removed; resolving residual...")
        }
        // Pass B: iteratively resolve residual co-annotation conflicts. Each
        // round, check() returns the minimal UNSAT-core terms; the weighted
        // MaxSAT runs over just those (tiny) and drops the lower-weight one.
        // This keeps the solver bounded — never the whole survivor set.
        int maxRounds = 100
        for (int round = 0; round < maxRounds; round++) {
            Set<String> survivors = genome.allGoTerms()
            def cr = satChecker.check(survivors)
            if (cr.consistent) { log.info("  minimal-flip: consistent after ${round} residual round(s)"); break }
            Set<String> conflict = (cr.violations.collectMany { it.involvedGoTerms ?: [] } as Set)
            conflict.retainAll(survivors)
            if (conflict.isEmpty()) break

            Map<String, Double> w = [:].withDefault { 0.0d }
            genome.proteins.each { Protein p ->
                p.annotations.goAnnotations().each {
                    if (conflict.contains(it.value)) w[it.value] += Math.max(0.0d, it.score)
                }
            }
            Set<String> remove = new TaxonConsistencyMaxSat(checker: satChecker, goOntology: goOntology)
                .termsToRemove(w)
            if (remove.isEmpty() && w) {
                // MaxSAT declined (e.g. all-equal weights); drop the single
                // lowest-weight conflict term to guarantee progress.
                remove = [w.min { it.value }.key] as Set
            }
            if (remove.isEmpty()) break
            genome.proteins.each { Protein p ->
                def victims = p.annotations.goAnnotations().findAll { remove.contains(it.value) }
                if (victims) {
                    result.violations++
                    applyAction(p, victims, result, 'remove', 'co-annotation taxon conflict', null)
                }
            }
        }
    }

    // --- organism-asserted per-term ------------------------------------------

    private void enforceOrganism(Genome genome, Result result) {
        Set<String> forbidden = computeForbiddenTerms()
        log.info("Organism ${satChecker.organismTaxon}: ${forbidden.size()} GO terms forbidden by taxon/anatomical constraints")
        String defaultReason = "taxon-inconsistent for ${satChecker.organismTaxon}"
        Map<String, Boolean> termForbiddenMemo = new HashMap<>()
        genome.proteins.each { Protein p ->
            List<Annotation> victims = p.annotations.goAnnotations().findAll { ann ->
                Boolean memo = termForbiddenMemo.get(ann.value)
                if (memo != null) return memo
                boolean isBad = forbidden.contains(ann.value) ||
                    (goOntology && goOntology.getAncestors(ann.value).any { forbidden.contains(it) })
                termForbiddenMemo.put(ann.value, isBad)
                isBad
            }
            if (victims) {
                result.violations++
                applyAction(p, victims, result, effectiveAction(), defaultReason, forbidden)
            }
        }
        if (effectiveAction() in ['remove', 'relax']) {
            resolveResidualCoAnnotationConflicts(genome, result, forbidden)
        }
    }

    private void resolveResidualCoAnnotationConflicts(Genome genome, Result result, Set<String> forbidden) {
        int maxRounds = 50
        for (int round = 0; round < maxRounds; round++) {
            Set<String> activeTerms = new LinkedHashSet<>()
            genome.proteins.each { Protein p ->
                p.annotations.goAnnotations().each { Annotation a ->
                    if (a.score >= 0.5d) activeTerms << a.value
                }
            }
            if (activeTerms.isEmpty()) break
            Set<String> checkSet = goOntology != null ? goOntology.propagateAnnotations(activeTerms) : activeTerms
            def cr = satChecker.check(checkSet)
            if (cr.consistent) break
            Set<String> conflict = (cr.violations.collectMany { it.involvedGoTerms ?: [] } as Set)
            if (conflict.isEmpty()) break

            Map<String, Double> w = [:].withDefault { 0.0d }
            genome.proteins.each { Protein p ->
                p.annotations.goAnnotations().each { Annotation a ->
                    if (a.score >= 0.5d && annotationViolates(a, conflict)) {
                        w[a.value] += Math.max(0.0d, a.score)
                    }
                }
            }
            if (w.isEmpty()) break
            Set<String> remove = new TaxonConsistencyMaxSat(checker: satChecker, goOntology: goOntology)
                .termsToRemove(w)
            if (remove.isEmpty()) {
                remove = [w.min { it.value }.key] as Set
            }
            Set<String> nextForbidden = new HashSet<>(forbidden ?: Collections.emptySet())
            nextForbidden.addAll(remove)
            nextForbidden.addAll(conflict)
            boolean acted = false
            genome.proteins.each { Protein p ->
                List<Annotation> victims = p.annotations.goAnnotations().findAll { Annotation a ->
                    a.score >= 0.5d && (remove.contains(a.value) || annotationViolates(a, remove))
                }
                if (victims) {
                    result.violations++
                    applyAction(p, victims, result, effectiveAction(), 'co-annotation taxon conflict', nextForbidden)
                    acted = true
                }
            }
            if (!acted) break
        }
    }

    private Set<String> computeForbiddenTerms() {
        Set<String> forbidden = new HashSet<>(satChecker.getForbiddenConstrainedTerms())
        if (anatomicalConstraints != null && satChecker.organismTaxon) {
            forbidden.addAll(anatomicalConstraints.computeForbiddenGoTerms(satChecker.organismTaxon, satChecker))
        }
        forbidden
    }

    // --- per-protein co-annotation -------------------------------------------

    private void enforceProtein(Protein p, Result result) {
        List<Annotation> go = p.annotations.goAnnotations()
        if (go.size() < 2) return
        Set<String> direct = go*.value as Set
        Set<String> toCheck = goOntology ? goOntology.propagateAnnotations(direct) : direct
        def cr = satChecker.check(toCheck)
        if (cr.consistent) return
        result.violations += cr.violations.size()
        Set<String> conflicting = cr.violations.collectMany { it.involvedGoTerms ?: [] } as Set
        if (!conflicting) return
        List<Annotation> victims = go.findAll { ann -> annotationViolates(ann, conflicting) }
        if (victims) applyAction(p, victims, result, effectiveAction(), 'co-annotation taxon conflict', conflicting)
    }

    private boolean annotationViolates(Annotation ann, Set<String> conflicting) {
        if (conflicting.contains(ann.value)) return true
        if (goOntology == null) return false
        goOntology.propagateAnnotations([ann.value] as Set).any { conflicting.contains(it) }
    }

    // --- actions + provenance ------------------------------------------------

    private String effectiveAction() {
        switch (mode?.toLowerCase()) {
            case 'downrank': return 'downrank'
            case 'flag': return 'flag'
            case 'relax': return 'relax'
            default: return 'remove'
        }
    }

    private void applyAction(Protein p, List<Annotation> victims, Result result, String action,
                             String defaultReason, Set<String> forbidden) {
        if (action == 'relax') {
            applyRelaxAction(p, victims, result, defaultReason, forbidden)
            return
        }
        if (action == 'remove') {
            // Batch removal — O(n) instead of O(n*victims) of per-element remove.
            p.annotations.annotations.removeAll(new HashSet<>(victims))
        }
        victims.each { Annotation a ->
            double before = a.score
            String reason = resolveReason(a.value, defaultReason)
            if (action == 'downrank') {
                a.score = a.score * downrankFactor
                a.evidence = tag(a.evidence)
            } else if (action == 'flag') {
                a.evidence = tag(a.evidence)
            }
            record(p, a, action, reason, before, action == 'remove' ? 0.0d : a.score)
        }
        result.proteinsAffected++
        result.annotationsAffected += victims.size()
    }

    private void applyRelaxAction(Protein p, List<Annotation> victims, Result result,
                                  String defaultReason, Set<String> forbidden) {
        p.annotations.annotations.removeAll(new HashSet<>(victims))
        Map<String, Annotation> existingGo = p.annotations.goAnnotations().collectEntries { [it.value, it] }
        AnatomicalConstraints relaxer = anatomicalConstraints ?: new AnatomicalConstraints()
        String orgTax = satChecker?.organismTaxon

        victims.each { Annotation a ->
            double before = a.score
            double relaxedScore = before * relaxFactor
            String baseReason = resolveReason(a.value, defaultReason)
            List<String> relaxedAncestors = (goOntology != null && orgTax) ?
                relaxer.findValidRelaxedAncestors(a.value, orgTax, satChecker, goOntology, forbidden) : Collections.emptyList()
            if (!relaxedAncestors.isEmpty()) {
                for (String relaxedGo : relaxedAncestors) {
                    Annotation existing = existingGo[relaxedGo]
                    if (existing != null) {
                        if (relaxedScore > existing.score) {
                            existing.score = relaxedScore
                        }
                        existing.evidence = tagRelaxed(existing.evidence, a.value)
                        record(p, existing, 'relax', "${baseReason} -> relaxed to ${relaxedGo}", before, existing.score)
                    } else {
                        Annotation replacement = new Annotation(
                            type: AnnotationType.GO,
                            value: relaxedGo,
                            score: relaxedScore,
                            source: a.source,
                            evidence: tagRelaxed(a.evidence, a.value),
                            goAspect: a.goAspect ?: (goOntology?.getAspect(relaxedGo)),
                            goLabel: goOntology?.getLabel(relaxedGo)
                        )
                        p.annotations.add(replacement)
                        existingGo[relaxedGo] = replacement
                        record(p, replacement, 'relax', "${baseReason} -> relaxed to ${relaxedGo}", before, relaxedScore)
                    }
                }
                result.annotationsRelaxed++
            } else {
                record(p, a, 'remove', "${baseReason} (no valid non-root ancestor)", before, 0.0d)
            }
        }
        result.proteinsAffected++
        result.annotationsAffected += victims.size()
    }

    private String resolveReason(String goId, String defaultReason) {
        if (anatomicalConstraints != null && satChecker?.organismTaxon) {
            def details = anatomicalConstraints.checkTerm(goId, satChecker.organismTaxon, satChecker, goOntology)
            if (!details.isEmpty()) {
                return "anatomical-inconsistent for ${satChecker.organismTaxon}: ${details[0].summary()}"
            }
        }
        defaultReason
    }

    private void record(Protein p, Annotation a, String action, String reason, double before, double after) {
        if (report == null) return
        def ea = new EnforcementAction(
            dimension: 'consistency', action: action,
            contigId: p.contig?.id, proteinId: p.id, term: a.value,
            reason: reason, basis: "${a.source ?: 'predictor'}@${String.format(Locale.ROOT, '%.3f', before)}",
            scoreBefore: before, scoreAfter: after)
        // For downrank/flag/relax the annotation survives -> attach trail; for remove only log.
        report.record(ea, action == 'remove' ? null : a)
    }

    private static String tag(String evidence) {
        evidence ? "${evidence};taxon-inconsistent" : 'taxon-inconsistent'
    }

    private static String tagRelaxed(String evidence, String fromGoId) {
        String suffix = "relaxed-from:${fromGoId}"
        evidence ? "${evidence};${suffix}" : suffix
    }
}
