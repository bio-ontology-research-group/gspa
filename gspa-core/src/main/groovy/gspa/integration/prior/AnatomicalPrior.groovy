package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.Prior
import gspa.ontology.AnatomicalConstraints
import gspa.ontology.AnatomicalConstraints.AnatomicalLink
import gspa.ontology.SatConsistencyChecker
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Genome-scale anatomical & developmental constraint prior.
 *
 * <p>Operates over GO-Uberon-PO-CL cross-product axioms ({@code results_in_development_of},
 * {@code results_in_morphogenesis_of}, {@code results_in_formation_of}, {@code occurs_in},
 * {@code part_of}) and the taxonomic range of each anatomical structure:</p>
 * <ul>
 *   <li><b>Anatomical penalty ({@code -alphaAnatPenalty}):</b> Downweights candidate GO
 *       claims whose target anatomical structure (or any GO ancestor's anatomical
 *       structure) cannot exist in the organism's clade ({@code never_in_taxon} or
 *       disjoint {@code only_in_taxon}).</li>
 *   <li><b>Anatomical module coherence boost ({@code +alphaAnatBoost}):</b> When an
 *       anatomical structure {@code S} is valid in the organism's clade AND strongly
 *       evidenced across the genome (&ge; {@link #minSupportProteins} distinct proteins
 *       with posterior &gt; {@link #supportThreshold} for processes acting on {@code S}),
 *       sub-threshold candidate claims linked to {@code S} receive a modest positive
 *       log-odds boost.</li>
 * </ul>
 */
class AnatomicalPrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(AnatomicalPrior)

    /** Soft penalty in log-odds for anatomically impossible GO terms. */
    double alphaAnatPenalty = 0.68d

    /** Soft boost in log-odds for candidate claims on an active, valid anatomical structure. */
    double alphaAnatBoost = 0.65d

    /** Minimum number of strongly-annotated proteins required to consider an anatomical structure active. */
    int minSupportProteins = 2

    /** Posterior probability threshold for counting a protein as strong evidence for a structure. */
    double supportThreshold = 0.50d

    /** If true, effectively remove anatomically forbidden claims (-1000 log-odds). */
    boolean hardFilter = false

    private Set<String> forbiddenGoTerms = Collections.emptySet()
    private Set<String> boostedGoTerms = Collections.emptySet()

    @Override
    String name() { 'anatomical' }

    @Override
    Set<String> inputs() { ['anatomicalConstraints', 'satConsistencyChecker', 'posteriors'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        forbiddenGoTerms = Collections.emptySet()
        boostedGoTerms = Collections.emptySet()

        AnatomicalConstraints ac = state.anatomicalConstraints
        SatConsistencyChecker checker = state.satConsistencyChecker
        if (ac == null || checker == null || !checker.organismTaxon) {
            log.debug('AnatomicalPrior disabled: missing anatomicalConstraints or organismTaxon')
            return
        }

        String orgTax = checker.organismTaxon
        Set<String> directForbidden = ac.computeForbiddenGoTerms(orgTax, checker)
        Set<String> allForbidden = new LinkedHashSet<>(directForbidden)

        // Check candidate GO terms in posteriorLogOdds for ancestor-inherited anatomical violations
        if (state.goOntology != null && !directForbidden.isEmpty()) {
            Set<String> candidateGo = new LinkedHashSet<>()
            for (String fkey : state.posteriorLogOdds.keySet()) {
                String[] parts = IntegrationState.splitFunctionKey(fkey)
                if (parts != null && parts[1] == 'GO') candidateGo.add(parts[2])
            }
            for (String goId : candidateGo) {
                if (allForbidden.contains(goId)) continue
                Set<String> ancs = state.goOntology.getAncestors(goId)
                if (ancs.any { directForbidden.contains(it) }) {
                    allForbidden.add(goId)
                }
            }
        }
        forbiddenGoTerms = allForbidden

        // Identify active, valid anatomical structures supported by >= minSupportProteins proteins
        Set<String> genericRoots = [
            'CL:0000000', 'UBERON:0001062', 'UBERON:0000061', 'UBERON:0000465',
            'UBERON:0000468', 'UBERON:0010000', 'UBERON:0000062', 'UBERON:0000479',
            'PO:0025131', 'PO:0009011'
        ] as Set
        Map<String, Set<String>> structureProteins = [:].withDefault { new LinkedHashSet<String>() }
        for (Map.Entry<String, Double> entry : state.posteriorLogOdds.entrySet()) {
            double likeProb = state.getLikelihoodProbability(entry.key)
            if (likeProb < supportThreshold) continue
            String[] parts = IntegrationState.splitFunctionKey(entry.key)
            if (parts == null || parts[1] != 'GO') continue
            String protId = parts[0]
            String goId = parts[2]
            if (forbiddenGoTerms.contains(goId)) continue
            List<AnatomicalLink> links = ac.goToAnatomy[goId]
            if (links) {
                for (AnatomicalLink link : links) {
                    if (!genericRoots.contains(link.anatomyId)) {
                        structureProteins[link.anatomyId].add(protId)
                        for (String pAnat : (ac.anatomyParents[link.anatomyId] ?: Collections.emptySet())) {
                            if (!genericRoots.contains(pAnat)) {
                                structureProteins[pAnat].add(protId)
                            }
                        }
                    }
                }
            }
        }

        Set<String> activeStructures = structureProteins.findAll { it.value.size() >= minSupportProteins }.keySet()
        if (!activeStructures.isEmpty()) {
            Set<String> boosted = new LinkedHashSet<>()
            for (Map.Entry<String, List<AnatomicalLink>> e : ac.goToAnatomy.entrySet()) {
                String goId = e.key
                if (forbiddenGoTerms.contains(goId)) continue
                if (e.value.any { link ->
                    activeStructures.contains(link.anatomyId) ||
                        (ac.anatomyParents[link.anatomyId] ?: Collections.emptySet()).any { activeStructures.contains(it) }
                }) {
                    boosted.add(goId)
                }
            }
            boostedGoTerms = boosted
        }

        log.debug("AnatomicalPrior (${orgTax}): ${forbiddenGoTerms.size()} forbidden GO terms, " +
            "${activeStructures.size()} active structures (${boostedGoTerms.size()} boosted GO terms)")
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (forbiddenGoTerms.isEmpty() && boostedGoTerms.isEmpty()) return 0.0d
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null || parts[1] != 'GO') return 0.0d
        String goId = parts[2]
        if (forbiddenGoTerms.contains(goId)) {
            return hardFilter ? -1000.0d : -alphaAnatPenalty
        }
        if (boostedGoTerms.contains(goId)) {
            return alphaAnatBoost
        }
        return 0.0d
    }
}
