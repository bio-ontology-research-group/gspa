package gspa.metrics

import gspa.ontology.GoOntology
import gspa.ontology.SatConsistencyChecker
import groovy.util.logging.Slf4j

/**
 * Infer an organism's taxon from its predicted GO annotations using the GO taxon
 * constraints — the Asaad et al. (genome-scale-pfp-adjust) idea applied in
 * reverse. Each predicted term carries {@code only_in_taxon} / {@code
 * never_in_taxon} restrictions; the organism's true taxon is the one under which
 * the prediction set is (most) consistent, and as <em>specific</em> as the
 * evidence allows.
 *
 * <p><b>Beyond the four domains.</b> The candidate taxa are not a fixed
 * kingdom list — they are the {@code only_in_taxon} targets that actually appear
 * among the confident predictions (a confident term whose function
 * {@code only_in_taxon T} is direct evidence the organism is within {@code T}),
 * together with the four domains as anchors. Any taxon named by the bundled GO
 * taxon constraints (Vertebrata, Fungi, Mammalia, <i>Saccharomyces cerevisiae</i>,
 * <i>Homo sapiens</i>, …) can therefore be inferred, not just Bacteria/Archaea/
 * Eukaryota/Viruses. Among the candidates that stay consistent we pick the most
 * specific (deepest in the hierarchy) — so a vertebrate proteome infers
 * Vertebrata/Mammalia rather than merely Eukaryota.</p>
 *
 * <p><b>Confidence matters.</b> A learned predictor over-predicts a long tail of
 * low-score, often cross-domain, terms; a raw count would then be biased toward
 * whichever taxon is most permissive of that noise. So inference runs over the
 * <i>high-confidence</i> predictions only ({@link #minScore}).</p>
 *
 * <p>Counting is done directly over the taxon hierarchy (lineage + explicit
 * disjointness) rather than one SAT solve per (candidate, term): for the
 * forced-unit Horn structure the SAT encoding uses, a term is taxon-inconsistent
 * under an assumed organism {@code T} iff (a) some {@code never_in} taxon lies on
 * {@code T}'s (or a required {@code only_in} target's) lineage, or (b) two taxa
 * forced true are explicitly disjoint. {@code TaxonInferenceSatEquivalenceSpec}
 * locks this against the real {@link SatConsistencyChecker}.</p>
 */
@Slf4j
class TaxonInference {

    SatConsistencyChecker checker
    GoOntology goOntology
    /** Only predictions at or above this score inform the inference. */
    double minScore = 0.5

    /** Domain anchors — always in the candidate pool so a gross-level call is always possible. */
    static final List<String> DOMAIN_TAXA = [
        'NCBITaxon_2',     // Bacteria
        'NCBITaxon_2157',  // Archaea
        'NCBITaxon_2759',  // Eukaryota
        'NCBITaxon_10239', // Viruses
    ]

    /**
     * Universal ancestors that are never a useful organism call — the root and
     * "cellular organisms" are consistent with everything (zero violations) and
     * would otherwise win by default whenever a specific lineage carries any
     * cross-taxon prediction noise. Excluded as candidates (their lineage still
     * informs support / consistency of real candidates).
     */
    static final Set<String> NON_ORGANISM_TAXA = ['NCBITaxon_1', 'NCBITaxon_131567'] as Set

    /**
     * A specific (sub-domain) taxon must draw at least this many supporting terms
     * to be eligible to win, so a single cross-taxon noise term can't crown an
     * unrelated genus/class. The four domains are always eligible.
     */
    static final int MIN_SPECIFIC_SUPPORT = 3

    /** NCBITaxon_<id> -> scientific name, loaded once from the bundled resource. */
    private static Map<String, String> LABELS = null

    static synchronized Map<String, String> labels() {
        if (LABELS != null) return LABELS
        Map<String, String> m = [:]
        def is = TaxonInference.getResourceAsStream('/taxon-constraints/ncbi-taxon-labels.tsv')
        if (is != null) {
            is.withReader('UTF-8') { r ->
                r.eachLine { line, n ->
                    if (n == 1 || !line?.trim() || line.startsWith('#')) return
                    def f = line.split('\t')
                    if (f.length >= 2) m[f[0].trim()] = f[1].trim()
                }
            }
        } else {
            log.warn("ncbi-taxon-labels.tsv not found on classpath — taxa shown by id")
        }
        LABELS = m
        m
    }

    /** Human-readable name for a taxon, or the raw id if unknown. */
    static String labelFor(String taxon) { labels()[taxon] ?: taxon }

    static class Candidate {
        String taxon
        String label
        int forbidden        // constraint-bearing present terms inconsistent under this taxon
        int support          // present terms whose only_in target lies on this taxon's lineage
        int depth            // size of the taxon's lineage (higher = more specific)
    }

    static class Result {
        String taxon         // inferred NCBITaxon_<id>, or null if undecidable
        String label         // human label, or 'Unknown'
        int constrainedPresent
        boolean confident
        List<Candidate> candidates = []   // ranked, best first
        List<String> lineage = []         // inferred taxon's ancestry (specific -> general), labelled
    }

    /**
     * Infer the taxon from per-term prediction scores (term &rarr; max score
     * across the proteome). Only terms scoring &ge; {@link #minScore} (and their
     * propagated ancestors) inform the call.
     */
    Result infer(Map<String, Double> termScore) {
        Set<String> highConf = termScore.findAll { it.value >= minScore }.keySet()
        Set<String> present = goOntology ? goOntology.propagateAnnotations(highConf) : highConf
        def tc = checker.taxonConstraints
        Set<String> constrainedKeys = (tc.onlyInTaxon.keySet() + tc.neverInTaxon.keySet()) as Set
        def relevant = present.findAll { constrainedKeys.contains(it) }.toList()

        def result = new Result(constrainedPresent: relevant.size())
        if (relevant.isEmpty()) {
            result.label = 'Unknown'
            log.info("Taxon inference: no constraint-bearing terms present — undecidable")
            return result
        }

        // Candidate pool: the only_in targets named by the confident predictions
        // (direct positive evidence "organism is within T") plus the domain
        // anchors. Union grouping nodes are not organisms, so they are excluded
        // as candidates (their lineage still informs support / consistency).
        Set<String> pool = new LinkedHashSet<>()
        relevant.each { t -> tc.onlyInTaxon[t]?.each { tax ->
            if (isRealTaxon(tax) && !NON_ORGANISM_TAXA.contains(tax)) pool.add(tax)
        } }
        pool.addAll(DOMAIN_TAXA)

        pool.each { tax ->
            def anc = checker.ancestorsWithSelf(tax)
            int forbidden = 0, support = 0
            relevant.each { term ->
                if (violates(tax, term)) {
                    forbidden++
                } else if (tc.onlyInTaxon.containsKey(term) &&
                           tc.onlyInTaxon[term].any { anc.contains(it) }) {
                    support++   // a confident only_in target on this taxon's lineage
                }
            }
            result.candidates << new Candidate(taxon: tax, label: labelFor(tax),
                forbidden: forbidden, support: support, depth: anc.size())
        }

        // Best = fewest violations, then most specific (deepest), then most support.
        result.candidates.sort { a, b ->
            (a.forbidden <=> b.forbidden) ?: (b.depth <=> a.depth) ?: (b.support <=> a.support)
        }
        // Eligible to win: the four domains always; any more specific taxon only if
        // it has real positive support (so one noise term can't crown a stray genus).
        def eligible = result.candidates.findAll {
            DOMAIN_TAXA.contains(it.taxon) || it.support >= MIN_SPECIFIC_SUPPORT
        }
        def best = (eligible ?: result.candidates)[0]
        result.taxon = best.taxon
        result.label = best.label
        result.lineage = lineageOf(best.taxon)

        // Confident when the winner is violation-free, has positive evidence, and
        // no equally-consistent candidate in a different lineage is at least as
        // specific (which would mean the data points two ways).
        def rivals = eligible.findAll {
            it.taxon != best.taxon && it.forbidden == best.forbidden &&
            it.depth >= best.depth && !checker.ancestorsWithSelf(best.taxon).contains(it.taxon) &&
            !checker.ancestorsWithSelf(it.taxon).contains(best.taxon)
        }
        result.confident = best.forbidden == 0 && best.support > 0 && rivals.isEmpty()

        log.info("Taxon inference: ${best.label} (${best.taxon}) from ${relevant.size()} " +
            "constraint-bearing terms; top=" +
            result.candidates.take(5).collect { "${it.label}[f=${it.forbidden},s=${it.support}]" }.join(' '))
        result
    }

    /** Is {@code true} when asserting organism={@code organism} makes {@code term} taxon-inconsistent. */
    private boolean violates(String organism, String term) {
        def tc = checker.taxonConstraints
        Set<String> only = (tc.onlyInTaxon[term] ?: ([] as Set)) as Set
        Set<String> never = (tc.neverInTaxon[term] ?: ([] as Set)) as Set
        if (only.isEmpty() && never.isEmpty()) return false

        // Taxa forced true: the organism's lineage plus every required only_in
        // target's lineage (only_in is encoded as a unit clause per target).
        Set<String> trueSet = new HashSet<>(checker.ancestorsWithSelf(organism))
        only.each { trueSet.addAll(checker.ancestorsWithSelf(it)) }

        // (b) never_in target forced false but present in the true set -> conflict.
        for (n in never) if (trueSet.contains(n)) return true

        // (a) two taxa forced true that are explicitly disjoint -> conflict.
        for (a in trueSet) {
            if (!checker.disjointWith.containsKey(a)) continue
            for (b in checker.disjointWith[a]) if (trueSet.contains(b)) return true
        }
        false
    }

    /** Labelled lineage of a taxon, most specific first, capped for display. */
    private List<String> lineageOf(String taxon) {
        def anc = checker.ancestorsWithSelf(taxon)
        // order specific -> general by lineage depth (#ancestors), descending
        anc.toList().sort { -checker.ancestorsWithSelf(it).size() }
            .findAll { isRealTaxon(it) }
            .collect { labelFor(it) }
    }

    private static boolean isRealTaxon(String tax) { tax && !tax.contains('Union') }
}
