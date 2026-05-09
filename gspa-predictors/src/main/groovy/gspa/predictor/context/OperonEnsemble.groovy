package gspa.predictor.context

import gspa.model.Genome
import gspa.model.Protein
import gspa.model.Strand
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Ensemble operon caller built on a small collection of independent
 * predictors, each scoring per adjacent gene-pair "co-operonic"
 * candidacy. Pair scores are combined with Noisy-OR over per-predictor
 * sensitivity θ, then transitive closure produces the final operons.
 *
 * <h3>Why an ensemble</h3>
 * Operons are not directly observable in a single genome. The classic
 * intergenic-distance + same-strand rule (the existing
 * {@link OperonPredictor}) is the most-used signal but has its blind spots:
 * very-short overlapping genes look co-operonic when they are not, and
 * functionally related genes separated by a hundred bp often miss the
 * cutoff. Combining several signals as independent evidence gives
 * sharper, calibrated co-operonic posteriors and lets us report
 * *confidence* per operon, not just a rule-vote.
 *
 * <h3>Predictors in this v1 ensemble</h3>
 * <ol>
 *   <li><b>Distance</b> — same strand, intergenic gap ≤ 300 bp.
 *       (Identical to {@link OperonPredictor}; the de-facto baseline.)</li>
 *   <li><b>Strict</b> — same strand, intergenic gap ≤ 50 bp. Higher
 *       precision; a "definitely co-transcribed" call.</li>
 *   <li><b>Functional</b> — same strand, gap ≤ 1000 bp, AND the two genes
 *       share at least one annotated GO biological-process term (excluding
 *       the very-broad terms in {@link #BP_DENYLIST}). Picks up
 *       longer-spaced operons whose genes belong to the same pathway.</li>
 * </ol>
 *
 * Default per-predictor sensitivity θ:
 * Distance = 0.85, Strict = 0.95, Functional = 0.70.
 * These are uncalibrated priors; a future Phase will tune them on a
 * curated operon ground-truth set (e.g. ODB4, RegulonDB).
 *
 * <h3>Output</h3>
 * {@link Operon} as before, with three extra fields the JSON serialiser
 * picks up: {@code supportSet} (which predictor names voted for the
 * pair), {@code minPairPosterior} (weakest co-operonic edge in the
 * operon — the limiting confidence), and {@code meanPairPosterior}.
 */
class OperonEnsemble {

    private static final Logger log = LoggerFactory.getLogger(OperonEnsemble)

    /** Per-predictor reliability θ ∈ (0, 1) used in the Noisy-OR step. */
    Map<String, Double> reliability = [
        distance:   0.85d,
        strict:     0.95d,
        functional: 0.70d,
    ]

    /** Posterior threshold for accepting a pair as co-operonic. */
    double pairPosteriorThreshold = 0.50d

    /** Min operon size (gene count) to emit an Operon. */
    int minOperonSize = 2

    /** Functional predictor: max intergenic distance (bp) at which functional
     *  agreement is allowed to imply co-operonicity. */
    int functionalMaxGap = 1000

    /** Distance predictor max gap (bp). */
    int distanceMaxGap = 300

    /** Strict predictor max gap (bp). */
    int strictMaxGap = 50

    /** GO BP terms too broad to count as "shared function" evidence. */
    static final Set<String> BP_DENYLIST = [
        'GO:0008150', 'GO:0008152', 'GO:0009987', 'GO:0050896', 'GO:0065007',
        'GO:0050789', 'GO:0050794', 'GO:0044238', 'GO:0044237', 'GO:0071704',
        'GO:0006807', 'GO:0019222',
    ] as Set

    /**
     * Run the ensemble across all contigs in the genome and return the
     * detected operons. Each Operon carries supportSet / pair-posterior
     * stats (see {@link Operon}).
     */
    List<Operon> detect(Genome genome) {
        List<Operon> operons = []
        for (def contig : genome.contigs) {
            List<Protein> sorted = contig.proteins.sort(false) { it.start }
            if (sorted.size() < 2) continue
            // Compute pair-level votes + posteriors for adjacent gene pairs.
            List<PairCall> pairs = []
            for (int i = 0; i < sorted.size() - 1; i++) {
                Protein a = sorted[i], b = sorted[i + 1]
                pairs << callPair(a, b)
            }
            // Walk pairs and group adjacent co-operonic pairs into operons.
            List<Protein> current = [sorted[0]]
            List<PairCall> currentEdges = []
            Set<String> currentSupport = new LinkedHashSet<>()
            for (int i = 0; i < pairs.size(); i++) {
                PairCall pc = pairs[i]
                if (pc.posterior >= pairPosteriorThreshold) {
                    current << sorted[i + 1]
                    currentEdges << pc
                    currentSupport.addAll(pc.support)
                } else {
                    if (current.size() >= minOperonSize) {
                        operons << buildOperon(contig.id, current, currentEdges, currentSupport)
                    }
                    current = [sorted[i + 1]]
                    currentEdges = []
                    currentSupport = new LinkedHashSet<>()
                }
            }
            if (current.size() >= minOperonSize) {
                operons << buildOperon(contig.id, current, currentEdges, currentSupport)
            }
        }
        log.info("OperonEnsemble: ${operons.size()} operons; ${operons.sum { it.size } ?: 0} CDS in operons")
        operons
    }

    /**
     * Score a single adjacent gene-pair across all 3 predictors and combine
     * via Noisy-OR. Returns the per-predictor vote set + the combined
     * posterior.
     */
    PairCall callPair(Protein a, Protein b) {
        Set<String> support = new LinkedHashSet<>()
        boolean sameStrand = a.strand != null && b.strand != null && a.strand.isSameStrand(b.strand)
        int gap = intergenicDistance(a, b)
        // Distance predictor
        if (sameStrand && gap <= distanceMaxGap) support << 'distance'
        // Strict predictor
        if (sameStrand && gap <= strictMaxGap) support << 'strict'
        // Functional predictor
        if (sameStrand && gap <= functionalMaxGap && sharesInformativeBp(a, b)) {
            support << 'functional'
        }
        // Noisy-OR over per-predictor sensitivities. p(co-operonic | votes) =
        // 1 - prod_{k in votes}(1 - θ_k). Predictors that didn't vote
        // contribute 0 evidence (factor 1).
        double oneMinusProd = 1.0d
        for (String pred : support) {
            double theta = reliability.getOrDefault(pred, 0.5d)
            oneMinusProd *= (1.0d - theta)
        }
        double posterior = 1.0d - oneMinusProd
        new PairCall(a: a, b: b, support: support, posterior: posterior, gap: gap)
    }

    private boolean sharesInformativeBp(Protein a, Protein b) {
        Set<String> aBp = informativeBpTerms(a)
        if (aBp.isEmpty()) return false
        Set<String> bBp = informativeBpTerms(b)
        if (bBp.isEmpty()) return false
        for (String t : aBp) if (bBp.contains(t)) return true
        return false
    }

    private Set<String> informativeBpTerms(Protein p) {
        Set<String> out = []
        if (p.annotations == null) return out
        for (def ann : p.annotations.goAnnotations()) {
            if (ann.goAspect != 'BP' && ann.goAspect != null) continue
            if (ann.value == null) continue
            if (BP_DENYLIST.contains(ann.value)) continue
            out << ann.value
        }
        out
    }

    private Operon buildOperon(String contigId, List<Protein> genes,
                                List<PairCall> edges, Set<String> support) {
        def op = new Operon(contigId: contigId, genes: new ArrayList<>(genes))
        op.supportSet = new LinkedHashSet<>(support)
        if (!edges.isEmpty()) {
            double minP = 1.0d, sumP = 0.0d
            for (def e : edges) {
                if (e.posterior < minP) minP = e.posterior
                sumP += e.posterior
            }
            op.minPairPosterior = minP
            op.meanPairPosterior = sumP / edges.size()
        } else {
            op.minPairPosterior = 0.0d
            op.meanPairPosterior = 0.0d
        }
        op
    }

    private static int intergenicDistance(Protein a, Protein b) {
        if (a.end < b.start) return b.start - a.end - 1
        if (b.end < a.start) return a.start - b.end - 1
        return -(Math.min(a.end, b.end) - Math.max(a.start, b.start) + 1)
    }

    /**
     * Per-pair call record produced by {@link #callPair}.
     */
    static class PairCall {
        Protein a, b
        Set<String> support
        double posterior
        int gap
    }
}
