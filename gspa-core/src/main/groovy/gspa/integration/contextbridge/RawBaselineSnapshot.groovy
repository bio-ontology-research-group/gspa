package gspa.integration.contextbridge

import gspa.integration.EvidenceClaim
import gspa.integration.cohort.CohortInput
import gspa.model.AnnotationType

import java.security.MessageDigest

/**
 * The pre-cohort baseline snapshot of a batch: for every genome, the direct EC calls with the
 * strongest RAW predictor score per (protein, EC). It is the Stage 0 evidence only.
 *
 * <p><b>Labels.</b> This is the {@code stage0_raw_max_score} baseline. It is not the legacy prior
 * baseline (the posterior after {@code evidence_integration}, which the audit ledger reports as
 * {@code legacy_prior_baseline}); posteriors, priors, floors and calibrated probabilities are never
 * read, so the cohort engine cannot consume the output of the heuristic stages it will be compared
 * with, nor its own. The score is an uncalibrated predictor score in [0, 1]; the cohort engine only
 * uses it to decide whether a neighbouring candidate is a confident anchor.</p>
 *
 * <p>The snapshot id is a SHA-256 over the sorted (genome, protein, EC, score bits) lines, so the same
 * evidence gives the same id whatever the order of genomes or claims.</p>
 */
class RawBaselineSnapshot {

    static final String LABEL = 'stage0_raw_max_score'
    /** Provenance of every baseline candidate (must not look like cohort/feedback output). */
    static final String PROVENANCE = 'stage0-raw-max-score'

    /** One genome's baseline. */
    static class GenomeBaseline {
        String genomeId
        final List<CohortInput.EcCall> calls = new ArrayList<>()
        /** normalized EC -> function keys the claims used, per protein: protein + TAB + EC. */
        final Map<String, Set<String>> functionKeys = new LinkedHashMap<>()
        int nClaims
        int nEcClaims
        int nNonEcClaimsIgnored
        int nInvalidScoreIgnored
        int nPartialEcIgnored
    }

    final Map<String, GenomeBaseline> genomes = new TreeMap<>()
    String id

    /** Strip the {@code EC:} prefix; the reaction graph keys ECs bare. Null for blank input. */
    static String normalizeEc(String ec) {
        if (ec == null) return null
        String e = ec.trim()
        if (e.regionMatches(true, 0, 'EC:', 0, 3)) e = e.substring(3)
        e ?: null
    }

    static GenomeBaseline forGenome(String genomeId, Collection<EvidenceClaim> claims) {
        GenomeBaseline out = new GenomeBaseline(genomeId: genomeId)
        Map<String, Double> best = new TreeMap<>()
        for (EvidenceClaim c : claims) {
            out.nClaims++
            if (c.isQualityObservation() || c.functionType != AnnotationType.EC) {
                out.nNonEcClaimsIgnored++
                continue
            }
            out.nEcClaims++
            String ec = normalizeEc(c.functionId)
            if (ec == null || ec.contains('-')) {
                out.nPartialEcIgnored++
                continue
            }
            if (!c.hasFunctionalScore() || !Double.isFinite(c.rawScore) || c.rawScore < 0d || c.rawScore > 1d) {
                out.nInvalidScoreIgnored++
                continue
            }
            String key = c.proteinId + '\t' + ec
            Double old = best.get(key)
            if (old == null || c.rawScore > old.doubleValue()) best.put(key, c.rawScore)
            out.functionKeys.computeIfAbsent(key) { new TreeSet<String>() }.add(c.functionKey())
        }
        best.each { String key, Double score ->
            String[] p = key.split('\t')
            out.calls << new CohortInput.EcCall(genomeId: genomeId, proteinId: p[0], ecNumber: p[1], confidence: score,
                provenance: PROVENANCE)
        }
        out
    }

    RawBaselineSnapshot add(GenomeBaseline baseline) {
        if (genomes.containsKey(baseline.genomeId)) throw new IllegalArgumentException("duplicate genome in baseline snapshot: ${baseline.genomeId}")
        genomes.put(baseline.genomeId, baseline)
        id = null
        this
    }

    String snapshotId() {
        if (id != null) return id
        MessageDigest md = MessageDigest.getInstance('SHA-256')
        md.update((LABEL + '/1').getBytes('UTF-8'))
        genomes.each { String g, GenomeBaseline b ->
            md.update(("\ngenome\t" + g).getBytes('UTF-8'))
            b.calls.each { CohortInput.EcCall c ->
                md.update(("\n" + c.proteinId + '\t' + c.ecNumber + '\t' + Long.toHexString(Double.doubleToLongBits(c.confidence))).getBytes('UTF-8'))
            }
        }
        id = LABEL + ':' + md.digest().collect { String.format('%02x', it & 0xff) }.join().substring(0, 24)
        id
    }
}
