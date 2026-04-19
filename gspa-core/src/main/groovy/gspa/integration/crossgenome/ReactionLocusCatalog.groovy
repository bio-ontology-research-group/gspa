package gspa.integration.crossgenome

import groovy.transform.Canonical
import groovy.transform.CompileStatic

/**
 * Phase 12 M2 cross-genome catalog: per-(orthogroup C, reaction R)
 * conditional-LR statistic.
 *
 * <p>For each pair (C, R) in the panel we record:</p>
 * <ul>
 *   <li>{@code nSigWith}: number of genomes where R-signature is present
 *       AND an ortholog from C lies in any R-signature window;</li>
 *   <li>{@code nSigTotal}: number of genomes where R-signature is present;</li>
 *   <li>{@code nBaseWith}: number of genomes containing any member of C;</li>
 *   <li>{@code nBaseTotal}: total panel size.</li>
 * </ul>
 *
 * <p>The conditional LR statistic is:</p>
 * <pre>
 *   LR(C, R) = [p(C in window | R present)] / [p(C present | baseline)]
 *            = (nSigWith / nSigTotal) / (nBaseWith / nBaseTotal)
 * </pre>
 *
 * <p>We compute the LR and a Jeffreys-prior beta-binomial 90% credible
 * interval at lookup time so the scorer can drop low-confidence LRs
 * (CI overlapping 1.0).</p>
 */
@CompileStatic
class ReactionLocusCatalog {

    @Canonical
    static class Entry {
        String orthogroupId
        String reactionId
        int nSigWith
        int nSigTotal
        int nBaseWith
        int nBaseTotal

        double rateSig() {
            nSigTotal == 0 ? 0.0d : (double) nSigWith / (double) nSigTotal
        }
        double rateBase() {
            nBaseTotal == 0 ? 0.0d : (double) nBaseWith / (double) nBaseTotal
        }
        double logLR(double epsilon = 1e-3d) {
            double num = Math.max(epsilon, rateSig())
            double den = Math.max(epsilon, rateBase())
            Math.log(num / den)
        }
        double lr(double epsilon = 1e-3d) {
            Math.exp(logLR(epsilon))
        }
        /**
         * 90% beta-binomial credible-interval half-width on logLR,
         * using Jeffreys prior α=β=0.5 on the two rates independently.
         * Returns the absolute log-space half-width; caller drops the
         * entry if this crosses 0 (CI includes LR=1).
         */
        double logLRCiHalfWidth(double epsilon = 1e-3d) {
            // Wilson-style interval on each rate; approximate via variance of logit.
            double pSig = Math.max(epsilon, rateSig())
            double pBase = Math.max(epsilon, rateBase())
            double varSig = pSig * (1 - pSig) / Math.max(1, nSigTotal)
            double varBase = pBase * (1 - pBase) / Math.max(1, nBaseTotal)
            // delta-method: Var(log p) ≈ Var(p) / p^2
            double varLog = varSig / (pSig * pSig) + varBase / (pBase * pBase)
            1.645d * Math.sqrt(varLog)   // 90% CI
        }
    }

    /** total number of genomes in the panel. */
    int panelSize = 0

    /** (orthogroupId, reactionId) -> Entry */
    final Map<Tuple2<String, String>, Entry> entries = new LinkedHashMap<>()

    /** Per-orthogroup baseline presence count (computed once across genomes). */
    final Map<String, Integer> baselinePresence = new LinkedHashMap<>()

    Entry get(String orthogroupId, String reactionId) {
        entries.get(new Tuple2(orthogroupId, reactionId))
    }

    void put(Entry e) {
        entries.put(new Tuple2(e.orthogroupId, e.reactionId), e)
    }

    int size() { entries.size() }

    /** Write a compact TSV: orthogroup_id, reaction_id, n_sig_with, n_sig_total, n_base_with, n_base_total. */
    void writeTo(File out) {
        out.withWriter { w ->
            w.writeLine('# ReactionLocusCatalog')
            w.writeLine('# panel_size=' + panelSize)
            w.writeLine('orthogroup_id\treaction_id\tn_sig_with\tn_sig_total\tn_base_with\tn_base_total')
            entries.values().each { Entry e ->
                w.writeLine([e.orthogroupId, e.reactionId,
                             e.nSigWith.toString(), e.nSigTotal.toString(),
                             e.nBaseWith.toString(), e.nBaseTotal.toString()].join('\t'))
            }
        }
    }

    static ReactionLocusCatalog readFrom(File tsv) {
        ReactionLocusCatalog cat = new ReactionLocusCatalog()
        tsv.withReader { r ->
            String line
            while ((line = r.readLine()) != null) {
                if (line.startsWith('#')) {
                    if (line.startsWith('# panel_size=')) {
                        cat.panelSize = Integer.parseInt(line.substring('# panel_size='.length()).trim())
                    }
                    continue
                }
                if (line.startsWith('orthogroup_id')) continue
                String[] parts = line.split('\t')
                if (parts.length < 6) continue
                cat.put(new Entry(
                    orthogroupId: parts[0],
                    reactionId: parts[1],
                    nSigWith: Integer.parseInt(parts[2]),
                    nSigTotal: Integer.parseInt(parts[3]),
                    nBaseWith: Integer.parseInt(parts[4]),
                    nBaseTotal: Integer.parseInt(parts[5]),
                ))
            }
        }
        cat
    }
}
