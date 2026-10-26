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
 * <p>Rates use Jeffreys smoothing. The log-ratio uncertainty is a
 * delta-method approximation using beta posterior variances, not an exact
 * credible interval. Small panels require independent calibration.</p>
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
            (nSigWith + 0.5d) / (nSigTotal + 1.0d)
        }
        double rateBase() {
            (nBaseWith + 0.5d) / (nBaseTotal + 1.0d)
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
         * Approximate 90% log-ratio interval half-width, using beta posterior
         * variances with Jeffreys prior α=β=0.5 on the two rates independently.
         * Returns the absolute log-space half-width; caller drops the
         * entry if this crosses 0 (CI includes LR=1).
         */
        double logLRCiHalfWidth(double epsilon = 1e-3d) {
            double pSig = Math.max(epsilon, rateSig())
            double pBase = Math.max(epsilon, rateBase())
            double varSig = pSig * (1 - pSig) / (nSigTotal + 2d)
            double varBase = pBase * (1 - pBase) / (nBaseTotal + 2d)
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
        if (e.nSigTotal < 0 || e.nBaseTotal < 0 || e.nSigWith < 0 || e.nBaseWith < 0 ||
                e.nSigWith > e.nSigTotal || e.nBaseWith > e.nBaseTotal) {
            throw new IllegalArgumentException('Catalog counts must satisfy 0 <= occurrences <= observations')
        }
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
