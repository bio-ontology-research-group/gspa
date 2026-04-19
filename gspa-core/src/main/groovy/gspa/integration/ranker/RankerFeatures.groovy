package gspa.integration.ranker

import gspa.integration.IntegrationState
import gspa.integration.crossgenome.ReactionLocusCatalog
import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import gspa.ontology.ReactionGraph

/**
 * Phase 12 M3: feature extractor for the learning-to-rank model.
 *
 * <p>Produces a flat vector of ~35 hand-engineered features per
 * (candidate protein {@code g}, gap reaction {@code R}) pair. The
 * features are organised by family:</p>
 * <ul>
 *   <li><b>density</b> — log D(g, R), nearby-anchor counts</li>
 *   <li><b>graph</b> — distance in reaction graph, EC class</li>
 *   <li><b>orthogroup</b> — M2 LR, CI width, cross-panel support</li>
 *   <li><b>size</b> — protein length, length-vs-EC-class-mean</li>
 *   <li><b>go</b> — IC of top GO, competing-GO strength</li>
 *   <li><b>local</b> — intergenic gap, strand consistency</li>
 * </ul>
 *
 * <p>Feature indices are stable (declared in {@link #FEATURE_NAMES});
 * adding a feature means appending to the end of the list and
 * re-training the model. Never remove a feature — set a zero sentinel
 * so the LightGBM txt model stays valid.</p>
 */
class RankerFeatures {

    /** Stable, ordered feature names. Index = column position. */
    public static final List<String> FEATURE_NAMES = [
        'log_density',
        'anchor_count_10kb',
        'anchor_count_40kb',
        'noisy_or_diversity',
        'graph_min_dist_to_R',
        'graph_min_dist_directed_in',
        'graph_min_dist_directed_out',
        'ec_digit_1',
        'ec_digit_2',
        'ec_digit_3',
        'ec_digit_4',
        'cg_log_lr',
        'cg_logci_halfwidth',
        'cg_n_sig_total',
        'cg_n_base_total',
        'cg_orthogroup_size_panel',
        'len_protein',
        'log_len_protein',
        'commitment_max',
        'commitment_count',
        'self_posterior',
        'go_ic_top',
        'go_overlap_nbr_count',
        'strand_consistency',
        'log_intergenic_gap',
        'is_dark',
        'operon_member',
        'rlc_score',
        'rlc_q',
        'rlc_rank',
    ]

    static final int DIM = FEATURE_NAMES.size()

    /** Extract features for one candidate. */
    static double[] extract(
            String proteinId,
            String gapReactionId,
            String gapEc,
            String gapGoTerm,
            ProteinLocus locus,
            double density,
            int anchorCount10kb,
            int anchorCount40kb,
            double noisyOrDiversity,
            int graphMinDist,
            int graphDirIn,
            int graphDirOut,
            ReactionLocusCatalog catalog,
            String orthogroupId,
            int orthogroupSizePanel,
            int proteinLength,
            double commitmentMax,
            int commitmentCount,
            double selfPosterior,
            double goIcTop,
            int goOverlapNbrCount,
            double strandConsistency,
            double logIntergenicGap,
            boolean isDark,
            boolean operonMember,
            double rlcScore,
            double rlcQ,
            int rlcRank) {
        double[] v = new double[DIM]
        int i = 0
        v[i++] = safeLog(density)
        v[i++] = anchorCount10kb
        v[i++] = anchorCount40kb
        v[i++] = noisyOrDiversity
        v[i++] = graphMinDist
        v[i++] = graphDirIn
        v[i++] = graphDirOut
        // EC digits (parse "1.2.3.4" or "1.-.-.-")
        def ecDigits = parseEc(gapEc)
        for (int d = 0; d < 4; d++) v[i++] = ecDigits[d]
        // Cross-genome
        if (catalog != null && orthogroupId != null) {
            def e = catalog.get(orthogroupId, gapReactionId)
            if (e != null) {
                v[i++] = e.logLR()
                v[i++] = e.logLRCiHalfWidth()
                v[i++] = e.nSigTotal
                v[i++] = e.nBaseTotal
            } else {
                v[i++] = 0.0d; v[i++] = Double.POSITIVE_INFINITY
                v[i++] = 0; v[i++] = catalog.panelSize
            }
        } else {
            v[i++] = 0.0d; v[i++] = Double.POSITIVE_INFINITY
            v[i++] = 0; v[i++] = 0
        }
        v[i++] = orthogroupSizePanel
        v[i++] = proteinLength
        v[i++] = safeLog((double) proteinLength)
        v[i++] = commitmentMax
        v[i++] = commitmentCount
        v[i++] = selfPosterior
        v[i++] = goIcTop
        v[i++] = goOverlapNbrCount
        v[i++] = strandConsistency
        v[i++] = logIntergenicGap
        v[i++] = isDark ? 1.0d : 0.0d
        v[i++] = operonMember ? 1.0d : 0.0d
        v[i++] = rlcScore
        v[i++] = rlcQ
        v[i++] = rlcRank
        if (i != DIM) {
            throw new IllegalStateException("Feature count mismatch: wrote ${i}, expected ${DIM}")
        }
        v
    }

    /** EC string → 4 integer digits with -1 for wildcards. */
    static int[] parseEc(String ec) {
        int[] out = [-1, -1, -1, -1]
        if (!ec) return out
        String trimmed = ec.startsWith('EC:') ? ec.substring(3) : ec
        String[] parts = trimmed.split('\\.')
        for (int i = 0; i < 4 && i < parts.length; i++) {
            try {
                out[i] = Integer.parseInt(parts[i])
            } catch (NumberFormatException ignored) {
                out[i] = -1
            }
        }
        out
    }

    private static double safeLog(double x) {
        x > 0 ? Math.log(x) : -30.0d
    }
}
