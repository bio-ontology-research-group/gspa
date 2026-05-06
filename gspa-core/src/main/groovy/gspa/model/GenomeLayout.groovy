package gspa.model

import groovy.transform.CompileStatic

/**
 * A genome's protein-on-contig layout.
 *
 * <p>Built by {@link GenomeLayoutLoader} from a layout TSV (protein_id,
 * contig, start, end, strand). Used by the Phase 12 Reaction-Local
 * Genomic Context framework for kernel-smoothed density queries and
 * by the cross-genome catalog builder for window extraction.</p>
 *
 * <p>Per-contig lists are kept sorted by midpoint for cheap windowed
 * queries. {@link #windowAround} is the hot path used per gap reaction
 * per candidate locus.</p>
 */
@CompileStatic
class GenomeLayout {

    /** genomeId carried on every ProteinLocus (may be empty / null). */
    String genomeId

    /** All known loci keyed by protein id. */
    final Map<String, ProteinLocus> byProtein = new LinkedHashMap<>()

    /** Per-contig loci sorted by midpoint ascending. */
    final Map<String, List<ProteinLocus>> byContig = new LinkedHashMap<>()

    GenomeLayout() {}

    GenomeLayout(String genomeId) {
        this.genomeId = genomeId
    }

    void add(ProteinLocus locus) {
        byProtein[locus.proteinId] = locus
        byContig.computeIfAbsent(locus.contig, { new ArrayList<ProteinLocus>() }) << locus
    }

    /** Call once after all add()s to finalise sort orders. */
    void finishLoading() {
        byContig.values().each { list ->
            list.sort { a, b -> Integer.compare(a.midpoint(), b.midpoint()) }
        }
    }

    int size() { byProtein.size() }

    ProteinLocus get(String proteinId) { byProtein[proteinId] }

    List<ProteinLocus> contig(String contigId) {
        byContig.getOrDefault(contigId, Collections.<ProteinLocus>emptyList())
    }

    /**
     * Proteins on {@code contig} whose midpoint falls within
     * {@code [center - halfWidth, center + halfWidth]}.
     * O(log N + W) via binary search + scan.
     */
    List<ProteinLocus> windowAround(String contigId, int center, int halfWidth) {
        List<ProteinLocus> loci = contig(contigId)
        if (loci.isEmpty() || halfWidth <= 0) return []
        int lo = center - halfWidth
        int hi = center + halfWidth
        int idx = binarySearchMidpoint(loci, lo)
        List<ProteinLocus> out = []
        for (int i = idx; i < loci.size(); i++) {
            int m = loci[i].midpoint()
            if (m > hi) break
            if (m >= lo) out << loci[i]
        }
        out
    }

    /** Index of first locus with midpoint >= target. */
    private static int binarySearchMidpoint(List<ProteinLocus> loci, int target) {
        int lo = 0, hi = loci.size()
        while (lo < hi) {
            int mid = (lo + hi).intdiv(2)
            if (loci[mid].midpoint() < target) lo = mid + 1
            else hi = mid
        }
        lo
    }

    /**
     * Intergenic-gap signal: distance (in bp) from {@code probe}'s CDS
     * boundary to its same-strand nearest neighbour on the same contig.
     * Returns {@link Integer#MAX_VALUE} if no same-strand neighbour is
     * found within a pragmatic upper bound (100 kb).
     */
    int sameStrandIntergenicGap(ProteinLocus probe, int maxBp = 100000) {
        List<ProteinLocus> loci = contig(probe.contig)
        if (loci == null || loci.isEmpty()) return Integer.MAX_VALUE
        // Bracket search: find probe's index via binary search, then walk
        // out in both directions until we leave maxBp or find a same-strand
        // neighbour. O(log N + k) where k is typically tiny.
        int idx = binarySearchMidpoint(loci, probe.midpoint())
        int best = Integer.MAX_VALUE
        // Walk right
        for (int i = idx; i < loci.size(); i++) {
            ProteinLocus other = loci[i]
            if (other.proteinId == probe.proteinId) continue
            int dist = other.start - probe.end
            if (dist > maxBp) break
            if (!other.strand.isSameStrand(probe.strand)) continue
            int gap = dist < 0 ? 0 : dist
            if (gap < best) best = gap
        }
        // Walk left
        for (int i = idx - 1; i >= 0; i--) {
            ProteinLocus other = loci[i]
            if (other.proteinId == probe.proteinId) continue
            int dist = probe.start - other.end
            if (dist > maxBp) break
            if (!other.strand.isSameStrand(probe.strand)) continue
            int gap = dist < 0 ? 0 : dist
            if (gap < best) best = gap
        }
        best
    }
}
