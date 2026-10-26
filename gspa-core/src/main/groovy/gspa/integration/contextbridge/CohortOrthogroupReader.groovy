package gspa.integration.contextbridge

import gspa.integration.cohort.CohortInput

/**
 * Reads family (orthogroup / subfamily) membership for the cohort engine.
 *
 * <p>Two formats, both strict about malformed rows:</p>
 * <ul>
 *   <li><b>Cohort file</b> (genome-aware, header required):
 *       {@code genome_id  protein_id  orthogroup_id  [subfamily_id]  [weight]}. A protein with several
 *       rows has soft membership; {@code -} or an empty field means "no subfamily" / weight 1.</li>
 *   <li><b>Per-genome file</b> (the existing {@code --orthogroups} format):
 *       {@code protein_id  cluster_id}, where cluster ids are cohort-wide (an MMseqs2 easy-cluster over
 *       all genomes). Rows for proteins that are not proteins of this genome are other genomes'
 *       members and are counted, not guessed at.</li>
 * </ul>
 *
 * <p>Protein ids are FAA sequence ids, the same ones the layout uses. A {@code namespace:} prefix on a
 * cluster id is stripped exactly as in the single-genome loader so both routes name a family the
 * same way.</p>
 */
class CohortOrthogroupReader {

    static class Rows {
        final Map<String, List<CohortInput.FamilyAssignment>> byGenome = new LinkedHashMap<>()
        int rows
        int rowsOfGenomesOutsideCohort
        int rowsOfProteinsOutsideGenome
        final Set<String> genomesOutsideCohort = new TreeSet<>()
    }

    static String normalizeFamily(String raw) {
        String rep = raw.trim()
        int colon = rep.indexOf(':')
        colon >= 0 ? rep.substring(colon + 1) : rep
    }

    /** Genome-aware cohort file. Rows of genomes outside {@code cohortGenomes} are counted and skipped. */
    static Rows readCohortFile(File file, Set<String> cohortGenomes) {
        if (file == null || !file.isFile()) throw new IllegalArgumentException("cohort orthogroups file does not exist: ${file}")
        Rows out = new Rows()
        int lineNo = 0
        boolean header = false
        file.eachLine('UTF-8') { String raw ->
            lineNo++
            String line = raw.trim()
            if (!line || line.startsWith('#')) return
            String[] f = raw.split('\t', -1)
            if (!header) {
                header = true
                if (f.length < 3 || f[0].trim().toLowerCase() != 'genome_id' || f[1].trim().toLowerCase() != 'protein_id' ||
                    !(f[2].trim().toLowerCase() in ['orthogroup_id', 'orthogroup'])) {
                    throw new IllegalArgumentException("${file}: the first line must be the header 'genome_id<TAB>protein_id<TAB>orthogroup_id[<TAB>subfamily_id<TAB>weight]'")
                }
                return
            }
            if (f.length < 3 || f.length > 5) throw new IllegalArgumentException("${file}:${lineNo}: expected 3 to 5 tab-separated fields, got ${f.length}")
            String genome = f[0].trim(), protein = f[1].trim(), group = normalizeFamily(f[2])
            if (!genome || !protein || !group) throw new IllegalArgumentException("${file}:${lineNo}: genome_id, protein_id and orthogroup_id are required")
            String sub = f.length > 3 && f[3].trim() && f[3].trim() != '-' ? f[3].trim() : null
            double weight = 1.0d
            if (f.length > 4 && f[4].trim() && f[4].trim() != '-') {
                try { weight = Double.parseDouble(f[4].trim()) } catch (NumberFormatException e) {
                    throw new IllegalArgumentException("${file}:${lineNo}: weight '${f[4]}' is not a number")
                }
            }
            out.rows++
            if (!cohortGenomes.contains(genome)) {
                out.rowsOfGenomesOutsideCohort++
                out.genomesOutsideCohort << genome
                return
            }
            out.byGenome.computeIfAbsent(genome) { new ArrayList<CohortInput.FamilyAssignment>() } <<
                new CohortInput.FamilyAssignment(proteinId: protein, orthogroupId: group, subfamilyId: sub, weight: weight)
        }
        if (!header) throw new IllegalArgumentException("${file} is empty")
        out
    }

    /** The existing two-column per-genome file; only rows of proteins in {@code knownProteins} belong to this genome. */
    static Rows readPerGenomeFile(File file, String genomeId, Set<String> knownProteins) {
        if (file == null || !file.isFile()) throw new IllegalArgumentException("orthogroups file does not exist: ${file}")
        Rows out = new Rows()
        int lineNo = 0
        file.eachLine('UTF-8') { String raw ->
            lineNo++
            String line = raw.trim()
            if (!line || line.startsWith('#')) return
            String[] f = raw.split('\t', -1)
            if (f.length < 2 || !f[0].trim() || !f[1].trim()) {
                throw new IllegalArgumentException("${file}:${lineNo}: expected 'protein_id<TAB>cluster_id'")
            }
            out.rows++
            String protein = f[0].trim()
            if (knownProteins != null && !knownProteins.contains(protein)) {
                out.rowsOfProteinsOutsideGenome++
                return
            }
            out.byGenome.computeIfAbsent(genomeId) { new ArrayList<CohortInput.FamilyAssignment>() } <<
                new CohortInput.FamilyAssignment(proteinId: protein, orthogroupId: normalizeFamily(f[1]), weight: 1.0d)
        }
        out
    }
}
