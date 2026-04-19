package gspa.ontology

import groovy.transform.CompileStatic

/**
 * Parse gapsmith's {@code seed_reactions.tsv} stoichiometry column
 * + {@code diffusion_mets.tsv} currency list, and optionally apply a
 * degree-based percentile threshold to detect high-connectivity
 * ("hub") metabolites to treat as currency.
 *
 * <p>Stoichiometry format produced by gapsmith:</p>
 * <pre>-1:cpd00001:0:0:"H2O";2:cpd00009:0:0:"Phosphate"</pre>
 *
 * <p>Semicolon-separated entries; each entry is
 * {@code <coef>:<compound>:<compartment>:<reserved>:"<name>"}.</p>
 */
@CompileStatic
class ReactionGraphLoader {

    /**
     * Load a ReactionGraph.
     *
     * @param reactionsTsv    seed_reactions.tsv
     * @param diffusionTsv    diffusion_mets.tsv (may be null)
     * @param currencyPercentile percentile threshold for automatic
     *        degree-based currency detection; 99.0 means metabolites in
     *        the top 1% by raw reaction-touching degree are treated as
     *        currency. Pass {@code 100.0} to disable.
     * @param ecAliasesTsv    optional seed_Enzyme_Class_Reactions_Aliases_unique.tsv —
     *        binds EC numbers to reaction sets so the suggester can
     *        resolve gap reactions by EC when the reactionId namespace
     *        differs (MetaCyc gap vs SEED graph).
     */
    static ReactionGraph load(File reactionsTsv, File diffusionTsv,
                              double currencyPercentile = 99.0d,
                              File ecAliasesTsv = null) {
        ReactionGraph g = new ReactionGraph()

        // 1. Explicit diffusion mets
        if (diffusionTsv != null && diffusionTsv.exists()) {
            diffusionTsv.withReader { r ->
                String line = r.readLine()  // header
                while ((line = r.readLine()) != null) {
                    String[] parts = line.split('\t')
                    if (parts.length >= 1 && parts[0].startsWith('cpd')) {
                        g.markCurrency(parts[0].trim())
                    }
                }
            }
        }

        // 2. Parse all reactions; first pass for graph build.
        // Also tally raw degree per metabolite for percentile threshold.
        Map<String, Integer> rawDegree = new HashMap<>()
        reactionsTsv.withReader { r ->
            String header = r.readLine()
            if (header == null) return
            List<String> cols = header.split('\t').collect { it.trim() }
            int idIdx = cols.indexOf('id')
            int stoichIdx = cols.indexOf('stoichiometry')
            if (stoichIdx < 0) stoichIdx = cols.indexOf('equation')
            int ecIdx = cols.indexOf('ec_numbers')
            if (ecIdx < 0) ecIdx = cols.indexOf('ec')
            if (idIdx < 0 || stoichIdx < 0) return

            String line
            while ((line = r.readLine()) != null) {
                String[] parts = line.split('\t', -1)
                if (parts.length <= stoichIdx) continue
                String rxnId = parts[idIdx].trim()
                String stoich = parts[stoichIdx]
                String ec = (ecIdx >= 0 && parts.length > ecIdx) ? parts[ecIdx].trim() : ''
                if (!rxnId || !stoich) continue
                Tuple2<Set<String>, Set<String>> sp = parseStoich(stoich)
                g.addReaction(new ReactionGraph.ReactionSpec(
                    rxnId: rxnId, ecNumber: ec,
                    substrates: sp.v1, products: sp.v2,
                ))
                for (String m : (Set<String>)(sp.v1 + sp.v2)) {
                    rawDegree.merge(m, 1, Integer::sum)
                }
            }
        }

        // 3. Degree-based currency detection (percentile of touched reactions)
        if (currencyPercentile < 100.0d && !rawDegree.isEmpty()) {
            List<Integer> sorted = rawDegree.values().sort()
            int idx = (int) Math.ceil(sorted.size() * currencyPercentile / 100.0d) - 1
            if (idx < 0) idx = 0
            if (idx >= sorted.size()) idx = sorted.size() - 1
            int threshold = sorted[idx]
            rawDegree.each { m, d ->
                if (d >= threshold) g.markCurrency(m)
            }
        }

        // 4. Optional EC aliases: binds EC → rxnId sets.
        if (ecAliasesTsv != null && ecAliasesTsv.exists()) {
            ecAliasesTsv.withReader { r ->
                String header = r.readLine()
                if (header == null) return
                List<String> cols = header.split('\t').collect { it.trim() }
                int msIdx = cols.findIndexOf { it.toLowerCase(Locale.ROOT).contains('ms id') } ?: 0
                int ecIdx = cols.findIndexOf { it.toLowerCase(Locale.ROOT).contains('external id') }
                if (ecIdx < 0) ecIdx = 2
                String line
                while ((line = r.readLine()) != null) {
                    String[] parts = line.split('\t', -1)
                    if (parts.length <= ecIdx) continue
                    String ec = parts[ecIdx].trim()
                    if (!ec || ec.contains('-')) continue   // skip partial ECs
                    String rxnSet = parts[msIdx]
                    if (!rxnSet) continue
                    for (String rxn : rxnSet.split('\\|')) {
                        rxn = rxn.trim()
                        if (rxn && g.reactions.containsKey(rxn)) {
                            g.bindEc(rxn, ec)
                        }
                    }
                }
            }
        }

        g.build()
        g
    }

    /** Returns (substrates, products). Coefficient sign decides side. */
    static Tuple2<Set<String>, Set<String>> parseStoich(String stoich) {
        Set<String> subs = new LinkedHashSet<>()
        Set<String> prods = new LinkedHashSet<>()
        if (stoich == null || stoich.isEmpty()) return new Tuple2(subs, prods)
        for (String entry : stoich.split(';')) {
            entry = entry.trim()
            if (entry.isEmpty()) continue
            // <coef>:<cpd>:<compartment>:<reserved>:"<name>"
            String[] fields = entry.split(':', 5)
            if (fields.length < 2) continue
            double coef
            try { coef = Double.parseDouble(fields[0]) } catch (NumberFormatException ignored) { continue }
            String cpd = fields[1].trim()
            if (cpd.isEmpty()) continue
            if (coef < 0) subs << cpd
            else if (coef > 0) prods << cpd
        }
        new Tuple2(subs, prods)
    }
}
