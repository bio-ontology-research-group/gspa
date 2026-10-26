package gspa.integration.contextbridge

/**
 * What the JVM needs to know about an SBML-FBC model to build valid experiment hypotheses: the gene
 * products (id and label), the reactions with their gene-product references, and, if the model came with the
 * refiner's {@code reactions.tsv}, each reaction's EC numbers, candidate genes and construction basis.
 * Read with a non-validating, DOCTYPE-free parser; nothing is solved or interpreted here.
 *
 * <p>The refiner builds gene-protein-reaction rules FROM the EC evidence it is given, so a gene being in
 * the model's rule for a reaction is the candidate assignment itself, not an independent observation. The
 * {@code basis} column ({@code flux-certified repair} vs {@code core/support}) is kept so a reaction that was
 * added to make growth feasible is never counted as independent evidence for a growth objective.</p>
 */
class SbmlModelInfo {

    /** One row of the refiner's reactions.tsv. */
    static class ReactionRow {
        String reactionId
        Set<String> ecs = new LinkedHashSet<>()
        Set<String> candidateGenes = new LinkedHashSet<>()
        String basis
    }

    File model
    String sha256
    /** gene product label -> id. */
    final Map<String, String> geneByLabel = new LinkedHashMap<>()
    final Set<String> geneIds = new LinkedHashSet<>()
    final Set<String> reactionIds = new LinkedHashSet<>()
    /** reaction -> gene product ids it references. */
    final Map<String, Set<String>> reactionGenes = new LinkedHashMap<>()
    /** reaction -> refiner row; empty when the model has no refiner table. */
    final Map<String, ReactionRow> table = new LinkedHashMap<>()
    /** species id -> exact compound cross-references ({@code namespace:id}) from its SBML annotation. */
    final Map<String, Set<String>> speciesXrefs = new LinkedHashMap<>()
    /** reaction -> species it consumes or produces. */
    final Map<String, Set<String>> reactionSpecies = new LinkedHashMap<>()

    /** Namespaces that identify a chemical compound independently of its compartment. */
    static final Set<String> COMPOUND_NAMESPACES = ['chebi', 'kegg.compound', 'bigg.metabolite', 'metanetx.chemical',
        'seed.compound', 'hmdb', 'pubchem.compound', 'inchikey'] as Set

    static SbmlModelInfo read(File model, File reactionsTsv = null) {
        if (model == null || !model.isFile()) throw new IllegalArgumentException("model file does not exist: ${model}")
        def info = new SbmlModelInfo(model: model, sha256: SidecarRun.sha256(model))
        def factory = javax.xml.parsers.DocumentBuilderFactory.newInstance()
        factory.setNamespaceAware(false)                 // prefixed names are kept as written (fbc:geneProduct)
        factory.setFeature('http://apache.org/xml/features/disallow-doctype-decl', true)
        factory.setFeature('http://xml.org/sax/features/external-general-entities', false)
        factory.setXIncludeAware(false)
        org.w3c.dom.Document doc
        try {
            doc = factory.newDocumentBuilder().parse(model)
        } catch (Exception e) {
            throw new IllegalArgumentException("model ${model} is not readable SBML: ${e.message}", e)
        }
        def attr = { org.w3c.dom.Element e, String... names -> names.collect { e.getAttribute(it) }.find { it } ?: '' }
        def elements = { org.w3c.dom.Node scope, String tag ->
            org.w3c.dom.NodeList l = scope instanceof org.w3c.dom.Document ? scope.getElementsByTagName(tag) : (scope as org.w3c.dom.Element).getElementsByTagName(tag)
            (0..<l.length).collect { l.item(it) as org.w3c.dom.Element }
        }
        elements(doc, 'fbc:geneProduct').each { org.w3c.dom.Element gp ->
            String id = attr(gp, 'fbc:id', 'id')
            String label = attr(gp, 'fbc:label', 'label')
            if (id) {
                info.geneIds << id
                if (label) info.geneByLabel.put(label, id)
            }
        }
        elements(doc, 'species').each { org.w3c.dom.Element sp ->
            String sid = sp.getAttribute('id')
            if (!sid) return
            elements(sp, 'rdf:li').each { org.w3c.dom.Element li ->
                String x = compoundXref(li.getAttribute('rdf:resource'))
                if (x != null) info.speciesXrefs.computeIfAbsent(sid) { new TreeSet<String>() } << x
            }
        }
        elements(doc, 'reaction').each { org.w3c.dom.Element rx ->
            String rid = rx.getAttribute('id')
            if (!rid) return
            info.reactionIds << rid
            Set<String> genes = new LinkedHashSet<>()
            elements(rx, 'fbc:geneProductRef').each { org.w3c.dom.Element ref -> genes << ref.getAttribute('fbc:geneProduct') }
            info.reactionGenes.put(rid, genes)
            Set<String> species = new TreeSet<>()
            elements(rx, 'speciesReference').each { org.w3c.dom.Element ref -> if (ref.getAttribute('species')) species << ref.getAttribute('species') }
            info.reactionSpecies.put(rid, species)
            elements(rx, 'rdf:li').each { org.w3c.dom.Element li ->
                def m = (li.getAttribute('rdf:resource') =~ /(?:ec-code[\/:]|enzyme[\/:]|ec[\/:])(\d+\.\d+\.\d+\.\d+)$/)
                if (m.find()) info.annotationEcs.computeIfAbsent(rid) { new LinkedHashSet<String>() } << m.group(1)
            }
        }
        if (info.reactionIds.isEmpty()) throw new IllegalArgumentException("model ${model} has no reactions")
        File tsv = reactionsTsv ?: new File(model.absoluteFile.parentFile, 'reactions.tsv')
        if (tsv.isFile()) info.readTable(tsv, reactionsTsv != null)
        info
    }

    private void readTable(File tsv, boolean required) {
        List<String> lines = tsv.readLines('UTF-8')
        if (lines.isEmpty()) return
        List<String> header = lines[0].split('\t', -1).collect { it.trim() }
        int id = header.indexOf('reaction_id'), ec = header.indexOf('EC'), genes = header.indexOf('candidate_genes'), basis = header.indexOf('basis')
        if (id < 0 || ec < 0 || genes < 0 || basis < 0) {
            if (required) throw new IllegalArgumentException("${tsv} lacks one of reaction_id, EC, candidate_genes, basis")
            return                                     // some other reactions.tsv: not the refiner's table
        }
        lines.drop(1).each { String line ->
            if (!line.trim()) return
            String[] f = line.split('\t', -1)
            if (f.length <= Math.max(Math.max(id, ec), Math.max(genes, basis))) return
            def row = new ReactionRow(reactionId: f[id].trim(), basis: f[basis].trim())
            f[ec].split(';').collect { it.trim() }.findAll { it }.each { row.ecs << it.replaceFirst(/^EC:/, '') }
            f[genes].split(';').collect { it.trim() }.findAll { it }.each { row.candidateGenes << it }
            table.put(row.reactionId, row)
        }
    }

    boolean hasTable() { !table.isEmpty() }

    /** reaction -> EC numbers named by the SBML annotation of the reaction (identifiers.org / ec-code URIs). */
    final Map<String, Set<String>> annotationEcs = new LinkedHashMap<>()

    /**
     * Model reactions validated for an EC number, with the source of each validation: the refiner's table
     * ({@code reactions_tsv}), the model's own annotation ({@code sbml_annotation}) or the reaction graph when
     * its reaction id IS a reaction id of the model ({@code reaction_graph_id}). No gene is consulted: an EC call
     * on a protein that is missing from the model's gene rules still maps to the reactions that EC names.
     * Empty = no validated mapping (the candidate is unavailable for experiments, never guessed).
     */
    Map<String, Set<String>> reactionsForEc(String ec, Collection<String> graphReactions = []) {
        Map<String, Set<String>> out = new TreeMap<>()
        def add = { String rid, String source -> if (reactionIds.contains(rid)) out.computeIfAbsent(rid) { new TreeSet<String>() } << source }
        table.values().each { if (it.ecs.contains(ec)) add(it.reactionId, 'reactions_tsv') }
        annotationEcs.each { String rid, Set<String> ecs -> if (ecs.contains(ec)) add(rid, 'sbml_annotation') }
        (graphReactions ?: []).each { add(it, 'reaction_graph_id') }
        out
    }

    /**
     * {@code namespace:id} of a compound identifier URI ({@code identifiers.org/chebi/CHEBI:15422},
     * {@code urn:miriam:chebi:CHEBI%3A15422}), or null for anything that is not a compound namespace.
     */
    static String compoundXref(String uri) {
        if (!uri) return null
        java.util.regex.Matcher m = null
        for (String re : ['identifiers\\.org/([A-Za-z_.]+)[/:](.+)$', '^urn:miriam:([A-Za-z_.]+):(.+)$']) {
            java.util.regex.Matcher candidate = java.util.regex.Pattern.compile(re).matcher(uri)
            if (candidate.find()) { m = candidate; break }
        }
        if (m == null) return null
        String ns = m.group(1).toLowerCase(Locale.ROOT)
        if (!COMPOUND_NAMESPACES.contains(ns)) return null
        String id = java.net.URLDecoder.decode(m.group(2), 'UTF-8').trim()
        id ? ns + ':' + id : null
    }

    /**
     * Declared exchanges whose metabolite IS one the reaction handles, by exact identity only: the same SBML
     * species id, or a shared compound cross-reference (so a cytosolic and an extracellular species of one
     * compound match when both carry the same ChEBI / KEGG / BiGG / SEED id). No name, formula or id-suffix
     * comparison is made: an exchange of an unannotated metabolite simply does not match. Result: exchange id
     * to the identity that matched ({@code same_species} or the shared xref).
     */
    Map<String, String> exchangesTouching(Collection<String> reactions, Collection<String> exchanges) {
        Map<String, String> out = new TreeMap<>()
        Set<String> species = new TreeSet<>()
        reactions.each { species.addAll(reactionSpecies.getOrDefault(it, [] as Set)) }
        Set<String> xrefs = new TreeSet<>()
        species.each { xrefs.addAll(speciesXrefs.getOrDefault(it, [] as Set)) }
        exchanges.each { String ex ->
            if (!reactionIds.contains(ex) || reactions.contains(ex)) return
            Set<String> exSpecies = reactionSpecies.getOrDefault(ex, [] as Set)
            if (exSpecies.size() != 1) return                // an exchange has exactly one metabolite
            String sp = exSpecies.first()
            if (species.contains(sp)) { out.put(ex, 'same_species:' + sp); return }
            String shared = (speciesXrefs.getOrDefault(sp, [] as Set)).find { xrefs.contains(it) }
            if (shared != null) out.put(ex, 'xref:' + shared)
        }
        out
    }

    boolean hasGene(String label) { geneByLabel.containsKey(label) }
}
