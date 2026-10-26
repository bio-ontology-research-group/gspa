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
        elements(doc, 'reaction').each { org.w3c.dom.Element rx ->
            String rid = rx.getAttribute('id')
            if (!rid) return
            info.reactionIds << rid
            Set<String> genes = new LinkedHashSet<>()
            elements(rx, 'fbc:geneProductRef').each { org.w3c.dom.Element ref -> genes << ref.getAttribute('fbc:geneProduct') }
            info.reactionGenes.put(rid, genes)
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

    /** Reactions that list {@code ec} AND the gene as a candidate and exist in the model; empty without a table. */
    List<String> reactionsFor(String ec, String gene) {
        table.values().findAll { it.ecs.contains(ec) && it.candidateGenes.contains(gene) && reactionIds.contains(it.reactionId) }
            *.reactionId.sort()
    }

    boolean hasGene(String label) { geneByLabel.containsKey(label) }
}
