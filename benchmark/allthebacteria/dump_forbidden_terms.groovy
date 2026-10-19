// Dump the GO terms the GSPA integrator (gspa-unified.jar) treats as forbidden for a given
// organism taxon, using exactly the bundled constraint tables the jar uses in
// `integrate --batch-manifest` mode (taxon-constraints/go-taxon-constraints.tsv,
// taxon-constraints/ncbi-taxon-hierarchy.tsv, anatomical-constraints/*).
//
// Usage: java -cp gspa-unified.jar groovy.ui.GroovyMain dump_forbidden_terms.groovy OUT.tsv TAXON [TAXON ...]
//   TAXON as in the batch manifest (NCBITaxon:2, NCBITaxon:562, ...).
// Output: OUT.tsv            taxon  in_hierarchy  taxon_set_id  anatomical_set_id
//         OUT.sets.tsv       set_id  kind  go_id      (identical sets are shared)
//   kind = taxon       : SatConsistencyChecker.getForbiddenConstrainedTerms() with plain taxon constraints
//          anatomical  : AnatomicalConstraints.computeForbiddenGoTerms()
// (the enforcer additionally treats any descendant of a listed term as forbidden;
//  the aggregator applies that closure with go-dag.tsv)
import gspa.ontology.TaxonConstraints
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.AnatomicalConstraints

def out = new File(args[0])
def taxa = args.drop(1) as List

File tmpFile(String res) {
    def s = this.class.getResourceAsStream(res) ?: TaxonConstraints.getResourceAsStream(res)
    File f = File.createTempFile('gspa', '.tsv'); f.deleteOnExit()
    f.withOutputStream { o -> s.transferTo(o) }
    f
}
def tcFile = tmpFile('/taxon-constraints/go-taxon-constraints.tsv')
def hierFile = tmpFile('/taxon-constraints/ncbi-taxon-hierarchy.tsv')
def tc = new TaxonConstraints(); tc.loadFromTsv(tcFile)
def checker = new SatConsistencyChecker(tc); checker.loadTaxonomyTsv(hierFile)
def ac = AnatomicalConstraints.loadBundled()

@groovy.transform.Field Map<List, Integer> setIds = [:]
@groovy.transform.Field List<List> sets = []
int sid(String kind, Set terms) {
    def key = [kind, terms.sort()]
    Integer id = setIds[key]
    if (id == null) { id = sets.size(); setIds[key] = id; sets << key }
    id
}
out.withWriter { w ->
    w.writeLine('taxon\tin_hierarchy\ttaxon_set_id\tanatomical_set_id')
    taxa.each { raw ->
        String t = raw.startsWith('NCBITaxon') ? raw.replace(':', '_') : "NCBITaxon_${raw}"
        def c = checker.copyForTaxon(t)
        boolean inH = c.ancestorsWithSelf(t).size() > 1
        int a = sid('taxon', c.getForbiddenConstrainedTerms() as Set)
        int b = sid('anatomical', ac.computeForbiddenGoTerms(t, c) as Set)
        w.writeLine("${raw}\t${inH}\t${a}\t${b}")
    }
}
new File(out.path.replaceAll(/\.tsv$/, '') + '.sets.tsv').withWriter { w ->
    w.writeLine('set_id\tkind\tgo_id')
    sets.eachWithIndex { k, i -> k[1].each { g -> w.writeLine("${i}\t${k[0]}\t${g}") } }
}
println "wrote ${out} (${sets.size()} distinct sets)"
