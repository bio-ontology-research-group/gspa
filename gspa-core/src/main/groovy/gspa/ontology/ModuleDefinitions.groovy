package gspa.ontology

import groovy.json.JsonSlurper
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * KEGG module definitions with their step logic, as written by
 * {@code benchmark/pathway_coherence/fetch_kegg_module_definitions.py}.
 *
 * <p>A module is a boolean expression over KEGG orthologs: required steps (AND), alternative
 * enzymes for a step (OR), multi-subunit complexes (COMPLEX, all components required, optional
 * components marked), steps KEGG has no ortholog for (GAP) and references to other modules.
 * Each ortholog maps to EC numbers, and through ec2go to GO molecular-function terms, which is
 * how the genome's annotation is matched. A step none of whose orthologs reaches a GO term is
 * <i>unmappable</i>: its absence from an annotation says nothing about the genome.</p>
 */
class ModuleDefinitions {

    private static final Logger log = LoggerFactory.getLogger(ModuleDefinitions)

    /** Module id (M00001) to definition. */
    Map<String, ModuleDef> modules = new LinkedHashMap<>()

    /** GO term to the number of modules in which it occurs in a required step (enzyme sharing). */
    Map<String, Integer> goModuleCount = new HashMap<>()

    String source = ''
    String fetchedUtc = ''

    static ModuleDefinitions load(File json, Map<String, String> ec2go) {
        def root = new JsonSlurper().parse(json)
        def defs = new ModuleDefinitions(source: root.source ?: '', fetchedUtc: root.fetched_utc ?: '')
        Map mods = (root.modules ?: [:]) as Map
        for (Map.Entry e : mods.entrySet()) {
            Map m = e.value as Map
            if (m.tree == null) continue
            def md = new ModuleDef(id: m.id as String, name: (m.name ?: '') as String, clazz: (m['class'] ?: '') as String,
                                   definition: (m.definition ?: '') as String, pathways: (m.pathways ?: []) as List<String>,
                                   completeOrganisms: ((m.complete_organisms ?: []) as List<String>) as Set<String>)
            Map orth = (m.orthology ?: [:]) as Map
            for (Map.Entry oe : orth.entrySet()) {
                Map o = oe.value as Map
                List<String> ecs = ((o.ec ?: []) as List).collect { "EC:${it}".toString() }
                Set<String> gos = new LinkedHashSet<>()
                for (String ec : ecs) {
                    // A partial EC (1.1.1.-) maps through ec2go to an enzyme-class term that most genomes carry
                    // somewhere; that is a coarse reference mapping, not evidence for the step, so it is ignored.
                    if (ec.contains('-')) continue
                    String g = ec2go?.get(ec)
                    if (g) gos << g
                }
                md.orthology[oe.key as String] = new KoInfo(ko: oe.key as String, ec: ecs, go: gos, desc: (o.desc ?: '') as String)
            }
            md.tree = Node.from(m.tree as Map)
            md.steps = md.tree instanceof AndNode ? ((AndNode) md.tree).args : [md.tree]
            defs.modules[md.id] = md
        }
        // enzyme sharing: in how many modules does each GO term appear in a required step
        for (ModuleDef md : defs.modules.values()) {
            Set<String> seen = new HashSet<>()
            for (Node step : md.steps) {
                for (String g : step.goTerms(md)) {
                    if (seen.add(g)) defs.goModuleCount.merge(g, 1, Integer::sum)
                }
            }
        }
        log.info("Loaded ${defs.modules.size()} module definitions from ${json.name}")
        defs
    }

    /** Specificity of a GO term: 1 / number of modules that require it (1.0 = unique to one module). */
    double specificity(String go) {
        int n = goModuleCount.getOrDefault(go, 1)
        1.0d / Math.max(1, n)
    }

    static class ModuleDef {
        String id
        String name
        String clazz
        String definition
        List<String> pathways = []
        Set<String> completeOrganisms = [] as Set
        Map<String, KoInfo> orthology = new LinkedHashMap<>()
        Node tree
        /** Top-level required steps (the AND arguments, or the whole tree when it is not an AND). */
        List<Node> steps = []

        boolean isSignature() { clazz?.startsWith('Signature') }

        /** Domain the module name restricts it to ('bacteria', 'archaea', 'eukaryote'), or null. KEGG names
         *  organism-specific variants of one chemistry ("..., eukaryotes", "..., plants", "..., archaea"). */
        String domainHint() {
            String n = (name ?: '').toLowerCase()
            String tail = n.contains(',') ? n.substring(n.lastIndexOf(',') + 1).trim() : ''
            if (tail =~ /^(eukaryotes?|plants?|animals?|fungi|yeast|mammals|vertebrates|insects)\b/) return 'eukaryote'
            if (tail =~ /^archaea\b/) return 'archaea'
            if (tail =~ /^(bacteria|prokaryotes)\b/) return 'bacteria'
            null
        }
    }

    static class KoInfo {
        String ko
        List<String> ec = []
        Set<String> go = [] as Set
        String desc = ''
    }

    /** Boolean step logic. */
    abstract static class Node {
        static Node from(Map m) {
            if (m.containsKey('ko')) return new KoNode(ko: m.ko as String)
            if (m.containsKey('module')) return new ModuleRefNode(module: m.module as String)
            String op = m.op as String
            List<Node> args = ((m.args ?: []) as List).collect { from(it as Map) }
            switch (op) {
                case 'AND': return new AndNode(args: args)
                case 'OR': return new OrNode(args: args)
                case 'COMPLEX': return new ComplexNode(args: args, optional: ((m.optional ?: []) as List).collect { from(it as Map) })
                case 'OPTIONAL': return new OptionalNode(args: args)
                case 'GAP': return new GapNode()
                default: throw new IllegalArgumentException("unknown node op ${op}")
            }
        }

        /** GO terms that can satisfy this node (any alternative, any complex component). */
        abstract Set<String> goTerms(ModuleDef md)

        /** Orthologs under this node. */
        abstract Set<String> kos()

        /** Short human-readable label. */
        abstract String label(ModuleDef md)
    }

    static class KoNode extends Node {
        String ko
        Set<String> goTerms(ModuleDef md) { md.orthology[ko]?.go ?: ([] as Set) }
        Set<String> kos() { [ko] as Set }
        String label(ModuleDef md) { def k = md.orthology[ko]; k?.desc ? "${ko} ${k.desc}".toString() : ko }
    }

    static class ModuleRefNode extends Node {
        String module
        Set<String> goTerms(ModuleDef md) { [] as Set }
        Set<String> kos() { [] as Set }
        String label(ModuleDef md) { "module ${module}".toString() }
    }

    static class GapNode extends Node {
        Set<String> goTerms(ModuleDef md) { [] as Set }
        Set<String> kos() { [] as Set }
        String label(ModuleDef md) { 'no ortholog assigned by KEGG' }
    }

    abstract static class ListNode extends Node {
        List<Node> args = []
        Set<String> goTerms(ModuleDef md) { Set<String> s = new LinkedHashSet<>(); args.each { s.addAll(it.goTerms(md)) }; s }
        Set<String> kos() { Set<String> s = new LinkedHashSet<>(); args.each { s.addAll(it.kos()) }; s }
    }

    static class AndNode extends ListNode {
        String label(ModuleDef md) { args.collect { it.label(md) }.join(' ; ') }
    }

    static class OrNode extends ListNode {
        String label(ModuleDef md) { args.collect { it.label(md) }.join(' | ') }
    }

    static class ComplexNode extends ListNode {
        List<Node> optional = []
        String label(ModuleDef md) { args.collect { it.label(md) }.join(' + ') + (optional ? ' (+ optional ' + optional.collect { it.label(md) }.join(', ') + ')' : '') }
    }

    static class OptionalNode extends ListNode {
        String label(ModuleDef md) { 'optional: ' + args.collect { it.label(md) }.join(' ') }
    }
}
