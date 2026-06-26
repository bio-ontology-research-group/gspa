package gspa.metrics

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import gspa.ontology.GoOntology
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * Demonstration (not a unit test): runs organism-level taxon-constraint
 * enforcement on a real DeepGO-PlusPlus-Light annotation dump, using the
 * bundled GO taxon constraints + NCBI disjointness backbone. Skipped unless
 * the predictions file is supplied:
 *
 *   ./gradlew :gspa-core:test --tests '*OrganismEnforcementDemoSpec' \
 *     -Dgspa.demo.annotations=/path/m_genitalium_annotations.tsv \
 *     -Dgspa.demo.taxon=NCBITaxon_2 \
 *     -Dgspa.demo.goobo=/path/go-basic.obo
 */
@IgnoreIf({ !System.getProperty('gspa.demo.annotations') })
class OrganismEnforcementDemoSpec extends Specification {

    @TempDir
    Path tmp

    def 'organism-level enforcement on real predictions'() {
        given:
        String annPath = System.getProperty('gspa.demo.annotations')
        String taxon = System.getProperty('gspa.demo.taxon', 'NCBITaxon_2')
        String goobo = System.getProperty('gspa.demo.goobo')

        def tc = new TaxonConstraints()
        tc.loadFromTsv(resource('/taxon-constraints/go-taxon-constraints.tsv'))
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(resource('/taxon-constraints/ncbi-taxon-hierarchy.tsv'))
        checker.organismTaxon = taxon

        GoOntology go = null
        if (goobo) { go = new GoOntology(); go.loadOwl(goobo) }

        def genome = loadGenome(new File(annPath))
        int totalGoBefore = genome.proteins.sum { it.annotations.goAnnotations().size() } as int

        when:
        def enforcer = new ConsistencyEnforcer(goOntology: go, satChecker: checker, mode: 'remove')
        def result = enforcer.enforce(genome)
        int totalGoAfter = genome.proteins.sum { it.annotations.goAnnotations().size() } as int

        then:
        println "================ organism-level enforcement demo ================"
        println "  organism taxon         : ${taxon}"
        println "  propagation (GO obo)   : ${goobo ?: 'none'}"
        println "  proteins               : ${genome.proteinCount}"
        println "  GO annotations before  : ${totalGoBefore}"
        println "  proteins flagged       : ${result.proteinsAffected}"
        println "  annotations removed    : ${result.annotationsAffected}"
        println "  GO annotations after   : ${totalGoAfter}"
        println "  violations found       : ${result.violations}"
        println "================================================================="
        totalGoAfter == totalGoBefore - result.annotationsAffected
    }

    private Genome loadGenome(File tsv) {
        def genome = new Genome(id: 'demo')
        def contig = new Contig(id: 'c1')
        genome.addContig(contig)
        Map<String, Protein> byId = [:]
        boolean header = true
        tsv.eachLine { line ->
            if (header) { header = false; return }
            def f = line.split('\t')
            if (f.length < 4 || f[1] != 'GO') return
            def p = byId.get(f[0])
            if (p == null) {
                p = new Protein(id: f[0], sequence: 'M')
                contig.addProtein(p)
                byId[f[0]] = p
            }
            p.annotations.add(new Annotation(type: AnnotationType.GO, value: f[2],
                score: (f[3] as double), source: 'deepgo-plusplus-light', evidence: 'IEA'))
        }
        genome
    }

    private File resource(String path) {
        def f = tmp.resolve(path.replaceAll('/', '_')).toFile()
        f.withOutputStream { out -> getClass().getResourceAsStream(path).withCloseable { it.transferTo(out) } }
        f
    }
}
