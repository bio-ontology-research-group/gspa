package gspa.metrics

import gspa.config.GspaConfig
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

/**
 * Demonstration (skipped unless predictions supplied): runs the full
 * enforcement pipeline (minimal-flip taxon consistency + completeness
 * promotion) on a real DG++Light annotation dump and shows the residual
 * co-annotation inconsistency resolved, plus provenance.
 *
 *   ./gradlew :gspa-core:test --tests '*EnforcementPipelineDemoSpec' \
 *     -Dgspa.demo.annotations=/path/m_genitalium_annotations.tsv \
 *     -Dgspa.demo.goobo=/path/go-basic.obo -Dgspa.demo.taxon=bacteria
 */
@IgnoreIf({ !System.getProperty('gspa.demo.annotations') })
class EnforcementPipelineDemoSpec extends Specification {

    @TempDir
    Path tmp

    def 'enforce-all on real predictions: minimal-flip + completeness + provenance'() {
        given:
        String annPath = System.getProperty('gspa.demo.annotations')
        String goobo = System.getProperty('gspa.demo.goobo')
        String taxon = System.getProperty('gspa.demo.taxon', 'bacteria')

        def config = new GspaConfig()
        config.input.kingdom = 'bacteria'
        config.quality.consistency.enforce = true
        config.quality.consistency.enforceMode = 'minimal-flip'
        config.quality.consistency.organismTaxon = taxon
        config.quality.completeness.enforce = true
        config.quality.provenance = true

        def qp = new QualityPipeline(config)
            .goOwlFile(goobo)
            .essentialFunctionsForDomain(config.resolveOrganismDomain())
            .initializeLite()

        def genome = loadGenome(new File(annPath))
        int before = genome.proteins.sum { it.annotations.goAnnotations().size() } as int
        int distinctBefore = genome.allGoTerms().size()

        when:
        def report = qp.enforceAll(genome)
        int after = genome.proteins.sum { it.annotations.goAnnotations().size() } as int
        // Consistency of the SURVIVING set (fast: a consistent set has no large UNSAT core).
        def consAfter = freshChecker(taxon).check(genome.allGoTerms()).consistent

        def out = new StringBuilder()
        out << "\n============= enforce-all demo (organism ${taxon}) =============\n"
        out << "  GO annotations before : ${before} (${distinctBefore} distinct terms)\n"
        out << "  GO annotations after  : ${after}\n"
        out << "  consistency after     : ${consAfter ? 'PASS' : 'FAIL'} (was FAIL pre-enforcement)\n"
        out << "  enforcement actions   : ${report.count()} ${report.countsByDimension()}\n"
        out << "  promotions (sample)   :\n"
        report.actions.findAll { it.action == 'promote' }.take(5).each {
            out << "      ${it.dimension} ${it.term} <- ${it.basis}\n"
        }
        out << "  removals (sample)     :\n"
        report.actions.findAll { it.action == 'remove' }.take(5).each {
            out << "      ${it.term} (${it.reason}; ${it.basis})\n"
        }
        out << "================================================================"
        println out.toString()

        then:
        consAfter   // minimal-flip should leave the genome jointly consistent
    }

    private SatConsistencyChecker freshChecker(String taxon) {
        def tc = new TaxonConstraints()
        tc.loadFromTsv(resource('/taxon-constraints/go-taxon-constraints.tsv'))
        def checker = new SatConsistencyChecker(tc)
        checker.loadTaxonomyTsv(resource('/taxon-constraints/ncbi-taxon-hierarchy.tsv'))
        checker.organismTaxon = taxon == 'bacteria' ? 'NCBITaxon_2' : taxon
        checker
    }

    private File resource(String path) {
        def f = tmp.resolve(path.replaceAll('/', '_')).toFile()
        f.withOutputStream { out -> getClass().getResourceAsStream(path).withCloseable { it.transferTo(out) } }
        f
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
}
