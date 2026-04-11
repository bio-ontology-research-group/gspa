package gspa.integration

import gspa.config.EssentialFunctions
import gspa.io.FastaReader
import gspa.io.GafReader
import gspa.io.GffReader
import gspa.io.GffWriter
import gspa.io.GenbankWriter
import gspa.metrics.QualityReportWriter
import gspa.metrics.HtmlReportWriter
import gspa.model.*
import gspa.ontology.SatConsistencyChecker
import gspa.ontology.TaxonConstraints
import gspa.ontology.PathwayLoader
import gspa.metrics.Completeness
import gspa.metrics.Coherence
import gspa.metrics.Consistency
import spock.lang.Specification
import spock.lang.Shared

/**
 * End-to-end integration tests exercising the full pipeline:
 * load genome -> assign annotations -> evaluate quality -> write output.
 */
class EndToEndSpec extends Specification {

    @Shared
    File tmpDir

    def setupSpec() {
        tmpDir = File.createTempDir("gspa_e2e_", "")
    }

    def cleanupSpec() {
        tmpDir?.deleteDir()
    }

    def "should load genome from GFF3, assign GAF annotations, and evaluate completeness"() {
        given: "A genome with features and GO annotations"
        def gffFile = new File(getClass().getResource('/test-genomes/small_genome.gff3').toURI())
        def gafFile = new File(getClass().getResource('/test-genomes/test_annotations.gaf').toURI())

        when: "Load genome and assign annotations"
        def genome = GffReader.readGff3(gffFile)
        def goAnnotations = GafReader.readGaf(gafFile)

        genome.contigs.each { contig ->
            contig.featuresOfType(FeatureType.CDS).each { feature ->
                def protein = new Protein(id: feature.id, sequence: 'MFAKE', sourceFeature: feature)
                goAnnotations[feature.id]?.each { protein.annotations.add(it) }
                contig.addProtein(protein)
            }
        }

        then: "Genome should have proteins with annotations"
        genome.proteinCount == 4
        genome.allGoTerms().size() > 0
        genome.allGoTerms().containsAll(['GO:0006412', 'GO:0006260', 'GO:0006270'])

        when: "Evaluate completeness"
        def mockGo = new gspa.ontology.GoOntology() {
            Set<String> propagateAnnotations(Set<String> terms) {
                Set<String> result = new HashSet<>(terms)
                if (terms.contains('GO:0006260')) result.add('GO:0006259')
                if (terms.contains('GO:0006270')) { result.add('GO:0006260'); result.add('GO:0006259') }
                return result
            }
        }
        def ef = new EssentialFunctions(profileName: 'test-mini')
        ef.functions['GO:0006259'] = 'Core'  // DNA metabolic process
        ef.functions['GO:0006412'] = 'Core'  // translation

        def completeness = new Completeness(mockGo, ef)
        def result = completeness.evaluate(genome)

        then: "Should have full completeness for these 2 terms"
        result.score == 1.0
        result.presentFunctions.containsAll(['GO:0006259', 'GO:0006412'])
    }

    def "should write annotated genome to all output formats"() {
        given: "A genome with annotations"
        def genome = buildAnnotatedGenome()

        when: "Write GFF3"
        def gffOut = new File(tmpDir, "test_output.gff3")
        GffWriter.writeGff3(genome, gffOut, true, true)

        then: "GFF3 should contain annotations in attributes"
        def gffText = gffOut.text
        gffText.contains('##gff-version 3')
        gffText.contains('Ontology_term=GO:0006412')
        gffText.contains('ec_number=EC:2.7.7.6')
        gffText.contains('Dbxref=Pfam:PF00001')

        when: "Write GenBank"
        def gbkOut = new File(tmpDir, "test_output.gbk")
        GenbankWriter.writeGenbank(genome, gbkOut)

        then: "GenBank should contain annotations as qualifiers"
        def gbkText = gbkOut.text
        gbkText.contains('LOCUS')
        gbkText.contains('/locus_tag="gene_001"')
        gbkText.contains('/EC_number="2.7.7.6"')
        gbkText.contains('GO: GO:0006412')
        gbkText.contains('Pfam: PF00001')
        gbkText.contains('ORIGIN')
    }

    def "should write JSON and TSV quality reports"() {
        given:
        def report = new QualityReport(
            genomeId: 'e2e_test',
            assessmentDate: '2026-04-10',
            completeness: 0.85,
            processCoherence: 0.92,
            pathwayCoherence: 0.78,
            complexCoherence: 1.0,
            consistent: true,
            meanIC: 6.5,
            annotatedProteinCount: 400,
            totalProteinCount: 500,
            presentEssentialFunctions: ['GO:0006412', 'GO:0006260'] as Set,
            missingEssentialFunctions: ['GO:0051301'] as Set,
        )

        when: "Write JSON report"
        def jsonFile = new File(tmpDir, "report.json")
        QualityReportWriter.writeJson(report, jsonFile)

        then:
        jsonFile.exists()
        def json = new groovy.json.JsonSlurper().parse(jsonFile)
        json.genome_id == 'e2e_test'
        json.completeness.score == 0.85
        json.coherence.process_coherence == 0.92
        json.consistency.consistent == true

        when: "Write TSV report"
        def tsvFile = new File(tmpDir, "report.tsv")
        QualityReportWriter.writeSummaryTsv([report], tsvFile)

        then:
        tsvFile.readLines().size() == 2
        tsvFile.readLines()[0].startsWith('genome_id')
        tsvFile.readLines()[1].startsWith('e2e_test')

        when: "Write HTML report"
        def htmlFile = new File(tmpDir, "report.html")
        HtmlReportWriter.writeHtml(report, htmlFile)

        then:
        htmlFile.exists()
        def html = htmlFile.text
        html.contains('e2e_test')
        html.contains('85.0%')
        html.contains('GSPA Quality Report')
        html.contains('GO:0051301')  // missing essential
    }

    def "should handle MAG with adjusted completeness"() {
        given:
        def genome = buildAnnotatedGenome()
        genome.mag = true
        genome.assemblyInfo = new AssemblyInfo(completeness: 60.0, contamination: 5.0)

        def mockGo = new gspa.ontology.GoOntology() {
            Set<String> propagateAnnotations(Set<String> terms) { new HashSet<>(terms) }
        }
        def ef = new EssentialFunctions(profileName: 'test')
        ef.functions['GO:0006412'] = 'Core'
        ef.functions['GO:0006260'] = 'Core'
        ef.functions['GO:0051301'] = 'Core'
        ef.functions['GO:0006351'] = 'Core'

        when:
        def completeness = new Completeness(mockGo, ef)
        def result = completeness.evaluateMAG(genome)

        then: "Raw score should be 1/4 = 0.25 (only GO:0006412 present)"
        result.score == 0.25
        result.adjustedScore != null
        result.magCompleteness == 60.0
        // Adjusted: 0.25 / 0.6 = 0.417
        Math.abs(result.adjustedScore - 0.25/0.6) < 0.01
    }

    def "should detect SAT consistency violation in end-to-end flow"() {
        given: "A genome with conflicting taxon constraints"
        def genome = new Genome(id: 'conflict_test')
        def contig = new Contig(id: 'c1')
        def p1 = new Protein(id: 'p1', sequence: 'MK', sourceFeature: new Feature(seqId: 'c1', type: FeatureType.CDS, start: 1, end: 300, strand: Strand.PLUS))
        def p2 = new Protein(id: 'p2', sequence: 'MA', sourceFeature: new Feature(seqId: 'c1', type: FeatureType.CDS, start: 301, end: 600, strand: Strand.PLUS))
        p1.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0009306', source: 'test'))
        p2.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0015979', source: 'test'))
        contig.addProtein(p1)
        contig.addProtein(p2)
        genome.addContig(contig)

        def constraints = new TaxonConstraints()
        constraints.onlyInTaxon['GO:0009306'] = ['NCBITaxon:2'] as Set
        constraints.neverInTaxon['GO:0015979'] = ['NCBITaxon:2'] as Set

        def checker = new SatConsistencyChecker(constraints)
        checker.loadTaxonomyHierarchy(['NCBITaxon:2': 'NCBITaxon:131567'])

        when:
        def consistency = new Consistency(checker)
        def result = consistency.evaluateWithAttribution(genome)

        then:
        !result.consistent
        result.violations.size() >= 1
        result.violations[0].involvedProteins.containsAll(['p1', 'p2'])
    }

    def "should produce multi-genome HTML comparison report"() {
        given:
        def reports = (1..3).collect { i ->
            new QualityReport(
                genomeId: "genome_${i}",
                assessmentDate: '2026-04-10',
                completeness: 0.5 + i * 0.1,
                processCoherence: 0.8 + i * 0.05,
                pathwayCoherence: 0.7,
                complexCoherence: 1.0,
                consistent: i != 2,
                meanIC: 5.0 + i,
                annotatedProteinCount: 100 * i,
                totalProteinCount: 150 * i,
            )
        }

        when:
        def htmlFile = new File(tmpDir, "comparison.html")
        HtmlReportWriter.writeMultiGenomeHtml(reports, htmlFile)

        then:
        htmlFile.exists()
        def html = htmlFile.text
        html.contains('genome_1')
        html.contains('genome_2')
        html.contains('genome_3')
        html.contains('FAIL')  // genome_2 is inconsistent
        html.contains('PASS')
        html.contains('GSPA Quality Comparison')
    }

    def "should round-trip GFF3 with annotations"() {
        given: "Genome with annotations"
        def genome = buildAnnotatedGenome()
        def gffOut = new File(tmpDir, "roundtrip.gff3")

        when: "Write and read back"
        GffWriter.writeGff3(genome, gffOut, true, true)
        def reloaded = GffReader.readGff3(gffOut)

        then: "Structure should be preserved"
        reloaded.contigs.size() == genome.contigs.size()
        reloaded.features.size() == genome.features.size()
        // Annotations should be in attributes
        def cdsFeature = reloaded.features.find { it.id == 'gene_001' }
        cdsFeature != null
        cdsFeature.attributes['Ontology_term']?.contains('GO:0006412')
        cdsFeature.attributes['ec_number']?.contains('EC:2.7.7.6')
    }

    // --- Helpers ---

    private Genome buildAnnotatedGenome() {
        def genome = new Genome(id: 'test_genome')
        def contig = new Contig(id: 'contig_1', sequence: 'ATGCATGC' * 100)

        def feature = new Feature(
            seqId: 'contig_1', source: 'prodigal', type: FeatureType.CDS,
            start: 1, end: 300, strand: Strand.PLUS
        )
        feature.id = 'gene_001'
        feature.attributes['locus_tag'] = ['gene_001']
        feature.product = 'Test protein'
        contig.addFeature(feature)

        def protein = new Protein(
            id: 'gene_001', sequence: 'MKAILKELLELDL',
            sourceFeature: feature
        )
        protein.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0006412', source: 'diamond', goAspect: 'BP'))
        protein.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0003735', source: 'diamond', goAspect: 'MF'))
        protein.annotations.add(new Annotation(type: AnnotationType.EC, value: 'EC:2.7.7.6', source: 'diamond'))
        protein.annotations.add(new Annotation(type: AnnotationType.PFAM, value: 'PF00001', source: 'interproscan'))
        protein.annotations.add(new Annotation(type: AnnotationType.INTERPRO, value: 'IPR000276', source: 'interproscan'))
        contig.addProtein(protein)

        genome.addContig(contig)
        genome
    }
}
