package gspa.predictor

import gspa.config.GspaConfig
import gspa.integration.EvidenceClaim
import gspa.integration.IntegrationRunAudit
import gspa.integration.workflow.IntegrationWorkflow
import gspa.integration.workflow.IntegrationWorkflowConfig
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.model.Protein
import gspa.predictor.genecalling.GeneCaller
import groovy.json.JsonSlurper
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Path

class PipelineEmptyInputSpec extends Specification {
    @TempDir Path tmp

    static class ProteinTool implements Predictor {
        int calls
        boolean available = true
        String getName() { 'protein-tool' }
        String getVersion() { 'fixture' }
        Set<AnnotationType> getOutputTypes() { [] as Set }
        boolean isAvailable() { available }
        List<Annotation> predict(Protein p) { calls++; [] }
        Map<String,List<Annotation>> predictBatch(List<Protein> proteins) { calls++; [:] }
    }

    static class GenomeTool extends ProteinTool implements GenomePredictor {
        String observedSequence
        String getName() { 'genome-tool' }
        Map<String,List<Annotation>> predictGenome(Genome genome) {
            calls++
            observedSequence = genome.contigs.first().sequence
            genome.metadata['genome-tool-scanned'] = true
            [:]
        }
    }

    private File fasta(String name, String sequence) {
        def f = tmp.resolve(name).toFile()
        f.text = ">input\n${sequence}\n"
        f
    }

    private AnnotationPipeline pipeline() {
        def cfg = new GspaConfig()
        cfg.integration.enabled = false
        def result = new AnnotationPipeline(cfg)
        result.parallel = false
        result
    }

    private void staleState(AnnotationPipeline p) {
        p.rawClaims = [new EvidenceClaim(proteinId: 'previous-protein')]
        p.lastIntegration = new IntegrationWorkflow.Result(genomeId: 'previous-genome')
        p.lastIntegrationAudit = new IntegrationRunAudit('previous-genome')
    }

    def 'nucleotide tools still scan a genome without proteins while protein tools explicitly skip'() {
        given:
        def p = pipeline()
        def protein = new ProteinTool()
        def genomic = new GenomeTool()
        p.registry.register(protein)
        p.registry.register(genomic)
        p.config.integration.intragenomeCluster.enabled = true

        when:
        def genome = p.annotate(fasta('input.fna', 'ACGTACGT'))

        then:
        genome.proteinCount == 0
        genomic.calls == 1
        genomic.observedSequence == 'ACGTACGT'
        genome.metadata['genome-tool-scanned'] == true
        protein.calls == 0
        p.rawClaims == []
        p.runStatuses.find { it.predictor == 'protein-tool' }.with {
            outcome.name() == 'SKIPPED' && !outcome.isSuccessful() && detail.contains('no protein')
        }
        p.runStatuses.find { it.predictor == 'genome-tool' }.outcome == PredictorOutcome.SUCCESS_ZERO

        when:
        p.writeOutput(genome, tmp.resolve('output').toFile())
        def manifest = new JsonSlurper().parse(tmp.resolve('output/input_predictor_status.json').toFile())

        then:
        manifest.summary.SKIPPED == 1
        manifest.summary.SUCCESS_ZERO == 1
    }

    def 'missing protein executable remains unavailable rather than being called on empty input'() {
        given:
        def p = pipeline()
        def predictor = new ProteinTool(available: false)
        p.registry.register(predictor)

        when:
        p.annotate(fasta('input.fna', 'ACGT'))

        then:
        predictor.calls == 0
        p.runStatuses*.outcome == [PredictorOutcome.UNAVAILABLE]
    }

    def 'zero gene-caller output does not suppress independent nucleotide tools'() {
        given:
        def p = pipeline()
        def genomic = new GenomeTool()
        p.registry.register(genomic)
        p.geneCaller = Stub(GeneCaller) {
            getName() >> 'zero-caller'
            callGenes(_) >> { Genome g -> }
        }

        when:
        p.annotate(fasta('input.fna', 'ACGT'))

        then:
        genomic.calls == 1
        p.runStatuses.find { it.predictor == 'gene-caller:zero-caller' }.outcome == PredictorOutcome.FAILED
        p.runStatuses.find { it.predictor == 'genome-tool' }.outcome == PredictorOutcome.SUCCESS_ZERO
    }

    def 'a new empty input cannot retain old claims or integration results'() {
        given:
        def p = pipeline()
        staleState(p)

        when:
        p.annotate(fasta('input.fna', 'ACGT'))

        then:
        p.rawClaims == []
        p.lastIntegration == null
        p.lastIntegrationAudit == null
    }

    def 'a rejected input also clears state from the previous genome'() {
        given:
        def p = pipeline()
        staleState(p)
        p.config.input.requireExistingGenes = true

        when:
        p.annotate(fasta('input.fna', 'ACGT'))

        then:
        thrown(IllegalArgumentException)
        p.rawClaims == null
        p.lastIntegration == null
        p.lastIntegrationAudit == null
    }

    def 'protein input still executes protein predictors and resets disabled integration state'() {
        given:
        def p = pipeline()
        staleState(p)
        def protein = new ProteinTool()
        p.registry.register(protein)
        p.config.input.proteinsOnly = true
        def proteins = fasta('input.faa', 'MKVG')
        p.config.input.proteinFasta = proteins.absolutePath

        when:
        def genome = p.annotate(proteins)

        then:
        genome.proteinCount == 1
        protein.calls == 1
        p.runStatuses*.outcome == [PredictorOutcome.SUCCESS_ZERO]
        p.rawClaims == []
        p.lastIntegration == null
        p.lastIntegrationAudit == null
    }

    def 'requested integration records the new empty genome instead of silently returning'() {
        given:
        def p = pipeline()
        staleState(p)
        p.config.integration.enabled = true
        p.config.integration.workflow = new IntegrationWorkflowConfig(lite: true, enabledPriors: 'none')
        p.integrationOutputDir = tmp.resolve('integrated').toFile()

        when:
        p.annotate(fasta('input.fna', 'ACGT'))

        then:
        p.lastIntegration.genomeId == 'input'
        p.lastIntegrationAudit.genomeId == 'input'
        p.lastIntegrationAudit.stages.evidence_integration.input_candidates == 0
        p.lastIntegrationAudit.stages.evidence_integration.output_candidates == 0
        p.rawClaims == []
        tmp.resolve('integrated/input_integrated.tsv.stages.json').toFile().isFile()
    }
}
