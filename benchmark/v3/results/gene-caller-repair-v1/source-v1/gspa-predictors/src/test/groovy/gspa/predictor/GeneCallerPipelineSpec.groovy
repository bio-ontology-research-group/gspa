package gspa.predictor

import gspa.config.GspaConfig
import gspa.model.Genome
import gspa.predictor.genecalling.GeneCaller
import spock.lang.Specification
import spock.lang.TempDir

class GeneCallerPipelineSpec extends Specification {
    @TempDir File dir
    static class ReportingCaller implements GeneCaller, StatusReportingPredictor {
        PredictorRunStatus lastRunStatus
        boolean throwError, stale, available = true
        int probes
        String getName() { 'reported' }
        boolean isAvailable() { probes++; available }
        String unavailableReason() { available ? null : 'fixture unavailable' }
        void callGenes(Genome genome) {
            if (throwError) throw new IllegalStateException('native failure')
            if (!stale) lastRunStatus = PredictorRunStatus.completed('gene-caller:reported', 0, 0)
        }
    }
    static class ConfigurablePipeline extends AnnotationPipeline {
        List<GeneCaller> candidates = []
        ConfigurablePipeline(GspaConfig cfg) { super(cfg) }
        @Override
        protected List<GeneCaller> geneCallerCandidates(boolean meta) { candidates }
    }
    File input() {
        def f = new File(dir, 'sample.fna'); f.text = '>input\nATGAAATAA\n'; f
    }
    AnnotationPipeline pipeline() {
        def cfg = new GspaConfig(); cfg.integration.enabled = false
        new AnnotationPipeline(cfg).tap { parallel = false }
    }

    def 'fresh typed zero status reaches the pipeline and independent nucleotide tool still runs'() {
        given:
        def p = pipeline()
        def tool = new PipelineEmptyInputSpec.GenomeTool()
        p.registry.register(tool)
        p.geneCaller = new ReportingCaller()
        when:
        p.annotate(input())
        then:
        tool.calls == 1
        p.runStatuses.find { it.predictor == 'gene-caller:reported' }.with {
            outcome == PredictorOutcome.SUCCESS_ZERO && exitCode == 0
        }
        p.runStatuses.find { it.predictor == 'genome-tool' }.outcome == PredictorOutcome.SUCCESS_ZERO
    }

    def 'exception and stale caller status never suppress nucleotide tools or fabricate success'() {
        given:
        def p = pipeline()
        def tool = new PipelineEmptyInputSpec.GenomeTool()
        p.registry.register(tool)
        p.geneCaller = new ReportingCaller(throwError: throwsError, stale: !throwsError,
            lastRunStatus: PredictorRunStatus.completed('gene-caller:reported', 4, 0))
        when:
        p.annotate(input())
        then:
        tool.calls == 1
        p.runStatuses.find { it.predictor == 'gene-caller:reported' }.outcome == PredictorOutcome.FAILED
        p.runStatuses.find { it.predictor == 'genome-tool' }.outcome == PredictorOutcome.SUCCESS_ZERO
        where:
        throwsError << [true, false]
    }

    def 'configured caller unavailability is reported only when gene calling is needed'() {
        given:
        def cfg = new GspaConfig(); cfg.integration.enabled = false
        def missing = new ReportingCaller(available:false)
        def p = new ConfigurablePipeline(cfg).tap { candidates = [missing]; parallel = false }
        p.configure()
        p.registry = new PredictorRegistry(); p.genomePredictors.clear()
        when:
        p.annotate(input())
        then:
        p.runStatuses.size() == 1
        p.runStatuses[0].predictor == 'gene-caller'
        p.runStatuses[0].outcome == PredictorOutcome.UNAVAILABLE
        p.runStatuses[0].detail.contains('fixture unavailable')
        when:
        def proteins = new File(dir, 'proteins.faa'); proteins.text = '>p1\nMMMM\n'
        cfg.input.proteinsOnly = true; cfg.input.proteinFasta = proteins.path
        p.configure()
        p.registry = new PredictorRegistry(); p.genomePredictors.clear()
        p.annotate(proteins)
        then:
        p.runStatuses.empty
        missing.probes == 1
    }
}
