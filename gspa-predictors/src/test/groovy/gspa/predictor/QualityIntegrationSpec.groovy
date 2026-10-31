package gspa.predictor

import gspa.config.GspaConfig
import gspa.integration.workflow.IntegrationWorkflowConfig
import gspa.model.Genome
import spock.lang.Specification

class QualityIntegrationSpec extends Specification {
    def 'standalone integration cannot reuse a previous genomes AntiFam scan outcome'() {
        given:
        def config = new GspaConfig()
        config.integration.workflow = new IntegrationWorkflowConfig(lite: true, enabledPriors: 'none')
        def pipeline = new AnnotationPipeline(config)
        pipeline.runStatuses << PredictorRunStatus.completed('antifam', 0)
        when:
        def unknown = pipeline.integrate(new Genome(id: 'new_genome'), [])
        def known = pipeline.integrate(new Genome(id: 'scanned_genome'), [], 'SUCCESS_ZERO')
        then:
        unknown.qualityReport.antifam_scan_status == 'UNKNOWN'
        !unknown.qualityReport.antifam_scan_complete
        known.qualityReport.antifam_scan_status == 'SUCCESS_ZERO'
        known.qualityReport.antifam_scan_complete
    }
}
