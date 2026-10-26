package gspa.cli

import gspa.config.GspaConfig
import gspa.integration.CalibrationTable
import gspa.integration.ClaimExtractor
import gspa.integration.workflow.IntegrationWorkflowConfig
import gspa.model.Contig
import gspa.model.Feature
import gspa.model.FeatureType
import gspa.model.Genome
import gspa.model.Protein
import gspa.model.Strand
import gspa.predictor.AnnotationPipeline
import picocli.CommandLine
import spock.lang.Specification
import spock.lang.TempDir

/**
 * The AnnotationPipeline route of the context workflow: the cohort stage and the feature export come
 * from the same typed `integration.workflow` block, and a run on the same claims, layout and peers
 * gives byte-identical feature rows to `gspa integrate`. Synthetic data; mechanism checks only.
 */
class ContextPipelineSpec extends Specification {
    @TempDir File dir
    ContextFixtures fx

    def setup() { fx = new ContextFixtures(dir) }

    private Genome genome(String id) {
        Genome g = new Genome(id: id)
        Contig contig = new Contig(id: 'c1')
        g.contigs << contig
        contig.genome = g
        (1..8).each { int i ->
            contig.addProtein(new Protein(id: "p${i}".toString(), sequence: 'M',
                sourceFeature: new Feature(seqId: 'c1', type: FeatureType.CDS, start: i * 1000, end: i * 1000 + 900, strand: Strand.PLUS)))
        }
        g
    }

    private File refs(Map files) {
        File f = new File(dir, 'refs.tsv')
        f.text = 'genome_id\tclaims_path\tlayout_path\torthogroups_path\n' +
            ['g2', 'g3'].collect { "${it}\t${files[it].claims}\t${files[it].layout}\t${files[it].orthogroups}" }.join('\n') + '\n'
        f
    }

    def 'the pipeline route runs the cohort stage and exports the same feature rows as gspa integrate'() {
        given:
        def files = ['g1', 'g2', 'g3'].collectEntries { [it, fx.genome(it)] }
        File meta = fx.metadata(['g1', 'g2', 'g3'])
        File viaPipeline = new File(dir, 'export-pipeline')
        File viaCli = new File(dir, 'export-cli')
        def config = new GspaConfig()
        config.integration.enabled = true
        config.integration.workflow = new IntegrationWorkflowConfig(lite: true, enabledPriors: 'none',
            reactionGraphFile: fx.reactions, reactionEcAliasesFile: fx.aliases, rlcCurrencyPercentile: 100.0d,
            cohortInference: true, cohortMetadataFile: meta, cohortManifestFile: refs(files),
            orthogroupsFile: files.g1.orthogroups, featureExportDir: viaPipeline, outFile: new File(dir, 'pipeline.tsv'))
        def pipeline = new AnnotationPipeline(config)
        def claims = new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(files.g1.claims)

        when:
        def result = pipeline.integrate(genome('g1'), claims)
        int code = new CommandLine(new IntegrateCommand()).execute((['--claims', files.g1.claims.path, '--out', new File(dir, 'cli.tsv').path,
            '--genome-id', 'g1', '--genome-layout', files.g1.layout.path, '--orthogroups', files.g1.orthogroups.path,
            '--cohort-manifest', refs(files).path, '--feature-export-dir', viaCli.path] + fx.cohortArgs(meta)) as String[])
        def stages = pipeline.lastIntegrationAudit.report().stages.collectEntries { [it.name, it.status] }

        then:
        code == 0
        stages.cohort_context == 'executed'
        stages.feature_export == 'executed'
        stages.learned_integration == 'skipped'
        result.context.cohort.containsKey('p3|EC|EC:2.2.2.2')
        new File(viaPipeline, 'candidate_features.jsonl').text == new File(viaCli, 'candidate_features.jsonl').text
        new File(viaPipeline, 'candidate_features.jsonl').text.contains('"cohort.w.available":1.0')
    }

    def 'a learned stage requested through the pipeline config without a model fails the call and is recorded'() {
        given:
        def files = fx.genome('g1')
        def config = new GspaConfig()
        config.integration.enabled = true
        config.integration.workflow = new IntegrationWorkflowConfig(lite: true, enabledPriors: 'none', learnedStage: true,
            outFile: new File(dir, 'pipeline-fail.tsv'))
        def pipeline = new AnnotationPipeline(config)

        when:
        pipeline.integrate(genome('g1'), new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(files.claims))

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('without a model')
        pipeline.lastIntegrationAudit.report().stages.find { it.name == 'reference_preparation' }.status == 'failed'
        pipeline.lastIntegration == null
    }
}
