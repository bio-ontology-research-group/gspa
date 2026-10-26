package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.config.GspaConfig
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.integration.workflow.IntegrationWorkflow
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
 * One enzyme-gap case (a glycolytic aldolase missing between two annotated neighbours, with a
 * candidate in the operon and a distant decoy) pushed through every execution route: the shared
 * executor, single-genome CLI, batch CLI and AnnotationPipeline. The routes must agree on
 * posteriors, suggestions and stage records because they are one implementation. These are
 * mechanism checks on a constructed case, not accuracy evidence.
 */
class WorkflowRouteParitySpec extends Specification {
    @TempDir File dir

    private File file(String name, String content) {
        File f = new File(dir, name)
        f.text = content
        f
    }

    private File claims, ec2go, pathway, reactions, layout, operons, gaps

    def setup() {
        def mapper = new ObjectMapper()
        claims = file('claims.jsonl', [
            [protein_id: 'upstream', function_id: 'GO:0004340', calibrated_prob: 0.95d],
            [protein_id: 'downstream', function_id: 'GO:0004618', calibrated_prob: 0.95d],
            [protein_id: 'candidate', function_id: 'GO:0004332', calibrated_prob: 0.15d],
            [protein_id: 'distant', function_id: 'GO:0004332', calibrated_prob: 0.15d]
        ].collect { mapper.writeValueAsString(it + [function_type: 'GO', source: 'diamond', go_aspect: 'MF']) }.join('\n') + '\n')
        ec2go = file('ec2go.tsv', 'EC:2.7.1.1\tGO:0004340\nEC:4.1.2.13\tGO:0004332\nEC:5.3.1.9\tGO:0004618\n')
        pathway = file('pathways.tsv', 'pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\tdepends_on_reaction\n' +
            'P\tTest carbon route\tGO:0006096\tr1\tEC:2.7.1.1\t-\n' +
            'P\tTest carbon route\tGO:0006096\tr2\tEC:4.1.2.13\tr1\n' +
            'P\tTest carbon route\tGO:0006096\tr3\tEC:5.3.1.9\tr2\n')
        reactions = file('reactions.tsv', 'id\tstoichiometry\tec_numbers\n' +
            'r1\t-1:a:0:0;1:b:0:0\t2.7.1.1\n' +
            'r2\t-1:b:0:0;1:c:0:0\t4.1.2.13\n' +
            'r3\t-1:c:0:0;1:d:0:0\t5.3.1.9\n')
        layout = file('layout.tsv', 'protein_id\tcontig\tstart\tend\tstrand\n' +
            'upstream\tchr\t1000\t1800\t+\ncandidate\tchr\t2000\t2800\t+\n' +
            'downstream\tchr\t3000\t3800\t+\ndistant\tchr\t90000\t91000\t+\n')
        operons = file('operons.tsv', 'upstream\tcandidate\tdownstream\n')
        gaps = file('gaps.jsonl', '{"pathway_id":"P","reaction_id":"r2","ec_number":"4.1.2.13","go_term":"GO:0004332"}\n')
    }

    private IntegrationWorkflowConfig baseConfig() {
        new IntegrationWorkflowConfig(lite: true, enabledPriors: 'none', pathwaysFile: pathway, ec2goFile: ec2go,
            reactionGraphFile: reactions, rlcCurrencyPercentile: 100.0d, rlcSuggester: true, darkMatter: true,
            operonsFile: operons, gapsFile: gaps)
    }

    private Genome genome(String id) {
        Genome g = new Genome(id: id)
        Contig contig = new Contig(id: 'chr')
        g.contigs << contig
        contig.genome = g
        [upstream: [1000, 1800], candidate: [2000, 2800], downstream: [3000, 3800], distant: [90000, 91000]].each { pid, span ->
            contig.addProtein(new Protein(id: pid, sequence: 'M',
                sourceFeature: new Feature(seqId: 'chr', type: FeatureType.CDS, start: span[0], end: span[1], strand: Strand.PLUS)))
        }
        g
    }

    private List<EvidenceClaim> readClaims(File f = claims) {
        new ClaimExtractor().readClaimsJsonl(f)
    }

    private static int execute(List<String> args, StringWriter errors = new StringWriter()) {
        new CommandLine(new IntegrateCommand()).setErr(new PrintWriter(errors)).execute(args as String[])
    }

    private static Map<String, Double> posteriorsFromTsv(File tsv) {
        tsv.readLines().drop(1).collectEntries { line ->
            def f = line.split('\t', -1)
            ["${f[0]}|${f[1]}|${f[2]}".toString(), Double.parseDouble(f[4])]
        }
    }

    private static List stageRecords(List<Map> stages) {
        stages.collect { [it.name, it.status] }
    }

    def 'the enzyme-gap case gives identical posteriors, suggestions and stages on every route'() {
        given:
        List<String> cli = ['--lite', '--enable-priors', 'none', '--pathways', pathway.path, '--ec2go', ec2go.path,
            '--reaction-graph', reactions.path, '--rlc-currency-pct', '100', '--rlc-suggester', '--dark-matter']

        // shared executor
        def directConfig = baseConfig()
        directConfig.genomeLayoutFile = layout
        directConfig.outFile = new File(dir, 'direct.tsv')
        directConfig.suggestionsOut = new File(dir, 'direct.suggestions.tsv')
        def workflow = new IntegrationWorkflow(directConfig)
        def direct = workflow.run(readClaims())
        workflow.writeArtifacts(direct)

        // single-genome CLI
        File single = new File(dir, 'single.tsv')
        File singleSuggestions = new File(dir, 'single.suggestions.tsv')

        // batch CLI
        File batchA = new File(dir, 'batchA.tsv')
        File batchB = new File(dir, 'batchB.tsv')
        File manifest = file('manifest.tsv',
            "genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\tgenome_layout_path\n" +
            "a\t${claims}\t${batchA}\t-\t-\t${operons}\t${gaps}\t-\t-\t-\t-\t${layout}\n" +
            "b\t${claims}\t${batchB}\t-\t-\t${operons}\t${gaps}\t-\t-\t-\t-\t${layout}\n")

        // AnnotationPipeline: raw claims supplied directly, genome keeps its identity,
        // loci come from the genome itself (no layout file).
        def pipelineConfig = new GspaConfig()
        pipelineConfig.integration.enabled = true
        pipelineConfig.integration.workflow = baseConfig()
        pipelineConfig.integration.workflow.suggestionsOut = new File(dir, 'pipeline.suggestions.tsv')
        def pipeline = new AnnotationPipeline(pipelineConfig)
        Genome g = genome('pipeline-genome')

        when:
        assert execute(cli + ['--claims', claims.path, '--out', single.path, '--operons', operons.path,
            '--gaps', gaps.path, '--genome-layout', layout.path, '--suggestions-out', singleSuggestions.path]) == 0
        assert execute(cli + ['--batch-manifest', manifest.path, '--batch-threads', '2', '--suggestions-out', 'unused.tsv']) == 0
        def result = pipeline.integrate(g, readClaims())
        Map<String, Double> pipelineScores = g.proteins.collectMany { p ->
            p.annotations.annotations.collect { a -> ["${p.id}|${a.type}|${a.value}".toString(), a.score] }
        }.collectEntries()
        def mapper = new ObjectMapper()

        then: 'posteriors are identical'
        Map<String, Double> expected = direct.integrated.provenance.collectEntries { k, v -> [k, v.finalProbability] }
        expected.size() >= 4
        posteriorsFromTsv(single) == expected
        posteriorsFromTsv(batchA) == expected
        posteriorsFromTsv(batchB) == expected
        pipelineScores == expected

        and: 'the gap candidate is the context-supported suggestion, the distant decoy is not'
        File directSuggestions = directConfig.suggestionsOut
        directSuggestions.text.contains('reaction_local_context')
        directSuggestions.text.contains('dark_matter')
        directSuggestions.text == singleSuggestions.text
        directSuggestions.text == pipelineConfig.integration.workflow.suggestionsOut.text

        and: 'all routes record the same executed/skipped stages'
        def reference = stageRecords(direct.audit.report().stages)
        reference.find { it[0] == 'reaction_context' }[1] == 'executed'
        reference.find { it[0] == 'dark_matter' }[1] == 'executed'
        reference.find { it[0] == 'metabolic_refinement' }[1] == 'skipped'
        stageRecords(mapper.readValue(new File(single.path + '.stages.json'), Map).stages) == reference
        stageRecords(mapper.readValue(new File(batchA.path + '.stages.json'), Map).stages) == reference
        stageRecords(mapper.readValue(new File(batchB.path + '.stages.json'), Map).stages) == reference
        stageRecords(pipeline.lastIntegrationAudit.report().stages) == reference

        and: 'the pipeline route keeps the genome identity and its provenance'
        result.genomeId == 'pipeline-genome'
        pipeline.lastIntegration.state.genome.is(g)
    }

    def 'dark-matter mode applies only dark-matter suggestions and leaves other suggestion sources on the channel'() {
        given:
        def config = new GspaConfig()
        config.integration.enabled = true
        config.integration.darkMatter.enabled = true
        config.integration.darkMatter.mode = 'top-q'
        config.integration.workflow = baseConfig()
        def pipeline = new AnnotationPipeline(config)
        Genome g = genome('modes')

        when:
        def result = pipeline.integrate(g, readClaims())
        def dm = result.integrated.suggestions.findAll { it.source == 'dark_matter' }
        def rlc = result.integrated.suggestions.findAll { it.source == 'reaction_local_context' }
        def applied = g.proteins.collectMany { it.annotations.annotations }.findAll { it.source == 'dark_matter_suggester' }

        then:
        dm
        rlc
        applied.size() == dm.size()
        applied.every { it.metadata.suggester == 'context_gap' }
    }

    def 'a requested stage without its data fails on every route and is recorded as failed'() {
        given:
        File ecClaims = file('ec.jsonl', '{"protein_id":"enzyme","function_type":"EC","function_id":"1.1.1.1","source":"diamond","calibrated_prob":0.7}\n')
        def config = new GspaConfig()
        config.integration.enabled = true
        config.integration.workflow = new IntegrationWorkflowConfig(lite: true, enabledPriors: 'none', enableFbaRefiner: true,
            outFile: new File(dir, 'pipeline-fail.tsv'))
        def pipeline = new AnnotationPipeline(config)
        Genome g = genome('fail')

        when:
        pipeline.integrate(g, readClaims(ecClaims))

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('FBA reference data directory is required')
        stageRecords(pipeline.lastIntegrationAudit.report().stages).find { it[0] == 'metabolic_refinement' }[1] == 'failed'
        stageRecords(pipeline.lastIntegrationAudit.report().stages).find { it[0] == 'agent_completion' }[1] == 'not_run'
        pipeline.lastIntegration == null
    }

    def 'a model stage with no EC predictions is a failure, not a silent success'() {
        given:
        File fbaData = new File(dir, 'fba')
        fbaData.mkdirs()
        def config = baseConfig()
        config.rlcSuggester = false
        config.darkMatter = false
        config.enableFbaRefiner = true
        config.fbaDataDir = fbaData
        config.outFile = new File(dir, 'noec.tsv')

        when:
        new IntegrationWorkflow(config).run(readClaims())    // GO claims only

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('no EC predictions')
    }

    def 'medium and biomass template reach the model of each genome and parallel genomes stay isolated'() {
        given:
        File stub = file('stub_sidecar.py', '''
import argparse, json, os
p = argparse.ArgumentParser()
for a in ['--fba-data-dir', '--input-json', '--output-json', '--export-model', '--gram', '--medium', '--tau-core', '--alpha-bio']:
    p.add_argument(a)
ns = p.parse_args()
json.dump({'n_core_rxns': 7, 'n_repaired_rxns': 0, 'n_uncertified_additions': 0},
          open(ns.output_json, 'w'))
open(os.path.join(os.path.dirname(ns.output_json), 'args.txt'), 'w').write(ns.gram + ' ' + ns.medium)
''')
        File fbaData = new File(dir, 'fba')
        fbaData.mkdirs()
        def mapper = new ObjectMapper()
        File ecA = file('ecA.jsonl', mapper.writeValueAsString([protein_id: 'a1', function_type: 'EC', function_id: '2.7.1.1', source: 'diamond', calibrated_prob: 0.9d]) + '\n')
        File ecB = file('ecB.jsonl', mapper.writeValueAsString([protein_id: 'b1', function_type: 'EC', function_id: '4.1.2.13', source: 'diamond', calibrated_prob: 0.8d]) + '\n')
        File outA = new File(dir, 'A.tsv')
        File outB = new File(dir, 'B.tsv')
        File manifest = file('fba-manifest.tsv',
            "genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\tgenome_layout_path\torthogroups_path\tcluster_consensus_path\tfba_medium\tfba_gram\n" +
            "A\t${ecA}\t${outA}\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tLB\tpos\n" +
            "B\t${ecB}\t${outB}\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tM9-glucose\tneg\n")
        File single = new File(dir, 'S.tsv')

        when:
        assert execute(['--batch-manifest', manifest.path, '--batch-threads', '2', '--lite', '--enable-priors', 'none',
            '--enable-fba-refiner', '--fba-data-dir', fbaData.path, '--fba-sidecar', stub.path]) == 0
        assert execute(['--claims', ecA.path, '--out', single.path, '--lite', '--enable-priors', 'none',
            '--enable-fba-refiner', '--fba-data-dir', fbaData.path, '--fba-sidecar', stub.path,
            '--fba-medium', 'LB', '--fba-gram', 'pos']) == 0
        def stagesA = new ObjectMapper().readValue(new File(outA.path + '.stages.json'), Map).stages.find { it.name == 'metabolic_refinement' }
        def stagesB = new ObjectMapper().readValue(new File(outB.path + '.stages.json'), Map).stages.find { it.name == 'metabolic_refinement' }

        then: 'each genome ran its own medium and template'
        new File(outA.path + '.metabolism/args.txt').text == 'pos LB'
        new File(outB.path + '.metabolism/args.txt').text == 'neg M9-glucose'
        new File(single.path + '.metabolism/args.txt').text == 'pos LB'

        and: 'the assumptions are part of the stage record'
        stagesA.status == 'executed'
        stagesA.medium == 'LB'
        stagesA.gram == 'pos'
        stagesB.medium == 'M9-glucose'
        stagesA.core_reactions == 7

        and: 'evidence did not cross genomes'
        posteriorsFromTsv(outA).keySet() == ['a1|EC|2.7.1.1'] as Set
        posteriorsFromTsv(outB).keySet() == ['b1|EC|4.1.2.13'] as Set
        posteriorsFromTsv(outA) == posteriorsFromTsv(single)
    }

    def 'a workflow config file reproduces the flag-built run and rejects unknown names'() {
        given:
        File viaFlags = new File(dir, 'flags.tsv')
        File viaFile = new File(dir, 'file.tsv')
        List<String> flags = ['--lite', '--enable-priors', 'none', '--pathways', pathway.path, '--ec2go', ec2go.path,
            '--reaction-graph', reactions.path, '--rlc-currency-pct', '100', '--rlc-suggester', '--dark-matter',
            '--operons', operons.path, '--gaps', gaps.path, '--genome-layout', layout.path]
        File workflowFile = new File(dir, 'workflow.json')
        def command = new IntegrateCommand()
        new CommandLine(command).parseArgs((flags + ['--claims', claims.path, '--out', viaFlags.path]) as String[])
        command.toWorkflowConfig().writeJson(workflowFile)
        def errors = new StringWriter()
        File bad = file('bad.json', '{"enabledPriors": "none", "enablePriors": "x"}')
        File badPrior = file('badprior.json', '{"enabledPriors": "essentiality,not_a_prior"}')

        when:
        assert execute(flags + ['--claims', claims.path, '--out', viaFlags.path]) == 0
        assert execute(['--workflow-config', workflowFile.path, '--claims', claims.path, '--out', viaFile.path]) == 0
        int unknownKey = execute(['--workflow-config', bad.path, '--claims', claims.path, '--out', new File(dir, 'x.tsv').path], errors)
        int unknownPrior = execute(['--workflow-config', badPrior.path, '--claims', claims.path, '--out', new File(dir, 'y.tsv').path], errors)

        then:
        viaFile.text == viaFlags.text
        unknownKey != 0
        unknownPrior != 0
        errors.toString().contains('not_a_prior')
        !new File(dir, 'x.tsv').exists()
        !new File(dir, 'y.tsv').exists()
    }

    def 'annotate maps its integration options onto the shared workflow configuration'() {
        given:
        File workflowFile = new File(dir, 'annotate-workflow.json')
        new IntegrationWorkflowConfig(enabledPriors: 'consistency', enforceStage2: true, fbaMedium: 'LB').writeJson(workflowFile)
        def command = new AnnotateCommand()
        new CommandLine(command).parseArgs('--integration-config', workflowFile.path, '--proteins', claims.path)

        when:
        def resolved = command.resolveConfiguration()
        def workflow = IntegrationWorkflowConfig.fromGspa(resolved)

        then:
        resolved.integration.enabled
        workflow.enabledPriors == 'consistency'
        workflow.enforceStage2
        workflow.fbaMedium == 'LB'

        and: 'integration stays off unless requested'
        def plain = new AnnotateCommand()
        new CommandLine(plain).parseArgs('--proteins', claims.path)
        !plain.resolveConfiguration().integration.enabled
    }

    def 'preset defaults rank below the workflow file, which ranks below explicit options'() {
        given:
        File workflowFile = file('precedence.json', '{"enabledPriors": "consistency"}')
        def command = new IntegrateCommand()
        def line = new CommandLine(command)
        line.parseArgs('--preset', 'bacteria', '--workflow-config', workflowFile.path, '--claims', claims.path,
            '--out', new File(dir, 'p.tsv').path, '--fba-medium', 'LB')
        command.applyPreset()

        when:
        def config = command.toWorkflowConfig()

        then: 'the file wins over the preset for priors; the preset fills what the file leaves open'
        config.enabledPriors == 'consistency'
        config.enforceStage2
        config.essentialProfile == 'bacteria'
        config.lite

        and: 'an explicit option wins over both'
        config.fbaMedium == 'LB'
    }

    def 'a requested reference that does not exist fails before any output is written'() {
        given:
        File out = new File(dir, 'missing.tsv')
        def errors = new StringWriter()

        when:
        int code = execute(['--claims', claims.path, '--out', out.path, '--lite', '--go-owl', new File(dir, 'absent.owl').path], errors)

        then:
        code != 0
        errors.toString().contains('goOwlFile does not exist')
        !out.exists()
    }
}
