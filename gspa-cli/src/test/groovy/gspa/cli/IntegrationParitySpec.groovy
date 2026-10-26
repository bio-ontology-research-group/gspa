package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import picocli.CommandLine
import spock.lang.Specification
import spock.lang.TempDir

/** End-to-end route parity with an enzyme gap and a distant alternative locus. */
class IntegrationParitySpec extends Specification {
    @TempDir File dir

    private File file(String name, String content) {
        File f = new File(dir, name)
        f.text = content
        f
    }

    private static void run(List<String> args) {
        def errors = new StringWriter()
        int code = new CommandLine(new IntegrateCommand()).setErr(new PrintWriter(errors)).execute(args as String[])
        assert code == 0: errors.toString()
    }

    def 'batch runs reaction context and dark matter with the same evidence and outputs as a single genome'() {
        given:
        def mapper = new ObjectMapper()
        File claims = file('claims.jsonl', [
            [protein_id: 'upstream', function_id: 'GO:0004340', calibrated_prob: 0.95d],
            [protein_id: 'downstream', function_id: 'GO:0004618', calibrated_prob: 0.95d],
            [protein_id: 'candidate', function_id: 'GO:0004332', calibrated_prob: 0.15d],
            [protein_id: 'distant', function_id: 'GO:0004332', calibrated_prob: 0.15d]
        ].collect { mapper.writeValueAsString(it + [function_type: 'GO', source: 'diamond', go_aspect: 'MF']) }.join('\n') + '\n')
        File ec2go = file('ec2go.tsv', 'EC:2.7.1.1\tGO:0004340\nEC:4.1.2.13\tGO:0004332\nEC:5.3.1.9\tGO:0004618\n')
        File pathway = file('pathways.tsv', 'pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\tdepends_on_reaction\n' +
            'P\tTest carbon route\tGO:0006096\tr1\tEC:2.7.1.1\t-\n' +
            'P\tTest carbon route\tGO:0006096\tr2\tEC:4.1.2.13\tr1\n' +
            'P\tTest carbon route\tGO:0006096\tr3\tEC:5.3.1.9\tr2\n')
        File reactions = file('reactions.tsv', 'id\tstoichiometry\tec_numbers\n' +
            'r1\t-1:a:0:0;1:b:0:0\t2.7.1.1\n' +
            'r2\t-1:b:0:0;1:c:0:0\t4.1.2.13\n' +
            'r3\t-1:c:0:0;1:d:0:0\t5.3.1.9\n')
        File layout = file('layout.tsv', 'protein_id\tcontig\tstart\tend\tstrand\n' +
            'upstream\tchr\t1000\t1800\t+\ncandidate\tchr\t2000\t2800\t+\n' +
            'downstream\tchr\t3000\t3800\t+\ndistant\tchr\t90000\t91000\t+\n')
        File operons = file('operons.tsv', 'upstream\tcandidate\tdownstream\n')
        File gaps = file('gaps.jsonl', '{"pathway_id":"P","reaction_id":"r2","ec_number":"4.1.2.13","go_term":"GO:0004332"}\n')
        File single = new File(dir, 'single.tsv')
        File singleSuggestions = new File(dir, 'single.suggestions.tsv')
        File batch = new File(dir, 'batch.tsv')
        File batch2 = new File(dir, 'batch2.tsv')
        File manifest = file('manifest.tsv',
            "genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\tgenome_layout_path\n" +
            "one\t${claims}\t${batch}\t-\t-\t${operons}\t${gaps}\t-\t-\t-\t-\t${layout}\n" +
            "two\t${claims}\t${batch2}\t-\t-\t${operons}\t${gaps}\t-\t-\t-\t-\t${layout}\n")
        List<String> common = ['--lite', '--enable-priors', 'none', '--pathways', pathway.path,
            '--ec2go', ec2go.path, '--reaction-graph', reactions.path, '--rlc-currency-pct', '100',
            '--rlc-suggester', '--dark-matter', '--suggestions-out', singleSuggestions.path]

        when:
        run(common + ['--claims', claims.path, '--out', single.path, '--operons', operons.path,
            '--gaps', gaps.path, '--genome-layout', layout.path])
        run(common + ['--batch-manifest', manifest.path, '--batch-threads', '2'])

        then:
        single.text == batch.text
        single.text == batch2.text
        singleSuggestions.readLines().size() > 1
        singleSuggestions.text.contains('reaction_local_context')
        singleSuggestions.text == new File(batch.path + '.suggestions.tsv').text
        singleSuggestions.text == new File(batch2.path + '.suggestions.tsv').text
        def stages = mapper.readValue(new File(batch.path + '.stages.json'), Map).stages
        stages.find { it.name == 'reaction_context' }.status == 'executed'
        stages.find { it.name == 'dark_matter' }.status == 'executed'
        stages.find { it.name == 'metabolic_refinement' }.status == 'skipped'
        stages.find { it.name == 'agent_completion' }.status == 'skipped'
    }

    def 'batch rejects duplicate output paths before workers can overwrite results'() {
        given:
        File claims = file('empty.jsonl', '')
        File out = new File(dir, 'same.tsv')
        File manifest = file('duplicates.tsv', "a\t${claims}\t${out}\nb\t${claims}\t${out}\n")
        def errors = new StringWriter()

        when:
        int code = new CommandLine(new IntegrateCommand()).setErr(new PrintWriter(errors)).execute(
            '--batch-manifest', manifest.path, '--enable-priors', 'none', '--lite')

        then:
        code != 0
        errors.toString().contains('Duplicate genome ID or output path')
        !out.exists()
    }

    def 'explicit FBA request with missing model data fails and writes a stage failure'() {
        given:
        File claims = file('claims.jsonl', '{"protein_id":"enzyme","function_type":"EC","function_id":"1.1.1.1","source":"diamond","calibrated_prob":0.7}\n')
        File out = new File(dir, 'failed.tsv')
        def errors = new StringWriter()

        when:
        int code = new CommandLine(new IntegrateCommand()).setErr(new PrintWriter(errors)).execute(
            '--claims', claims.path, '--out', out.path, '--enable-priors', 'none', '--enable-fba-refiner', '--lite')

        then:
        code != 0
        !out.exists()
        def stages = new ObjectMapper().readValue(new File(out.path + '.stages.json'), Map).stages
        stages.find { it.name == 'metabolic_refinement' }.status == 'failed'
        stages.find { it.name == 'agent_completion' }.status == 'not_run'
    }
}
