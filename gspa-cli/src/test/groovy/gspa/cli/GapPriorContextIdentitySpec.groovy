package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import picocli.CommandLine
import spock.lang.Specification

import java.nio.file.Files

class GapPriorContextIdentitySpec extends Specification {
    File dir
    ContextFixtures fx
    final ObjectMapper mapper = new ObjectMapper()

    def setup() {
        File baseTmp = new File(System.getenv('TMPDIR') ?: '/home/leechuck/tmp')
        baseTmp.mkdirs()
        dir = Files.createTempDirectory(baseTmp.toPath(), 'gspa-gap-prior-ctx-').toFile()
        fx = new ContextFixtures(dir)
    }

    def cleanup() {
        dir?.deleteDir()
    }

    private int execute(List<String> args, StringWriter err = new StringWriter()) {
        new CommandLine(new IntegrateCommand()).setErr(new PrintWriter(err)).execute(args as String[])
    }

    private Map normalizedLedger(File outFile) {
        Map ledger = mapper.readValue(new File(outFile.path + '.stages.json'), Map)
        (ledger.stages as List<Map>).each { Map st -> st.remove('elapsed_seconds') }
        ledger
    }

    def 'unchanged code path (no option) and --gap-prior-context-mode none are byte-for-byte identical'() {
        given:
        def files = fx.genome('g1')
        File gaps = new File(dir, 'g1.gaps.jsonl')
        gaps.text = '{"pathway_id":"P","reaction_id":"r2","ec_number":"EC:2.2.2.2","go_term":"GO:0002222"}\n'
        File ec2go = new File(dir, 'ec2go.tsv')
        ec2go.text = 'EC:1.1.1.1\tGO:0001111\nEC:2.2.2.2\tGO:0002222\nEC:3.3.3.3\tGO:0003333\nEC:4.4.4.4\tGO:0004444\n'

        File outUnchanged = new File(dir, 'unchanged.integrated.tsv')
        File outNone = new File(dir, 'none.integrated.tsv')

        List<String> commonArgs = [
            '--claims', files.claims.path,
            '--genome-id', 'g1',
            '--lite',
            '--gaps', gaps.path,
            '--ec2go', ec2go.path,
            '--reaction-graph', fx.reactions.path,
            '--reaction-ec-aliases', fx.aliases.path,
            '--genome-layout', files.layout.path,
        ]

        when:
        int codeUnchanged = execute(commonArgs + ['--out', outUnchanged.path])
        int codeNone = execute(commonArgs + ['--out', outNone.path, '--gap-prior-context-mode', 'none'])

        then:
        codeUnchanged == 0
        codeNone == 0
        outUnchanged.isFile()
        outNone.isFile()
        outUnchanged.bytes == outNone.bytes
        outUnchanged.text == outNone.text
        outUnchanged.text.contains('gap_filling:')
        normalizedLedger(outUnchanged) == normalizedLedger(outNone)
        !new File(outUnchanged.path + '.gap-prior-context.json').exists()
        !new File(outNone.path + '.gap-prior-context.json').exists()
    }

    def 'non-none gap-prior context modes run without starting cohort_context and write .gap-prior-context.json'() {
        given:
        def files = fx.genome('g1')
        // Add GO claims for anchor proteins so r1 and r3 are anchored in GO space
        files.claims << mapper.writeValueAsString([
            protein_id: 'p2', function_type: 'GO', function_id: 'GO:0001111',
            source: 'diamond', evidence_type: 'SEQUENCE_SIMILARITY', raw_score: 0.9d, calibrated_prob: 0.9d
        ]) + '\n'
        File gaps = new File(dir, 'g1.gaps.jsonl')
        gaps.text = '{"pathway_id":"P","reaction_id":"r2","ec_number":"EC:2.2.2.2","go_term":"GO:0002222"}\n'
        File ec2go = new File(dir, 'ec2go.tsv')
        ec2go.text = 'EC:1.1.1.1\tGO:0001111\nEC:2.2.2.2\tGO:0002222\nEC:3.3.3.3\tGO:0003333\nEC:4.4.4.4\tGO:0004444\n'
        File pathways = new File(dir, 'pathway-groups.tsv')
        pathways.text = 'pathway_id\tgroup_id\treaction_id\tsource\nP\tG1\tr2\ttest\nP\tG2\tr1\ttest\n'
        File orgRxns = new File(dir, 'organism-reactions.txt')
        orgRxns.text = '# plausible reactions\nr1\nr3\n'

        when:
        ['universal', 'pathway', 'organism'].each { String mode ->
            File out = new File(dir, "${mode}.integrated.tsv")
            List<String> args = [
                '--claims', files.claims.path,
                '--out', out.path,
                '--genome-id', 'g1',
                '--lite',
                '--gaps', gaps.path,
                '--ec2go', ec2go.path,
                '--reaction-graph', fx.reactions.path,
                '--reaction-ec-aliases', fx.aliases.path,
                '--genome-layout', files.layout.path,
                '--rlc-currency-pct', '100',
                '--gap-prior-context-mode', mode,
                '--gap-prior-context-window-bp', '5000',
            ]
            if (mode == 'pathway') args += ['--gap-prior-pathway-groups', pathways.path]
            if (mode == 'organism') args += ['--gap-prior-organism-reactions', orgRxns.path]
            assert execute(args) == 0
        }

        then:
        ['universal', 'pathway', 'organism'].every { String mode ->
            File out = new File(dir, "${mode}.integrated.tsv")
            File reportFile = new File(out.path + '.gap-prior-context.json')
            assert reportFile.isFile()
            Map report = mapper.readValue(reportFile, Map)
            Map ledger = mapper.readValue(new File(out.path + '.stages.json'), Map)
            Map cohortStage = (ledger.stages as List<Map>).find { it.name == 'cohort_context' }
            Map evStage = (ledger.stages as List<Map>).find { it.name == 'evidence_integration' }
            report.mode == mode &&
                report.window_bp == 5000 &&
                report.gaps_evaluated == 1 &&
                report.gaps_with_context == 1 &&
                cohortStage.status == 'skipped' &&
                evStage.gap_prior_context_mode == mode &&
                evStage.gap_prior_context_window_bp == 5000
        }
    }
}
