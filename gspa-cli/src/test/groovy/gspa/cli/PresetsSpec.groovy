package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.config.ConfigLoader
import gspa.config.WorkflowPresets
import picocli.CommandLine
import spock.lang.Specification
import spock.lang.TempDir

class PresetsSpec extends Specification {
    @TempDir File tempDir

    private Map invoke(Object command, List<String> args) {
        def output = new StringWriter()
        def errors = new StringWriter()
        def cli = new CommandLine(command).setOut(new PrintWriter(output)).setErr(new PrintWriter(errors))
        int code = cli.execute(args as String[])
        [code: code, output: output.toString(), errors: errors.toString()]
    }

    def 'catalog and previews need no files or external tools'() {
        expect:
        invoke(new GspaMain(), ['presets']).output.contains('metagenome')

        when:
        def result = invoke(new AnnotateCommand(), ['--preset', name, '--print-config'])
        then:
        result.code == 0
        ConfigLoader.loadFromString(result.output).input.type == type

        where:
        name          | type
        'prokaryote'  | 'genome'
        'bacteria'    | 'genome'
        'archaea'     | 'genome'
        'eukaryote'   | 'genome'
        'mag'         | 'mag'
        'metagenomics' | 'metagenome'
    }

    def 'integration preset assumptions match workflow and explicit switches win'() {
        when:
        def command = new IntegrateCommand()
        def cli = new CommandLine(command)
        cli.parseArgs('--preset', name)
        command.applyPreset()

        then:
        command.enforceStage2 == stage2
        command.essentialProfile == profile
        command.enabledPriors.split(',').contains('essentiality') == essentiality
        command.lite
        !command.enableAgent
        !command.enableFbaRefiner

        where:
        name          | stage2 | profile     | essentiality
        'prokaryote'   | true   | null        | false
        'bacteria'     | true   | 'bacteria'  | true
        'archaea'      | true   | 'archaea'   | true
        'eukaryote'    | true   | 'eukaryote' | true
        'mag'         | false  | null        | false
        'metagenomics' | false  | null        | false
    }

    def 'negative boolean flags and explicit false override integration defaults in either order'() {
        when:
        def command = new IntegrateCommand()
        def cli = new CommandLine(command)
        cli.parseArgs(args as String[])
        command.applyPreset()

        then:
        !command.enforceStage2
        !command.lite
        command.organismTaxon == '562'
        command.enabledPriors == 'consistency'
        command.essentialProfile == 'archaea'

        where:
        args << [
            ['--preset', 'bacteria', '--no-enforce-stage2', '--no-lite', '--taxon', '562', '--enable-priors', 'consistency', '--essential-profile', 'archaea'],
            ['--enforce-stage2=false', '--lite=false', '--taxon', '562', '--enable-priors', 'consistency', '--essential-profile', 'archaea', '--preset', 'bacteria']
        ]
    }

    def 'annotation precedence includes output ontology and disabling enforcement'() {
        given:
        File config = new File(tempDir, 'user.yaml')
        config.text = 'output: {outputDir: from-yaml}\nquality: {goOwlFile: custom-go.owl}'
        def command = new AnnotateCommand()
        new CommandLine(command).parseArgs('--preset', 'bacteria', '-c', config.path,
            '--no-enforce-consistency', '--taxon', '562', '--db', 'ref.dmnd')

        when:
        def resolved = command.resolveConfiguration()
        then:
        resolved.output.outputDir == 'from-yaml'
        resolved.quality.goOwlFile == 'custom-go.owl'
        !resolved.quality.consistency.enforce
        resolved.quality.consistency.organismTaxon == '562'
        resolved.predictors.similarity.enabled
        resolved.predictors.similarity.database == 'ref.dmnd'
    }

    def 'explicit taxon inference replaces preset taxon'() {
        given:
        def command = new AnnotateCommand()
        new CommandLine(command).parseArgs('--preset', 'bacteria', '--infer-taxon')
        expect:
        command.resolveConfiguration().quality.consistency.organismTaxon == null
    }

    def 'integration taxon inference is not shadowed by the preset domain'() {
        given:
        def command = new IntegrateCommand()
        new CommandLine(command).parseArgs('--preset', 'bacteria', '--infer-taxon')
        command.applyPreset()
        expect:
        command.inferTaxon
        command.organismTaxon == null
    }

    def 'batch rows retain overrides on top of a workflow preset'() {
        given:
        File claims = new File(tempDir, 'batch-claims.jsonl')
        claims.text = '{"protein_id":"p1","function_type":"GO","function_id":"GO:0003824","source":"diamond","calibrated_prob":0.9}\n'
        File out = new File(tempDir, 'batch.tsv')
        File manifest = new File(tempDir, 'manifest.tsv')
        manifest.text = "genome_id\tclaims_path\tout_path\ttaxon\tessential_profile\toperons_path\tgaps_path\tmax_iter\tenforce_stage2\tenforce_mode\tenable_priors\n" +
            "bin1\t${claims.path}\t${out.path}\tarchaea\tarchaea\t-\t-\t1\tfalse\tremove\tconsistency\n"

        when:
        def result = invoke(new IntegrateCommand(), ['--preset', 'mag', '--batch-manifest', manifest.path])
        then:
        result.code == 0
        out.text.contains('p1')
    }

    def 'invalid presets and unsupported raw eukaryotic gene calling fail helpfully'() {
        expect:
        invoke(new AnnotateCommand(), ['--preset', 'typo', '--print-config']).code != 0
        invoke(new IntegrateCommand(), ['--preset', 'typo', '--print-config']).code != 0
        invoke(new AnnotateCommand(), ['--preset', 'eukaryote', '-i', 'raw.fna']).errors.contains('splice-aware')
    }

    def 'integration preview is JSON and does not create output'() {
        given:
        File out = new File(tempDir, 'absent.tsv')
        when:
        def result = invoke(new IntegrateCommand(), ['--preset', 'metagenomics', '--out', out.path, '--print-config'])
        def settings = new ObjectMapper().readValue(result.output, Map)
        then:
        result.code == 0
        settings['--enable-priors'] == ''
        settings['--enforce-stage2'] == false
        !out.exists()
    }

    def 'no preset keeps integration defaults'() {
        given:
        def command = new IntegrateCommand()
        new CommandLine(command).parseArgs()
        command.applyPreset()
        expect:
        !command.lite
        !command.enforceStage2
        command.essentialProfile == null
        command.enabledPriors == 'essentiality,coherence,consistency,anatomical,gap_filling,genomic_context'
    }

    def 'annotation runs with a controlled sidecar and respects YAML output and disabled operons'() {
        given:
        // A tiny deterministic executable exercises the real sidecar protocol without
        // model downloads or external bioinformatics tools. It emits one EC claim.
        File runner = new File(tempDir, 'fake-python')
        runner.text = '''#!/bin/sh
if [ "$1" = "--version" ]; then echo fixture; exit 0; fi
while [ "$#" -gt 0 ]; do
  if [ "$1" = "--manifest" ]; then shift; manifest="$1"; fi
  shift
done
outdir=$(awk -F '\\t' 'NR == 2 {print $3}' "$manifest")
printf 'protein_id\\tterm\\tscore\\tannotation_type\\np1\\t1.1.1.1\\t0.9\\tEC\\n' > "$outdir/query.clean.tsv"
'''
        runner.setExecutable(true)
        File sidecar = new File(tempDir, 'sidecar.py')
        sidecar.text = '# fixture'
        File fasta = new File(tempDir, 'proteins.faa')
        fasta.text = '>p1\nMAAAAACCCCCC\n'
        File out = new File(tempDir, 'annotation-output')
        File yaml = new File(tempDir, 'annotation.yaml')
        yaml.text = """
predictors:
  domains:
    enabled: false
  neural:
    pythonExecutable: '${runner.path}'
    sidecarScript: '${sidecar.path}'
    clean:
      enabled: true
      modelDir: '${tempDir.path}'
output:
  outputDir: '${out.path}'
"""

        when:
        def result = invoke(new AnnotateCommand(), ['--preset', 'eukaryote', '--proteins', fasta.path, '-c', yaml.path])

        then:
        result.code == 0
        out.listFiles().find { it.name.endsWith('_annotations.tsv') }.text.contains('1.1.1.1')
        !new File(out, 'operons.tsv').exists()
    }

    def 'each integration preset executes with bundled resources and produces annotations'() {
        given:
        File claims = new File(tempDir, 'claims.jsonl')
        claims.text = '{"protein_id":"p1","function_type":"GO","function_id":"GO:0003824","source":"diamond","raw_score":0.9,"calibrated_prob":0.9}\n'
        File out = new File(tempDir, "${name}.tsv")

        when:
        def result = invoke(new IntegrateCommand(), ['--preset', name, '--claims', claims.path, '--out', out.path])

        then:
        result.code == 0
        out.text.contains('p1')
        out.text.contains('GO:0003824')

        where:
        name << WorkflowPresets.DESCRIPTIONS.keySet().toList()
    }
}
