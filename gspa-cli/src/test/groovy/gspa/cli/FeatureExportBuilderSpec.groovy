package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import picocli.CommandLine
import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.TempDir

/**
 * A real export, written by the production route, is validated and joined by
 * benchmark/v3/prepare_training_data.py in a subprocess (the dataset builder, not the trainer). Guards the
 * builder against stale or invented export layouts. Needs only a Python 3 (stdlib); set GSPA_TEST_PYTHON to
 * override the interpreter.
 */
@IgnoreIf({ !FeatureExportBuilderSpec.pythonPath().exists() })
class FeatureExportBuilderSpec extends Specification {
    @TempDir File dir

    static File pythonPath() {
        for (String c : [System.getenv('GSPA_TEST_PYTHON'), '/usr/bin/python3', '/tmp/gspa-v3-testenv/bin/python']) {
            if (c && new File(c).exists()) return new File(c)
        }
        new File('/nonexistent-python')
    }

    def 'the dataset builder accepts rows exported by the production route and refuses a drifted layout'() {
        given:
        def fx = new ContextFixtures(dir)
        def ids = ['g1', 'g2', 'g3']
        def files = ids.collectEntries { [it, fx.genome(it)] }
        File export = new File(dir, 'export')
        int code = new CommandLine(new IntegrateCommand()).execute((['--batch-manifest', fx.manifest(ids, files).path,
            '--feature-export-dir', export.path] + fx.cohortArgs(fx.metadata(ids))) as String[])
        File script = new File(dir, 'check.py')
        script.text = FeatureExportBuilderSpec.getResourceAsStream('/export_builder_check.py').getText('UTF-8')
        File v3 = new File('../benchmark/v3').canonicalFile

        when:
        def proc = new ProcessBuilder(pythonPath().path, '-I', script.path, v3.path, export.path).redirectErrorStream(true).start()
        String output = proc.inputStream.getText('UTF-8')
        int exit = proc.waitFor()

        then:
        code == 0
        exit == 0
        def summary = new ObjectMapper().readValue(output.readLines().last(), Map)
        summary.rows_read > 0
        summary.table_rows == 3
        summary.contract_id == 'gspa.context-features/2'

        when: 'a stale manifest schema is introduced'
        File manifest = new File(export, 'g1/export_manifest.json')
        manifest.text = manifest.text.replace('gspa.context-feature-export/1', 'gspa.context-feature-export/0')
        def bad = new ProcessBuilder(pythonPath().path, '-I', script.path, v3.path, export.path).redirectErrorStream(true).start()
        String badOutput = bad.inputStream.getText('UTF-8')

        then:
        bad.waitFor() != 0
        badOutput.contains('export manifest schema')
    }
}
