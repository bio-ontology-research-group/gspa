package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import picocli.CommandLine
import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.TempDir

/**
 * A real export, written by the command from claims we supply, is read by the Python trainer's own
 * reader in a subprocess. Needs numpy; set GSPA_TEST_PYTHON or provide ~/tmp/gspa-v3-testenv. Rows carry
 * no labels, and the trainer refuses them until a (test-only) truth join supplies them.
 */
@IgnoreIf({ !FeatureExportTrainerSpec.pythonPath().exists() })
class FeatureExportTrainerSpec extends Specification {
    @TempDir File dir

    static File pythonPath() { new File(System.getenv('GSPA_TEST_PYTHON') ?: new File(System.getProperty('user.home'), 'tmp/gspa-v3-testenv/bin/python').path) }

    def 'rows exported from the production route are readable by the trainer, bit-exact, and carry no labels'() {
        given:
        def fx = new ContextFixtures(dir)
        def ids = ['g1', 'g2', 'g3']
        def files = ids.collectEntries { [it, fx.genome(it)] }
        File export = new File(dir, 'export')
        int code = new CommandLine(new IntegrateCommand()).execute((['--batch-manifest', fx.manifest(ids, files).path,
            '--feature-export-dir', export.path] + fx.cohortArgs(fx.metadata(ids))) as String[])
        File script = new File(dir, 'check.py')
        script.text = FeatureExportTrainerSpec.getResourceAsStream('/export_trainer_check.py').getText('UTF-8')

        when:
        def proc = new ProcessBuilder(pythonPath().path, '-I', script.path, new File('../scripts').canonicalPath,
            new File(export, 'g1').path).redirectErrorStream(true).start()
        String output = proc.inputStream.getText('UTF-8')
        int exit = proc.waitFor()

        then:
        code == 0
        exit == 0
        def summary = new ObjectMapper().readValue(output.readLines().last(), Map)
        summary.rows > 0
        summary.cohort_rows > 0
        summary.contract_id == 'gspa.context-features/3'
    }
}
