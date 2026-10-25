package gspa.cli

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import java.util.concurrent.Callable

/** Ship the same stdlib-only downloader in the standalone jar and source checkout. */
class ResourceRunner {
    static int run(String python, List<String> args) {
        File app = File.createTempFile('gspa-resources-', '.pyz')
        try {
            def stream = ResourceRunner.getResourceAsStream('/downloads/gspa-resources.pyz')
            if (stream == null) throw new IllegalStateException('Resource manager missing; rebuild the GSPA jar.')
            stream.withCloseable { input -> app.withOutputStream { it << input } }
            def process = new ProcessBuilder([python, app.absolutePath] + args).inheritIO().start()
            process.waitFor()
        } catch (IOException e) {
            System.err.println("Resource manager could not run '${python}'. Use --python with a Python 3.9+ interpreter. ${e.message}")
            return 1
        } finally {
            app.delete()
        }
    }
}

@Command(name = 'resources', mixinStandardHelpOptions = true,
    description = 'List downloadable references/models, or show one resource and its setup requirements.')
class ResourcesCommand implements Callable<Integer> {
    @Parameters(index = '0', arity = '0..1', paramLabel = 'ID', description = 'Resource ID for full details')
    String resourceId
    @Option(names = '--json', description = 'Print the entire resource catalog as JSON')
    boolean json
    @Option(names = '--python', defaultValue = 'python3', description = 'Python 3.9+ interpreter (standard library only)')
    String python

    Integer call() {
        List<String> args = ['list']
        if (resourceId) args << resourceId
        if (json) args << '--json'
        ResourceRunner.run(python, args)
    }
}

@Command(name = 'download', mixinStandardHelpOptions = true,
    description = 'Download/export GSPA resources and models. Use gspa resources for IDs; --all includes optional resources.')
class DownloadCommand implements Callable<Integer> {
    @Parameters(arity = '0..*', paramLabel = 'ID', description = 'Resource IDs (dependencies are included)')
    List<String> ids = []
    @Option(names = '--all', description = 'Attempt every resource; report unavailable/licensed/custom resources as failures')
    boolean all
    @Option(names = '--group', description = 'Resource group: core, annotation, models, databases, bundled, deepgo-plusplus (repeatable)')
    List<String> groups = []
    @Option(names = ['--dir', '--dest'], description = 'Destination (default: GSPA_DATA_DIR or ~/.local/share/gspa)')
    File destination
    @Option(names = '--dry-run', description = 'Show downloads, commands and prerequisites without network access or filesystem changes')
    boolean dryRun
    @Option(names = '--force', description = 'Download again and replace managed resources after successful verification')
    boolean force
    @Option(names = '--source-url', description = 'Override a one-file resource URL, or supply a custom resource bundle')
    String sourceUrl
    @Option(names = '--sha256', description = 'Expected SHA-256 of --source-url')
    String sha256
    @Option(names = '--python', defaultValue = 'python3', description = 'Python 3.9+ interpreter (standard library only)')
    String python

    Integer call() {
        List<String> args = ['download'] + ids
        if (all) args << '--all'
        groups.each { args.addAll(['--group', it]) }
        if (destination) args.addAll(['--dir', destination.path])
        if (dryRun) args << '--dry-run'
        if (force) args << '--force'
        if (sourceUrl) args.addAll(['--source-url', sourceUrl])
        if (sha256) args.addAll(['--sha256', sha256])
        ResourceRunner.run(python, args)
    }
}
