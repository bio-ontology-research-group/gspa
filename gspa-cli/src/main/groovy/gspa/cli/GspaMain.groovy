package gspa.cli

import picocli.CommandLine
import picocli.CommandLine.Command

/**
 * GSPA - Genome-Scale Protein Annotation
 *
 * A comprehensive genome annotation pipeline combining multiple prediction
 * strategies with built-in quality assessment (completeness, coherence, consistency).
 *
 * Subcommands live in sibling files in this package: AnnotateCommand,
 * EvaluateCommand, CompareCommand, ReportCommand, IntegrateCommand.
 */
@Command(
    name = 'gspa',
    description = 'Genome-Scale Protein Annotation pipeline',
    versionProvider = VersionProvider,
    mixinStandardHelpOptions = true,
    subcommands = [
        AnnotateCommand,
        EvaluateCommand,
        CompareCommand,
        ReportCommand,
        IntegrateCommand,
    ]
)
class GspaMain implements Runnable {

    @Override
    void run() {
        CommandLine.usage(this, System.out)
    }

    static void main(String[] args) {
        int exitCode = new CommandLine(new GspaMain()).execute(args)
        System.exit(exitCode)
    }
}
