package gspa.cli

import gspa.config.WorkflowPresets
import picocli.CommandLine.Command
import picocli.CommandLine.Spec
import picocli.CommandLine.Model.CommandSpec

@Command(name = 'presets', mixinStandardHelpOptions = true,
    description = 'List biological workflow presets for annotate and integrate.')
class PresetsCommand implements Runnable {
    @Spec CommandSpec spec

    void run() {
        def out = spec.commandLine().out
        WorkflowPresets.DESCRIPTIONS.each { name, description ->
            out.println("${name.padRight(13)} ${description}")
        }
        out.println('\nInspect settings: gspa annotate --preset NAME --print-config')
        out.println('                  gspa integrate --preset NAME --print-config')
        out.println('Presets do not install predictors or download databases. No preset enables FBA or the agent.')
        out.println('No --preset: existing command defaults apply. Explicit options override preset values.')
        out.println('For MAGs, use one bin per run (or one claims file per batch row); metagenomics is an alias of metagenome.')
    }
}
