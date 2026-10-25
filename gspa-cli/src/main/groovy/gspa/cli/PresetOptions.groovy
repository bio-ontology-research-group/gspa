package gspa.cli

import picocli.CommandLine.Option
import picocli.CommandLine.ITypeConverter
import picocli.CommandLine.TypeConversionException
import gspa.config.WorkflowPresets

class PresetOptions {
    @Option(names = '--preset', paramLabel = 'NAME', converter = NameConverter,
        description = 'Workflow defaults: prokaryote, bacteria, archaea, eukaryote, mag, metagenome (alias: metagenomics). See gspa presets. Explicit options override defaults.')
    String name

    @Option(names = '--print-config',
        description = 'Print resolved settings and exit without loading input or running tools.')
    boolean printConfig

    static class NameConverter implements ITypeConverter<String> {
        String convert(String value) {
            try {
                WorkflowPresets.canonicalName(value)
            } catch (IllegalArgumentException e) {
                throw new TypeConversionException(e.message)
            }
        }
    }
}
