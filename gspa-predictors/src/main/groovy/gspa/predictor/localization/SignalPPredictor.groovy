package gspa.predictor.localization

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * SignalP predictor for signal peptide detection.
 * Identifies secreted proteins and their signal peptide types
 * (Sec/SPI, Sec/SPII for lipoproteins, Tat/SPI).
 *
 * Supports SignalP 5 and 6 output formats.
 */
class SignalPPredictor extends AbstractToolPredictor {

    /** Organism type: gram+, gram-, arch, eukarya */
    String organismType = 'gram-'

    /** Output format: short, long */
    String outputFormat = 'short'

    @Override
    String getName() { 'signalp' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.SIGNAL_PEPTIDE, AnnotationType.SUBCELLULAR_LOCALIZATION] as Set
    }

    @Override
    String getExecutable() { 'signalp' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        [
            executablePath ?: executable,
            '--fasta', inputFasta.absolutePath,
            '--org', organismType,
            '--output_dir', outputDir.absolutePath,
            '--format', outputFormat,
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // SignalP 6 output: prediction_results.txt
        def resultFile = outputDir.listFiles()?.find {
            it.name.contains('prediction_results') || it.name.endsWith('_summary.signalp5')
        }
        if (resultFile == null) {
            // Try SignalP 5 format
            resultFile = outputDir.listFiles()?.find { it.name.endsWith('.signalp5') || it.name.endsWith('.out') }
        }
        if (resultFile == null) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        resultFile.eachLine { line ->
            if (line.startsWith('#') || line.startsWith('ID') || line.trim().isEmpty()) return
            def fields = line.split(/\s+/)
            if (fields.length < 3) return

            String proteinId = fields[0]
            String prediction = fields[1]    // SP(Sec/SPI), LIPO(Sec/SPII), TAT(Tat/SPI), OTHER

            // SignalP6 format: ID Prediction SP_prob LIPO_prob TAT_prob OTHER_prob [CS_pos]
            // Pick the probability column matching the prediction
            double probability
            if (fields.length >= 6) {
                probability = switch (prediction) {
                    case 'SP'    -> safeDouble(fields[2])
                    case 'LIPO'  -> safeDouble(fields[3])
                    case 'TAT'   -> safeDouble(fields[4])
                    case 'OTHER' -> safeDouble(fields[5])
                    default      -> safeDouble(fields[2])
                }
            } else {
                probability = safeDouble(fields[2])
            }

            if (prediction != 'OTHER' && probability > 0.5) {
                results[proteinId] << new Annotation(
                    type: AnnotationType.SIGNAL_PEPTIDE,
                    value: prediction,
                    source: name,
                    score: probability,
                    metadata: [signalType: prediction]
                )

                // Infer subcellular localization
                String localization = switch (prediction) {
                    case 'SP' -> 'extracellular/periplasm'
                    case 'LIPO' -> 'membrane (lipoprotein)'
                    case 'TAT' -> 'extracellular (TAT-secreted)'
                    default -> 'secreted'
                }
                results[proteinId] << new Annotation(
                    type: AnnotationType.SUBCELLULAR_LOCALIZATION,
                    value: localization,
                    source: name,
                    score: probability,
                )
            }
        }
        results
    }

    private static double safeDouble(String s) {
        try { return s as double } catch (Exception e) { return 0.0 }
    }
}
