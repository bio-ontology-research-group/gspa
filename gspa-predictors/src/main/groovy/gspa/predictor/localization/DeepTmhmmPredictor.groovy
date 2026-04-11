package gspa.predictor.localization

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * DeepTMHMM predictor for transmembrane domain prediction.
 * Identifies transmembrane helices and topology.
 */
class DeepTmhmmPredictor extends AbstractToolPredictor {

    @Override
    String getName() { 'deeptmhmm' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.TRANSMEMBRANE, AnnotationType.SUBCELLULAR_LOCALIZATION] as Set
    }

    @Override
    String getExecutable() { 'deeptmhmm' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        [
            executablePath ?: executable,
            '--fasta', inputFasta.absolutePath,
            '-o', outputDir.absolutePath,
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // DeepTMHMM outputs predicted_topologies.3line
        def topoFile = new File(outputDir, 'predicted_topologies.3line')
        if (!topoFile.exists()) {
            // Also check for biolib output format
            topoFile = outputDir.listFiles()?.find { it.name.endsWith('.3line') }
        }
        if (topoFile == null) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        def lines = topoFile.readLines()
        for (int i = 0; i + 2 < lines.size(); i += 3) {
            String header = lines[i]
            if (!header.startsWith('>')) continue
            String proteinId = header.substring(1).split(/\s+/)[0]
            // lines[i+1] = sequence
            String topology = lines[i + 2]  // topology string: i=inside, o=outside, M=membrane, S=signal

            int tmCount = countTransmembraneHelices(topology)

            if (tmCount > 0) {
                results[proteinId] << new Annotation(
                    type: AnnotationType.TRANSMEMBRANE,
                    value: "${tmCount} TM helices",
                    source: name,
                    score: 0.9,
                    metadata: [
                        tmHelixCount: tmCount,
                        topology: topology,
                    ]
                )

                results[proteinId] << new Annotation(
                    type: AnnotationType.SUBCELLULAR_LOCALIZATION,
                    value: 'integral membrane',
                    source: name,
                    score: 0.9,
                )
            }
        }
        results
    }

    private static int countTransmembraneHelices(String topology) {
        int count = 0
        boolean inTM = false
        topology.each { ch ->
            if (ch == 'M' && !inTM) {
                count++
                inTM = true
            } else if (ch != 'M') {
                inTM = false
            }
        }
        count
    }
}
