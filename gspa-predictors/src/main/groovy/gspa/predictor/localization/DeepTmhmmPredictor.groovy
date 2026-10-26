package gspa.predictor.localization

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

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
        // DeepTMHMM writes one 3-line record per input protein (also non-TM ones), so a
        // missing or empty topology file is a failed run, not "no TM proteins".
        if (topoFile == null || !topoFile.isFile()) {
            throw new ToolOutputContractException(
                "deeptmhmm: no predicted_topologies.3line (or *.3line) after a clean exit")
        }
        if (topoFile.length() == 0) {
            throw new ToolOutputContractException("deeptmhmm: ${topoFile.name} is empty; DeepTMHMM writes a record for every protein")
        }

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        def lines = topoFile.readLines()
        for (int i = 0; i + 2 < lines.size(); i += 3) {
            String header = lines[i]
            if (!header.startsWith('>')) continue
            String proteinId = header.substring(1).split(/\s+/)[0]
            // lines[i+1] = sequence
            String topology = lines[i + 2]  // topology string: i=inside, o=outside, M=membrane, S=signal

            def helices = transmembraneHelixRegions(topology)

            if (!helices.isEmpty()) {
                results[proteinId] << new Annotation(
                    type: AnnotationType.TRANSMEMBRANE,
                    value: "${helices.size()} TM helices",
                    source: name,
                    score: 0.9,
                    metadata: [
                        tmHelixCount: helices.size(),
                        topology: topology,
                    ]
                )

                results[proteinId] << new Annotation(
                    type: AnnotationType.SUBCELLULAR_LOCALIZATION,
                    value: 'integral membrane',
                    source: name,
                    score: 0.9,
                )

                // One region annotation per helix (1-based inclusive residue
                // indices read straight off the topology string).
                helices.each { int[] span ->
                    results[proteinId] << new Annotation(
                        type: AnnotationType.TRANSMEMBRANE,
                        value: 'tm_helix',
                        source: name,
                        score: 0.9,
                        evidenceType: EvidenceType.SEQUENCE_REGION_ML,
                        regionStart: span[0],
                        regionEnd: span[1],
                        regionType: 'tm_helix',
                    )
                }
            }
        }
        results
    }

    /**
     * Walk the topology string and collect contiguous runs of {@code 'M'}
     * as 1-based inclusive (start, end) residue spans.
     */
    protected static List<int[]> transmembraneHelixRegions(String topology) {
        List<int[]> out = []
        int start = -1
        for (int i = 0; i < topology.length(); i++) {
            char ch = topology.charAt(i)
            if (ch == (char) 'M') {
                if (start < 0) start = i
            } else if (start >= 0) {
                int[] span = [start + 1, i] as int[]
                out.add(span)
                start = -1
            }
        }
        if (start >= 0) {
            int[] span = [start + 1, topology.length()] as int[]
            out.add(span)
        }
        out
    }
}
