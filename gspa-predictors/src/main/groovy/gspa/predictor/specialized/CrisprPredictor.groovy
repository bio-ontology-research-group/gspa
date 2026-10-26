package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor
import gspa.model.Protein

/**
 * CRISPR-Cas system predictor using minced (CRISPR finder).
 * Identifies CRISPR arrays and associated Cas genes.
 */
class CrisprPredictor extends AbstractToolPredictor implements GenomePredictor {

    /** Minimum number of repeats to call a CRISPR array */
    int minRepeats = 3

    @Override
    String getName() { 'crispr' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.CRISPR] as Set }

    @Override
    String getExecutable() { 'minced' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def outputFile = new File(outputDir, 'minced_results.gff')
        [
            executablePath ?: executable,
            '-minNR', minRepeats.toString(),
            '-gffFull',
            inputFasta.absolutePath,
            outputFile.absolutePath,
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // MinCED opens its output file unconditionally (UNVERIFIED here: whether it is zero bytes
        // or carries a '##gff-version' line when no array is found, hence allowEmpty).
        def outputFile = requireOutput(outputDir, 'minced_results.gff',
            'minced writes the GFF output file on every successful run', true)

        Map<String, List<Annotation>> results = [:].withDefault { [] }
        String currentContig = null
        int arrayCount = 0

        outputFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t')
            if (fields.length < 9) return

            String contigId = fields[0]
            String featureType = fields[2]

            if (featureType == 'repeat_region') {
                arrayCount++
                String arrayId = "CRISPR_${contigId}_${arrayCount}"
                // Assign to the contig as a whole
                results[contigId] << new Annotation(
                    type: AnnotationType.CRISPR,
                    value: "CRISPR array",
                    source: name,
                    score: 0.9,
                    metadata: [
                        start: fields[3] as int,
                        end: fields[4] as int,
                        arrayId: arrayId,
                    ]
                )
            }
        }
        results
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        runGenomeTool(genome, 'genome.fna')
    }

    /** MinCED searches nucleotide sequence; protein FASTA is meaningless to it. */
    @Override
    Map<String, List<Annotation>> predictBatch(List<Protein> proteins) {
        throw new UnsupportedOperationException('crispr (MinCED) is a genome-level predictor; call predictGenome(Genome)')
    }
}
