package gspa.predictor.specialized

import gspa.integration.EvidenceType
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.predictor.AbstractGenomicRegionSidecarPredictor
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor
import gspa.predictor.PredictorOutcome
import gspa.predictor.PredictorRunStatus
import gspa.predictor.ToolOutputContractException
import gspa.model.Protein

/**
 * CRISPR array predictor using MinCED (CRISPR finder).
 *
 * <p>MinCED is a nucleotide tool: it reports arrays as {@code repeat_region} rows with
 * contig coordinates. Arrays are kept as whole genomic regions in
 * {@code genome.metadata["genomic_regions:crispr"]} (exported to GFF3) and projected onto
 * the proteins whose CDS overlap them through the shared
 * {@link AbstractGenomicRegionSidecarPredictor#projectRegions} step, exactly like the
 * geNomad/CheckV/PhiSpy and antiSMASH regions. Contig-keyed results are never handed to
 * the protein assignment step, which would drop them silently.
 *
 * <p>Array membership is an unscored observation ({@code score_available=false}): MinCED
 * reports repeat counts, not a probability, and a CDS overlapping an array is not thereby
 * assigned a function.
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

    /**
     * Parse MinCED GFF into contig-keyed region annotations. Verified against MinCED
     * 0.4.2: with no array the output file is created with zero bytes and exit 0.
     */
    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def outputFile = requireOutput(outputDir, 'minced_results.gff',
            'minced writes the GFF output file on every successful run (zero bytes when no array is found)', true)

        Map<String, List<Annotation>> results = [:].withDefault { [] }
        int lineNo = 0
        int arrayCount = 0

        outputFile.eachLine { line ->
            lineNo++
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t', -1)
            if (fields.length < 9) {
                throw new ToolOutputContractException("${name}: GFF line ${lineNo} has ${fields.length} columns, GFF3 has 9")
            }

            String contigId = fields[0]
            String featureType = fields[2]
            if (featureType != 'repeat_region') return   // -gffFull also lists repeat units

            int start, end
            try {
                start = fields[3] as int
                end = fields[4] as int
            } catch (NumberFormatException ignored) {
                throw new ToolOutputContractException("${name}: GFF line ${lineNo} has non-numeric coordinates for ${contigId}")
            }
            Integer repeats = fields[5] ==~ /\d+/ ? (fields[5] as Integer) : null
            Map<String, String> attrs = fields[8].split(';').findAll { it.contains('=') }
                .collectEntries { def kv = it.split('=', 2); [(kv[0]): kv[1]] }
            arrayCount++
            String arrayId = attrs.ID ?: "CRISPR_${contigId}_${arrayCount}"
            results[contigId] << new Annotation(
                type: AnnotationType.CRISPR,
                value: 'CRISPR array',
                source: name,
                // Unscored observation: MinCED reports a repeat count, not a confidence.
                score: 0d,
                evidence: 'IEA',
                evidenceType: EvidenceType.GENOMIC_CONTEXT,
                contigId: contigId,
                genomicStart: start,
                genomicEnd: end,
                regionType: 'crispr_array',
                metadata: [
                    arrayId: arrayId,
                    repeatCount: repeats,
                    repeatUnitSequence: attrs.rpt_unit_seq,
                    score_available: false,
                    score_calibrated: false,
                    score_model: 'unscored_region_membership/1',
                    coordinateSystem: 'contig_1_based_inclusive',
                ],
            )
        }
        results
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        lastRunStatus = null
        String key = "genomic_regions:${name}".toString()
        genome.metadata.remove(key)
        long t0 = System.currentTimeMillis()
        def tmpDir = File.createTempDir("gspa_${name}_", '')
        try {
            def input = new File(tmpDir, 'genome.fna')
            def export = writeContigFasta(genome, input)
            if (export.written.isEmpty()) {
                recordStatus(PredictorRunStatus.failed(name,
                    'no contig carries nucleotide sequence; a nucleotide tool cannot run on protein-only input'), t0)
                return [:]
            }
            def result = runOrFail(buildCommand(input, tmpDir), tmpDir, t0)
            if (result == null) return [:]
            Map<String, List<Annotation>> regions = parseGuarded(t0) { parseOutput(tmpDir) }

            def projection = AbstractGenomicRegionSidecarPredictor.projectRegions(genome, regions)
            genome.metadata[key] = projection.validRegions
            int n = projection.validRegions.size()
            int overlaps = countAnnotations(projection.byProtein)
            String note = "${n} CRISPR arrays; ${overlaps} CDS/array overlaps; membership is unscored and asserts no function"

            List<String> problems = []
            List<String> uncovered = []
            if (export.skipped || export.unplacedProteins) {
                problems << describeSkipped(export, genome.contigs.size())
                uncovered.addAll(export.uncovered())
            }
            if (projection.invalidRegions) {
                problems << "${projection.invalidRegions.size()} array row(s) rejected for an unknown contig or invalid coordinates " +
                    "(outside the contig or end < start): ${projection.invalidRegions.take(5).join('; ')}"
                uncovered.addAll(projection.invalidRegionContigs)
            }
            if (projection.invalidCds) {
                problems << "${projection.invalidCds.size()} CDS have invalid or missing coordinates; array membership " +
                    "could not be determined for: ${projection.invalidCds.take(10).join(', ')}${projection.invalidCds.size() > 10 ? ', ...' : ''}"
                uncovered.addAll(projection.invalidCds)
            }
            if (problems) {
                recordStatus(PredictorRunStatus.partial(name, n, note + ' | ' + problems.join(' | '), uncovered.unique(), result.exitCode), t0)
            } else {
                recordStatus(PredictorRunStatus.completed(name, n, result.exitCode, note), t0)
            }
            projection.byProtein
        } catch (Exception e) {
            genome.metadata.remove(key)
            if (lastRunStatus == null || lastRunStatus.outcome != PredictorOutcome.FAILED) {
                recordStatus(PredictorRunStatus.failed(name, "${e.class.simpleName}: ${e.message}"), t0)
            }
            throw e
        } finally {
            tmpDir.deleteDir()
        }
    }

    /** MinCED searches nucleotide sequence; protein FASTA is meaningless to it. */
    @Override
    Map<String, List<Annotation>> predictBatch(List<Protein> proteins) {
        throw new UnsupportedOperationException('crispr (MinCED) is a genome-level predictor; call predictGenome(Genome)')
    }
}
