package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor
import gspa.predictor.ToolOutputContractException
import gspa.model.Protein

/**
 * antiSMASH predictor for biosynthetic gene cluster (BGC) detection.
 * Identifies secondary metabolite clusters: NRPS, PKS, terpenes, etc.
 * Operates at the genome level (needs full nucleotide context).
 */
class AntiSmashPredictor extends AbstractToolPredictor implements GenomePredictor {

    int threads = Runtime.runtime.availableProcessors()

    /** Strictness level: strict, relaxed, loose */
    String strictness = 'relaxed'

    @Override
    String getName() { 'antismash' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.SECONDARY_METABOLITE] as Set }

    @Override
    String getExecutable() { 'antismash' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        [
            executablePath ?: executable,
            inputFasta.absolutePath,
            '--output-dir', outputDir.absolutePath,
            '--cpus', threads.toString(),
            '--hmmdetection-strictness', strictness,
            '--minimal',  // skip resource-heavy analyses for speed
            '--genefinding-tool', 'prodigal',
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // antiSMASH outputs a JSON file with region/cluster information
        // antiSMASH always writes <record>.json (zero regions = empty 'areas'/'regions' list).
        def jsonFiles = outputDir.listFiles()?.findAll { it.name.endsWith('.json') && it.name != 'regions.js' }
        if (!jsonFiles) {
            throw new ToolOutputContractException('antismash: no results .json in the output directory after a clean exit')
        }

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        jsonFiles.each { jsonFile ->
            try {
                def data = new groovy.json.JsonSlurper().parse(jsonFile)
                def records = data.records ?: []

                records.each { record ->
                    def regions = record.areas ?: record.regions ?: []
                    regions.each { region ->
                        String bgcType = region.products?.join(', ') ?: region.product ?: 'unknown'
                        def cdses = region.cds_members ?: region.genes ?: []

                        cdses.each { cds ->
                            String geneId = cds.locus_tag ?: cds.gene_id ?: cds.toString()
                            results[geneId] << new Annotation(
                                type: AnnotationType.SECONDARY_METABOLITE,
                                value: bgcType,
                                source: name,
                                score: 0.8,
                                metadata: [regionType: bgcType]
                            )
                        }
                    }
                }
            } catch (Exception e) {
                // an unreadable results file must not look like "no clusters"
                throw new ToolOutputContractException("antismash: cannot parse ${jsonFile.name}: ${e.message}")
            }
        }
        results
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        runGenomeTool(genome, 'genome.fna') { File tmp -> new File(tmp, 'output') }
    }

    /** antiSMASH needs nucleotide records (it re-calls genes itself); protein FASTA is invalid input. */
    @Override
    Map<String, List<Annotation>> predictBatch(List<Protein> proteins) {
        throw new UnsupportedOperationException('antismash is a genome-level predictor; call predictGenome(Genome)')
    }

    /**
     * UNVERIFIED: antiSMASH is run with --genefinding-tool prodigal, so its CDS locus tags need
     * not match the GSPA protein ids; unmatched ids are dropped by the pipeline's id join.
     */
    @Override
    protected String successNote() {
        'region members are keyed by antiSMASH locus tags; unmatched ids are not mapped to GSPA proteins (unverified)'
    }
}
