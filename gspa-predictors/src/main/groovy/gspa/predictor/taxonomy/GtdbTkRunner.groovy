package gspa.predictor.taxonomy

import gspa.model.Genome
import gspa.model.OrganismDomain
import gspa.model.TaxonomyInfo
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Runs GTDB-Tk for taxonomic classification of genomes/MAGs.
 * Parses the classify output to assign GTDB lineage.
 */
class GtdbTkRunner {

    private static final Logger log = LoggerFactory.getLogger(GtdbTkRunner)

    String executablePath
    String gtdbDataPath
    int threads = Runtime.runtime.availableProcessors()

    boolean isAvailable() {
        try {
            def proc = [executablePath ?: 'gtdbtk', '--version'].execute()
            proc.waitForOrKill(10000)
            return proc.exitValue() == 0
        } catch (Exception e) {
            return false
        }
    }

    /**
     * Run GTDB-Tk classify_wf on genome FASTA files.
     */
    Map<String, TaxonomyInfo> classify(File inputDir, File outputDir) {
        outputDir.mkdirs()

        def cmd = [
            executablePath ?: 'gtdbtk',
            'classify_wf',
            '--genome_dir', inputDir.absolutePath,
            '--out_dir', outputDir.absolutePath,
            '--cpus', threads.toString(),
            '-x', 'fna',
        ]
        if (gtdbDataPath) {
            cmd.addAll(['--mash_db', gtdbDataPath])
        }

        log.info("Running GTDB-Tk: ${cmd.join(' ')}")
        def proc = new ProcessBuilder(cmd).start()
        proc.waitFor(240, java.util.concurrent.TimeUnit.MINUTES)

        if (proc.exitValue() != 0) {
            log.error("GTDB-Tk failed: ${proc.errorStream.text}")
            return [:]
        }

        // Parse both bacterial and archaeal summary files
        Map<String, TaxonomyInfo> results = [:]
        def bacFile = new File(outputDir, 'classify/gtdbtk.bac120.summary.tsv')
        def arFile = new File(outputDir, 'classify/gtdbtk.ar53.summary.tsv')
        if (bacFile.exists()) results.putAll(parseSummary(bacFile))
        if (arFile.exists()) results.putAll(parseSummary(arFile))
        results
    }

    /**
     * Parse GTDB-Tk summary TSV output.
     * Format: user_genome\tclassification\t...
     * Classification: d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;...
     */
    static Map<String, TaxonomyInfo> parseSummary(File summaryFile) {
        if (!summaryFile.exists()) return [:]

        Map<String, TaxonomyInfo> results = [:]

        summaryFile.eachLine { line ->
            if (line.startsWith('user_genome')) return  // header
            def fields = line.split('\t')
            if (fields.length < 2) return

            String genomeName = fields[0].trim()
            String classification = fields[1].trim()

            if (classification && classification != 'N/A') {
                results[genomeName] = TaxonomyInfo.fromGtdb(classification)
            }
        }

        log.info("GTDB-Tk: classified ${results.size()} genomes")
        results
    }

    /**
     * Parse and assign GTDB-Tk taxonomy to a genome.
     */
    static void assignTaxonomy(Genome genome, File summaryFile) {
        def results = parseSummary(summaryFile)
        def taxonomy = results[genome.id] ?: results.values().find { true }
        if (taxonomy) {
            genome.taxonomy = taxonomy
            log.info("Assigned taxonomy to ${genome.id}: ${taxonomy.gtdbLineage}")
        }
    }

    /**
     * Extract the NCBI taxonomy hierarchy for taxa mentioned in the lineage.
     * Returns a parent map suitable for the SAT consistency checker.
     *
     * This builds a simplified hierarchy from the GTDB lineage string.
     */
    static Map<String, String> buildTaxonomyHierarchy(TaxonomyInfo taxonomy) {
        if (!taxonomy?.gtdbLineage) return [:]

        Map<String, String> parentMap = [:]
        def ranks = taxonomy.gtdbLineage.split(';')

        // Build parent chain: species -> genus -> family -> ... -> domain
        for (int i = ranks.length - 1; i > 0; i--) {
            String child = gtdbToFakeTaxonId(ranks[i])
            String parent = gtdbToFakeTaxonId(ranks[i - 1])
            if (child && parent) {
                parentMap[child] = parent
            }
        }
        // Top-level domain -> root
        if (ranks.length > 0) {
            parentMap[gtdbToFakeTaxonId(ranks[0])] = 'NCBITaxon:131567'  // cellular organisms
        }

        parentMap
    }

    private static String gtdbToFakeTaxonId(String rank) {
        // Convert GTDB rank (e.g., "d__Bacteria") to a pseudo-taxon ID
        // In a real pipeline, this would map to actual NCBI Taxonomy IDs
        if (!rank || rank.length() <= 3) return null
        String name = rank.substring(3)
        if (!name) return null
        "GTDB:${name}"
    }
}
