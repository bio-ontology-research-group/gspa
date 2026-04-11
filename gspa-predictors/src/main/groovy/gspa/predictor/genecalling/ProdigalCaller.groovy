package gspa.predictor.genecalling

import gspa.io.GffReader
import gspa.io.FastaReader
import gspa.model.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Prodigal gene caller for prokaryotic genomes.
 * Standard gene caller for bacteria and archaea.
 * Supports single mode (default) and meta mode (for MAGs/metagenomes).
 */
class ProdigalCaller implements GeneCaller {

    private static final Logger log = LoggerFactory.getLogger(ProdigalCaller)

    /** Use meta mode for metagenomes/MAGs */
    boolean metaMode = false

    /** Translation table (11 = standard bacterial, 4 = Mycoplasma) */
    int translationTable = 11

    /** Path to prodigal executable */
    String executablePath

    @Override
    String getName() { 'prodigal' }

    @Override
    boolean isAvailable() {
        try {
            def proc = [executablePath ?: 'prodigal', '-v'].execute()
            proc.waitForOrKill(10000)
            return true
        } catch (Exception e) {
            return false
        }
    }

    @Override
    void callGenes(Genome genome) {
        def tmpDir = File.createTempDir("gspa_prodigal_", '')
        try {
            def inputFasta = new File(tmpDir, 'input.fna')
            def gffOutput = new File(tmpDir, 'prodigal.gff')
            def proteinOutput = new File(tmpDir, 'proteins.faa')

            // Write genome FASTA
            inputFasta.withWriter { writer ->
                genome.contigs.each { contig ->
                    if (contig.sequence) {
                        writer.writeLine(">${contig.id}")
                        writer.writeLine(contig.sequence)
                    }
                }
            }

            // Run prodigal
            def cmd = [
                executablePath ?: 'prodigal',
                '-i', inputFasta.absolutePath,
                '-o', gffOutput.absolutePath,
                '-a', proteinOutput.absolutePath,
                '-f', 'gff',
                '-g', translationTable.toString(),
            ]
            if (metaMode) {
                cmd.addAll(['-p', 'meta'])
            }

            log.info("Running prodigal: ${cmd.join(' ')}")
            def pb = new ProcessBuilder(cmd)
            pb.directory(tmpDir)
            pb.redirectErrorStream(true)  // merge stderr into stdout
            def proc = pb.start()

            // Consume output to prevent process blocking
            def output = proc.inputStream.text
            proc.waitFor(30, java.util.concurrent.TimeUnit.MINUTES)

            if (proc.exitValue() != 0) {
                log.error("Prodigal failed (exit ${proc.exitValue()}): ${output}")
                return
            }

            if (!gffOutput.exists() || gffOutput.length() == 0) {
                log.warn("Prodigal produced no GFF output")
                return
            }

            // Parse GFF output to get features
            def gffGenome = GffReader.readGff3(gffOutput)
            def proteins = FastaReader.readProteins(proteinOutput)

            // Build protein lookup map
            // Prodigal FASTA IDs are like "contigId_N" (e.g., NC_001416.1_1)
            // Prodigal GFF IDs are like "N_N" (e.g., 1_1)
            // We match by position: iterate features and proteins in order per contig
            Map<String, List<Protein>> proteinsByContig = [:].withDefault { [] }
            proteins.each { p ->
                // Extract contig ID from protein ID: everything before the last _N
                String contigId = p.id.replaceAll(/_\d+$/, '')
                proteinsByContig[contigId] << p
            }

            int proteinCounter = 1
            gffGenome.contigs.each { gffContig ->
                def genomeContig = genome.findContig(gffContig.id)
                if (!genomeContig) return

                // Get proteins for this contig, sorted by the trailing index
                def contigProteins = proteinsByContig[gffContig.id] ?: []
                int protIdx = 0

                gffContig.features.sort()
                gffContig.features.each { feature ->
                    if (feature.type != FeatureType.CDS) return

                    feature.source = 'prodigal'
                    genomeContig.addFeature(feature)

                    // Match protein by position order
                    Protein protein = protIdx < contigProteins.size() ? contigProteins[protIdx] : null
                    protIdx++

                    if (protein) {
                        String locusTag = "${genome.id}_${String.format('%05d', proteinCounter++)}"
                        protein.id = locusTag
                        protein.sourceFeature = feature
                        feature.id = locusTag
                        feature.attributes['locus_tag'] = [locusTag]
                        genomeContig.addProtein(protein)
                    }
                }
                genomeContig.sortFeatures()
            }

            log.info("Prodigal called ${genome.proteinCount} genes")
        } finally {
            tmpDir.deleteDir()
        }
    }
}
