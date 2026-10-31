package gspa.predictor.genecalling

import gspa.io.GffReader
import gspa.io.FastaReader
import gspa.model.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Pyrodigal gene caller — Python binding for Prodigal.
 * Faster alternative with identical results.
 * Falls back to standard Prodigal if pyrodigal is not available.
 */
class PyrodigalCaller implements GeneCaller {

    private static final Logger log = LoggerFactory.getLogger(PyrodigalCaller)

    boolean metaMode = false
    int translationTable = 11
    String executablePath

    @Override
    String getName() { 'pyrodigal' }

    @Override
    boolean isAvailable() {
        try {
            def proc = [executablePath ?: 'pyrodigal', '-V'].execute()
            proc.waitForOrKill(10000)
            return proc.exitValue() == 0
        } catch (Exception e) {
            return false
        }
    }

    @Override
    void callGenes(Genome genome) {
        def tmpDir = File.createTempDir("gspa_pyrodigal_", '')
        try {
            def inputFasta = new File(tmpDir, 'input.fna')
            def gffOutput = new File(tmpDir, 'output.gff')
            def proteinOutput = new File(tmpDir, 'proteins.faa')

            inputFasta.withWriter { writer ->
                genome.contigs.each { contig ->
                    if (contig.sequence) {
                        writer.writeLine(">${contig.id}")
                        writer.writeLine(contig.sequence)
                    }
                }
            }

            def cmd = [
                executablePath ?: 'pyrodigal',
                '-i', inputFasta.absolutePath,
                '-o', gffOutput.absolutePath,
                '-a', proteinOutput.absolutePath,
                '-f', 'gff',
                '-g', translationTable.toString(),
            ]
            if (metaMode) cmd.addAll(['-p', 'meta'])

            log.info("Running pyrodigal: ${cmd.join(' ')}")
            def pb = new ProcessBuilder(cmd)
            pb.directory(tmpDir)
            pb.redirectErrorStream(true)
            def proc = pb.start()
            def output = proc.inputStream.text
            proc.waitFor(30, java.util.concurrent.TimeUnit.MINUTES)

            if (proc.exitValue() != 0) {
                log.error("Pyrodigal failed (exit ${proc.exitValue()}): ${output}")
                return
            }

            if (!gffOutput.exists() || gffOutput.length() == 0) {
                log.warn("Pyrodigal produced no GFF output")
                return
            }

            def gffGenome = GffReader.readGff3(gffOutput)
            def proteins = FastaReader.readProteins(proteinOutput)

            // Match proteins to features by order per contig (same as ProdigalCaller)
            Map<String, List<Protein>> proteinsByContig = [:].withDefault { [] }
            proteins.each { p ->
                String contigId = p.id.replaceAll(/_\d+$/, '')
                proteinsByContig[contigId] << p
            }

            int counter = 1
            gffGenome.contigs.each { gffContig ->
                def genomeContig = genome.findContig(gffContig.id)
                if (!genomeContig) return

                def contigProteins = proteinsByContig[gffContig.id] ?: []
                int protIdx = 0

                gffContig.features.sort()
                gffContig.features.each { feature ->
                    if (feature.type != FeatureType.CDS) return
                    feature.source = 'pyrodigal'
                    genomeContig.addFeature(feature)

                    Protein protein = protIdx < contigProteins.size() ? contigProteins[protIdx] : null
                    protIdx++

                    if (protein) {
                        String tag = "${genome.id}_${String.format('%05d', counter++)}"
                        protein.id = tag
                        protein.sourceFeature = feature
                        feature.id = tag
                        feature.attributes['locus_tag'] = [tag]
                        genomeContig.addProtein(protein)
                    }
                }
                genomeContig.sortFeatures()
            }

            log.info("Pyrodigal called ${genome.proteinCount} genes")
        } finally {
            tmpDir.deleteDir()
        }
    }
}
