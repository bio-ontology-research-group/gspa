package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.model.*
import gspa.predictor.specialized.AntiSmashPredictor
import picocli.CommandLine.Command
import picocli.CommandLine.Option

import java.security.MessageDigest
import java.util.concurrent.Callable

/** Import native regions against the caller's actual protein IDs, without rerunning gene calling. */
@Command(name = 'antismash-claims', mixinStandardHelpOptions = true,
    description = 'Import antiSMASH regions as unscored CDS overlaps; no tool execution or enzyme inference.')
class AntiSmashClaimsCommand implements Callable<Integer> {
    @Option(names = '--results', required = true) File results
    @Option(names = '--genome', required = true) File fasta
    @Option(names = '--gff3', required = true) File gff
    @Option(names = '--proteins', required = true) File proteins
    @Option(names = '--output', required = true) File output

    @Override
    Integer call() {
        File receipt = new File(output.path + '.import.json')
        if (output.exists() || receipt.exists()) throw new IllegalArgumentException('output or import receipt already exists')
        Map<File, String> inputs = [results, fasta, gff, proteins].collectEntries { [(it): sha256(it)] }
        Genome genome = loadCoordinates(fasta, gff, proteins)
        def predictor = new AntiSmashPredictor()
        def imported = predictor.importGenomeResults(genome, results)
        def extractor = new ClaimExtractor()
        List<EvidenceClaim> claims = []
        imported.projection.byProtein.each { id, calls -> calls.each { claims << extractor.fromAnnotation(id, it) } }
        inputs.each { file, hash ->
            if (sha256(file) != hash) throw new IllegalStateException('antiSMASH import input changed during reading')
        }
        output.absoluteFile.parentFile.mkdirs()
        File stagedClaims = File.createTempFile('antismash-claims-', '.jsonl', output.absoluteFile.parentFile)
        File stagedReceipt = File.createTempFile('antismash-receipt-', '.json', output.absoluteFile.parentFile)
        boolean published = false
        try {
            ClaimExtractor.writeClaimsJsonl(claims, stagedClaims)
            stagedReceipt.text = new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString([
                schema: 'gspa.antismash-import/1', action: 'import_existing_output',
                results_sha256: inputs[results], genome_sha256: inputs[fasta],
                gff3_sha256: inputs[gff], proteins_sha256: inputs[proteins], output_sha256: sha256(stagedClaims),
                protein_count: genome.proteinCount, annotation_count: claims.size(),
                region_count: imported.projection.validRegions.collect { it.metadata.region_id }.unique().size(),
                coverage: imported.skipped ? 'partial' : 'complete', unassessed_inputs: imported.skipped,
                records: predictor.recordObservations,
                regions: imported.projection.validRegions.collect { region ->
                    [contig_id: region.contigId, start: region.genomicStart, end: region.genomicEnd,
                     metadata: region.metadata]
                },
                confidence: 'unscored_region_membership',
                limitation: 'Imported predictions, not a new antiSMASH run, verified enzyme functions, expression or biological accuracy.']) + '\n'
            java.nio.file.Files.move(stagedClaims.toPath(), output.toPath())
            published = true
            java.nio.file.Files.move(stagedReceipt.toPath(), receipt.toPath())
        } catch (Exception e) {
            if (published) output.delete()
            throw e
        } finally {
            stagedClaims.delete()
            stagedReceipt.delete()
        }
        0
    }

    /** Explicit GFF IDs, or Prodigal's header ID plus coordinates; never pair proteins by row order. */
    static Genome loadCoordinates(File fasta, File gff, File proteins) {
        gspa.io.GeneModelReader.read(fasta, gff, proteins)
    }

    private static String sha256(File file) {
        MessageDigest digest = MessageDigest.getInstance('SHA-256')
        file.withInputStream { stream ->
            byte[] buffer = new byte[65536]
            int count
            while ((count = stream.read(buffer)) != -1) digest.update(buffer, 0, count)
        }
        digest.digest().encodeHex().toString()
    }
}
