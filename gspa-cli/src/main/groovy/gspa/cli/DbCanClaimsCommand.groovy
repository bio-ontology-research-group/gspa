package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.io.FastaReader
import gspa.predictor.specialized.DbCanPredictor
import picocli.CommandLine.Command
import picocli.CommandLine.Option

import java.security.MessageDigest
import java.util.concurrent.Callable

/** Import native dbCAN observations through the same parser used by annotate. */
@Command(name = 'dbcan-claims', mixinStandardHelpOptions = true,
    description = 'Import dbCAN overview families as unscored claims; no inference or tool execution.')
class DbCanClaimsCommand implements Callable<Integer> {
    @Option(names = '--overview', required = true) File overview
    @Option(names = '--proteins', required = true) File proteins
    @Option(names = '--output', required = true) File output

    @Override
    Integer call() {
        if (!proteins.isFile()) throw new IllegalArgumentException('protein FASTA is missing')
        if (output.exists() || new File(output.path + '.import.json').exists()) {
            throw new IllegalArgumentException('output or import receipt already exists')
        }
        String overviewHash = sha256(overview), proteinHash = sha256(proteins)
        Set<String> declared = [] as Set
        proteins.withInputStream { stream ->
            InputStream decoded = proteins.name.endsWith('.gz') ? new java.util.zip.GZIPInputStream(stream) : stream
            new BufferedReader(new InputStreamReader(decoded, 'UTF-8')).eachLine { line ->
                if (line.trim().startsWith('>')) {
                    String id = line.trim().substring(1).trim().split(/\s+/)[0]
                    if (!id || !declared.add(id)) throw new IllegalArgumentException('empty or duplicate protein FASTA identifier')
                }
            }
        }
        Set<String> ids = FastaReader.readProteins(proteins)*.id as Set
        if (declared != ids) throw new IllegalArgumentException('protein FASTA readers disagree on identifiers')
        def annotations = new DbCanPredictor().parseOverview(overview, ids)
        def extractor = new ClaimExtractor()
        List<EvidenceClaim> claims = []
        annotations.each { id, calls -> calls.each { claims << extractor.fromAnnotation(id, it) } }
        if (sha256(overview) != overviewHash || sha256(proteins) != proteinHash) {
            throw new IllegalStateException('dbCAN import inputs changed during reading')
        }
        output.parentFile?.mkdirs()
        File receipt = new File(output.path + '.import.json')
        File stagedClaims = File.createTempFile('dbcan-claims-', '.jsonl', output.absoluteFile.parentFile)
        File stagedReceipt = File.createTempFile('dbcan-receipt-', '.json', output.absoluteFile.parentFile)
        boolean claimsPublished = false
        try {
            ClaimExtractor.writeClaimsJsonl(claims, stagedClaims)
            stagedReceipt.text = new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString([
                schema: 'gspa.dbcan-import/1', action: 'import_existing_output',
                overview_sha256: overviewHash, proteins_sha256: proteinHash,
                output_sha256: sha256(stagedClaims), protein_count: ids.size(), annotation_count: claims.size(),
                confidence: 'unscored_family_observations',
                limitation: 'Imported tool output; this command does not establish successful dbCAN execution or biological accuracy.']) + '\n'
            java.nio.file.Files.move(stagedClaims.toPath(), output.toPath())
            claimsPublished = true
            java.nio.file.Files.move(stagedReceipt.toPath(), receipt.toPath())
        } catch (Exception e) {
            if (claimsPublished) output.delete()
            throw e
        } finally {
            stagedClaims.delete()
            stagedReceipt.delete()
        }
        0
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
