package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.predictor.neural.CleanPredictor
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import java.security.MessageDigest
import java.util.concurrent.Callable

@Command(name='clean-claims', mixinStandardHelpOptions=true,
    description='Import typed CLEAN output with its bound FASTA and receipt; no model execution.')
class CleanClaimsCommand implements Callable<Integer> {
    @Option(names='--tsv', required=true) File tsv
    @Option(names='--receipt', required=true) File receipt
    @Option(names='--proteins', required=true) File proteins
    @Option(names='--score-kind', required=true) String scoreKind
    @Option(names='--min-score', defaultValue='0') double minScore
    @Option(names='--output', required=true) File output

    Integer call() {
        File imported = new File(output.path+'.import.json')
        if (output.exists() || imported.exists()) throw new IllegalArgumentException('output or import receipt exists')
        List<File> inputs = [tsv, receipt, proteins]
        Map<String,String> hashes = inputs.collectEntries { [(it.absolutePath):sha(it)] }
        def annotations = new CleanPredictor(scoreKind:scoreKind,minScore:minScore).parseFiles(tsv,receipt,proteins)
        def extractor = new ClaimExtractor()
        List<EvidenceClaim> claims=[]
        annotations.each { id, calls -> calls.each { claims << extractor.fromAnnotation(id,it) } }
        if (inputs.any { sha(it)!=hashes[it.absolutePath] }) throw new IllegalStateException('CLEAN inputs changed during import')
        output.absoluteFile.parentFile.mkdirs()
        File staged=File.createTempFile('clean-claims-','.jsonl',output.absoluteFile.parentFile)
        File stagedReceipt=File.createTempFile('clean-receipt-','.json',output.absoluteFile.parentFile)
        boolean published=false
        try {
            ClaimExtractor.writeClaimsJsonl(claims,staged)
            stagedReceipt.text=new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString([
                schema:'gspa.clean-import/1', action:'import_existing_output',
                inputs:hashes, output_sha256:sha(staged), annotation_count:claims.size(), score_kind:scoreKind,
                limitation:'Receipt-validated import; no assertion of actual model execution, calibration or biological accuracy.'])+'\n'
            java.nio.file.Files.move(staged.toPath(),output.toPath());published=true
            java.nio.file.Files.move(stagedReceipt.toPath(),imported.toPath())
        } catch (Exception failure) {
            if (published) output.delete()
            throw failure
        } finally {
            staged.delete();stagedReceipt.delete()
        }
        0
    }

    private static String sha(File file) {
        MessageDigest digest=MessageDigest.getInstance('SHA-256')
        file.withInputStream { stream ->
            byte[] buffer=new byte[65536];int count
            while ((count=stream.read(buffer))!=-1) digest.update(buffer,0,count)
        }
        digest.digest().encodeHex().toString()
    }
}
