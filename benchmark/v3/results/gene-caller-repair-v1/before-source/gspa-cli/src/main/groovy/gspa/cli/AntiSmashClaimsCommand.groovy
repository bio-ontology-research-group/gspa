package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.io.FastaReader
import gspa.io.GffReader
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
        Map<String, String> nucleotideHeaders = headers(fasta)
        Genome genome = FastaReader.readGenome(fasta)
        if (!genome.contigs || (genome.contigs*.id as Set) != nucleotideHeaders.keySet() || genome.contigs.any { !it.sequence }) {
            throw new IllegalArgumentException('empty or inconsistent genome FASTA')
        }
        Map<String, Contig> contigs = genome.contigs.collectEntries { [(it.id): it] }
        Map<String, String> proteinHeaders = headers(proteins)
        List<Protein> proteinList = FastaReader.readProteins(proteins)
        if ((proteinList*.id as Set) != proteinHeaders.keySet() || proteinList.any { !it.sequence }) {
            throw new IllegalArgumentException('empty or inconsistent protein sequence')
        }
        List<Feature> features = []
        gff.eachLine { line ->
            if (line == '##FASTA') throw new IllegalArgumentException('use a separate genome FASTA, not embedded GFF sequences')
            if (!line.trim() || line.startsWith('#')) return
            if (line.split('\t', -1).length != 9) throw new IllegalArgumentException('malformed GFF row')
            Feature feature = GffReader.parseFeatureLine(line)
            if (feature.type != FeatureType.CDS) return
            def contig = contigs[feature.seqId]
            if (!contig || feature.start < 1 || feature.end < feature.start || feature.end > contig.length ||
                    !(feature.strand in [Strand.PLUS, Strand.MINUS]) || !(feature.phase in [0, 1, 2])) {
                throw new IllegalArgumentException('invalid CDS coordinate, strand or phase')
            }
            features << feature
            contig.addFeature(feature)
        }
        Map<String, List<Feature>> aliases = [:].withDefault { [] }
        features.each { feature ->
            ['protein_id', 'ID', 'locus_tag'].collectMany { feature.attributes[it] ?: [] }.unique().each {
                aliases[it] << feature
            }
        }
        Set<Feature> assigned = Collections.newSetFromMap(new IdentityHashMap<Feature, Boolean>())
        proteinList.each { protein ->
            String header = proteinHeaders[protein.id]
            def prodigal = header =~ /^\S+\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#\s+(.*)$/
            List<Feature> candidates = new ArrayList<>(aliases[protein.id])
            if (prodigal.matches()) {
                Map attributes = GffReader.parseAttributes(prodigal[0][4])
                (attributes['ID'] ?: []).each { candidates.addAll(aliases[it]) }
            }
            // A feature can be named by both its protein ID and its native caller ID.
            Set<Feature> unique = Collections.newSetFromMap(new IdentityHashMap<Feature, Boolean>())
            unique.addAll(candidates)
            if (unique.size() != 1) throw new IllegalArgumentException("missing or ambiguous CDS identity for ${protein.id}")
            Feature feature = unique.first()
            if (prodigal.matches() && (feature.start != prodigal[0][1].toInteger() ||
                    feature.end != prodigal[0][2].toInteger() ||
                    feature.strand != (prodigal[0][3] == '1' ? Strand.PLUS : Strand.MINUS) ||
                    protein.id.replaceAll(/_\d+$/, '') != feature.seqId)) {
                throw new IllegalArgumentException("Prodigal header disagrees with GFF for ${protein.id}")
            }
            if (!assigned.add(feature)) throw new IllegalArgumentException('multiple proteins map to one CDS')
            protein.sourceFeature = feature
            contigs[feature.seqId].addProtein(protein)
        }
        if (assigned.size() != features.size()) throw new IllegalArgumentException('GFF CDS lacks a protein sequence')
        genome
    }

    private static Map<String, String> headers(File file) {
        Map<String, String> out = [:]
        file.withInputStream { stream ->
            InputStream decoded = file.name.endsWith('.gz') ? new java.util.zip.GZIPInputStream(stream) : stream
            new BufferedReader(new InputStreamReader(decoded, 'UTF-8')).eachLine { line ->
                if (!line.trim().startsWith('>')) return
                String header = line.trim().substring(1).trim()
                String id = header.split(/\s+/)[0]
                if (!id || out.containsKey(id)) throw new IllegalArgumentException('empty or duplicate FASTA identifier')
                out[id] = header
            }
        }
        out
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
