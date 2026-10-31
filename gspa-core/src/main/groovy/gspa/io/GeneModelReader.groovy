package gspa.io

import gspa.model.*

/**
 * Strict one-to-one reader for separate nucleotide FASTA, single-row CDS GFF and
 * protein FASTA. Native Prodigal headers may supply additional identity evidence.
 * This is not a spliced/eukaryotic transcript assembler. Builds an isolated model;
 * malformed or incomplete output is never attached to a caller-supplied genome.
 */
class GeneModelReader {
    static Genome read(File fasta, File gff, File proteins) {
        Map<String, String> nucleotideHeaders = readHeaders(fasta)
        Genome genome = FastaReader.readGenome(fasta)
        if (!genome.contigs || (genome.contigs*.id as Set) != nucleotideHeaders.keySet() || genome.contigs.any { !it.sequence }) {
            throw new IllegalArgumentException('empty or inconsistent genome FASTA')
        }
        Map<String, Contig> contigs = genome.contigs.collectEntries { [(it.id): it] }
        Map<String, String> proteinHeaders = readHeaders(proteins)
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
        Map<String, List<Feature>> aliases = indexAliases(features)
        Set<Feature> assigned = Collections.newSetFromMap(new IdentityHashMap<Feature, Boolean>())
        proteinList.each { protein ->
            Set<Feature> unique = matchingFeatures(protein.id, proteinHeaders[protein.id], aliases)
            if (unique.size() != 1) throw new IllegalArgumentException("missing or ambiguous CDS identity for ${protein.id}")
            Feature feature = unique.first()
            if (!assigned.add(feature)) throw new IllegalArgumentException('multiple proteins map to one CDS')
            protein.sourceFeature = feature
            contigs[feature.seqId].addProtein(protein)
        }
        if (assigned.size() != features.size()) throw new IllegalArgumentException('GFF CDS lacks a protein sequence')
        genome
    }

    /** Index all explicit aliases; never silently select the first duplicate. */
    static Map<String, List<Feature>> indexAliases(Collection<Feature> features) {
        Map<String, List<Feature>> aliases = [:].withDefault { [] }
        features.each { feature ->
            ['protein_id', 'ID', 'locus_tag'].collectMany { feature.attributes[it] ?: [] }.unique().each {
                aliases[it] << feature
            }
        }
        aliases
    }

    /**
     * Resolve identity evidence only; the caller decides whether incomplete or
     * ambiguous supplied annotations are fatal. Native metadata constrains even
     * an otherwise unique alias, since numeric contig IDs can collide.
     */
    static Set<Feature> matchingFeatures(String proteinId, String header, Map<String, List<Feature>> aliases) {
        def prodigal = header =~ /^\S+\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#\s+(.*)$/
        List<Feature> candidates = new ArrayList<>((aliases[proteinId] ?: []))
        if (prodigal.matches()) {
            Map attributes = GffReader.parseAttributes(prodigal[0][4])
            (attributes['ID'] ?: []).each { candidates.addAll((aliases[it] ?: [])) }
        }
        // A feature can be named by both its protein ID and its native caller ID.
        Set<Feature> unique = Collections.newSetFromMap(new IdentityHashMap<Feature, Boolean>())
        unique.addAll(candidates)
        // Native metadata IDs can name another contig's CDS (e.g. contigs 2,1).
        // Resolve the union using the header coordinates BEFORE requiring uniqueness.
        if (prodigal.matches()) {
            int start = prodigal[0][1].toInteger()
            int end = prodigal[0][2].toInteger()
            Strand strand = prodigal[0][3] == '1' ? Strand.PLUS : Strand.MINUS
            String seqId = proteinId.replaceAll(/_\d+$/, '')
            unique.removeIf { it.start != start || it.end != end || it.strand != strand || it.seqId != seqId }
        }
        unique
    }

    /** Read identity-bearing headers without discarding native caller metadata. */
    static Map<String, String> readHeaders(File file) {
        Map<String, String> out = [:]
        file.withInputStream { stream ->
            InputStream decoded = file.name.endsWith('.gz') ? new java.util.zip.GZIPInputStream(stream) : stream
            new BufferedReader(new InputStreamReader(decoded, 'UTF-8')).eachLine { line ->
                if (!line.trim().startsWith('>')) {
                    if (!out && line.trim() && !line.trim().startsWith(';')) {
                        throw new IllegalArgumentException('sequence data before FASTA header')
                    }
                    return
                }
                String header = line.trim().substring(1).trim()
                String id = header.split(/\s+/)[0]
                if (!id || out.containsKey(id)) throw new IllegalArgumentException('empty or duplicate FASTA identifier')
                out[id] = header
            }
        }
        out
    }

}
