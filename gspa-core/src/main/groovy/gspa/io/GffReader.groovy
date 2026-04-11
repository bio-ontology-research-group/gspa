package gspa.io

import gspa.model.*

/**
 * Reads GFF3 format genome annotation files.
 * Parses features and optionally embedded FASTA sequences (##FASTA directive).
 */
class GffReader {

    /**
     * Read a GFF3 file and populate a Genome with features.
     * If the GFF3 contains ##FASTA, contig sequences are also loaded.
     */
    static Genome readGff3(File gffFile, Genome genome = null) {
        if (genome == null) {
            genome = new Genome(id: gffFile.name.replaceAll(/\.(gff3?|gff)(\.gz)?$/, ''))
        }

        boolean inFasta = false
        String fastaId = null
        StringBuilder fastaSeq = new StringBuilder()

        def reader = gffFile.name.endsWith('.gz') ?
            new BufferedReader(new InputStreamReader(new java.util.zip.GZIPInputStream(new FileInputStream(gffFile)))) :
            gffFile.newReader()

        try {
            reader.eachLine { line ->
                if (inFasta) {
                    if (line.startsWith('>')) {
                        if (fastaId != null) {
                            setContigSequence(genome, fastaId, fastaSeq.toString())
                        }
                        fastaId = line.substring(1).trim().split(/\s+/)[0]
                        fastaSeq = new StringBuilder()
                    } else {
                        fastaSeq.append(line.trim())
                    }
                    return
                }

                if (line == '##FASTA') {
                    inFasta = true
                    return
                }

                if (line.startsWith('#') || line.trim().isEmpty()) return

                def feature = parseFeatureLine(line)
                if (feature) {
                    def contig = genome.findContig(feature.seqId)
                    if (!contig) {
                        contig = new Contig(id: feature.seqId)
                        genome.addContig(contig)
                    }
                    contig.addFeature(feature)
                }
            }
            // Handle last FASTA entry
            if (fastaId != null) {
                setContigSequence(genome, fastaId, fastaSeq.toString())
            }
        } finally {
            reader.close()
        }

        // Sort features on each contig
        genome.contigs.each { it.sortFeatures() }
        genome
    }

    private static void setContigSequence(Genome genome, String contigId, String sequence) {
        def contig = genome.findContig(contigId)
        if (contig) {
            contig.sequence = sequence
        }
    }

    static Feature parseFeatureLine(String line) {
        def fields = line.split('\t')
        if (fields.length < 9) return null

        def feature = new Feature(
            seqId: fields[0],
            source: fields[1],
            type: FeatureType.fromGff3(fields[2]),
            start: fields[3] as int,
            end: fields[4] as int,
            featureScore: fields[5] == '.' ? null : fields[5] as double,
            strand: Strand.fromSymbol(fields[6]),
            phase: fields[7] == '.' ? -1 : fields[7] as int,
            attributes: parseAttributes(fields[8])
        )
        feature
    }

    /**
     * Parse GFF3 attribute column (col 9).
     * Format: key1=value1,value2;key2=value3
     */
    static Map<String, List<String>> parseAttributes(String attrString) {
        Map<String, List<String>> attrs = [:]
        if (!attrString || attrString == '.') return attrs

        attrString.split(';').each { pair ->
            def parts = pair.split('=', 2)
            if (parts.length == 2) {
                def key = URLDecoder.decode(parts[0].trim(), 'UTF-8')
                def values = parts[1].split(',').collect { URLDecoder.decode(it.trim(), 'UTF-8') }
                attrs[key] = values
            }
        }
        attrs
    }
}
