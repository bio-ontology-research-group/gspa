package gspa.io

import gspa.model.Contig
import gspa.model.Feature
import gspa.model.FeatureType
import gspa.model.Genome
import gspa.model.Protein
import gspa.model.Strand
import groovy.util.logging.Slf4j

/**
 * Translates CDS features into proteins, given contig nucleotide sequences.
 *
 * This is the bridge that lets {@code gspa annotate} accept a genome (or a
 * multi-sequence/metagenome FASTA) together with a GFF3 annotation: rather
 * than re-calling genes, GSPA reuses the supplied CDS coordinates and
 * translates them in place. The same path serves a single GFF3 that carries
 * its contig sequences in a {@code ##FASTA} block.
 *
 * Translation uses the standard genetic code; the amino-acid assignment for
 * sense codons is identical to the NCBI bacterial table (11), so the same
 * map serves prokaryotes and eukaryotes. Alternative initiator codons
 * (GTG/TTG/CTG/ATT/ATC/ATA) are emitted as Methionine when they start a CDS.
 */
@Slf4j
class CdsTranslator {

    /** Standard genetic code (NCBI table 1; sense codons match table 11). */
    static final Map<String, Character> CODON_TABLE = buildCodonTable()

    /** Codons treated as Methionine when they begin a CDS (bacterial table 11). */
    static final Set<String> ALT_START_CODONS =
        ['ATG', 'GTG', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA'] as Set

    private static Map<String, Character> buildCodonTable() {
        def aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        def bases = ['T', 'C', 'A', 'G']
        def table = [:]
        int i = 0
        for (b1 in bases) {
            for (b2 in bases) {
                for (b3 in bases) {
                    table[(b1 + b2 + b3)] = aas.charAt(i++)
                }
            }
        }
        table
    }

    /** Reverse-complement a nucleotide string (IUPAC-aware for the common bases). */
    static String reverseComplement(String seq) {
        def sb = new StringBuilder(seq.length())
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(complement(seq.charAt(i)))
        }
        sb.toString()
    }

    private static char complement(char c) {
        switch (Character.toUpperCase(c)) {
            case 'A': return 'T'
            case 'T': return 'A'
            case 'U': return 'A'
            case 'G': return 'C'
            case 'C': return 'G'
            case 'R': return 'Y'
            case 'Y': return 'R'
            case 'S': return 'S'
            case 'W': return 'W'
            case 'K': return 'M'
            case 'M': return 'K'
            case 'B': return 'V'
            case 'V': return 'B'
            case 'D': return 'H'
            case 'H': return 'D'
            default: return 'N'
        }
    }

    /**
     * Translate a coding nucleotide sequence (already in the correct
     * orientation, 5'->3') into an amino-acid string.
     *
     * @param cds                 coding sequence, 5'->3'
     * @param phase               bases to skip before the first complete codon (0/1/2; -1 treated as 0)
     * @param translateStartAsMet emit M for an alternative initiator codon
     * @param trimStop            drop a single trailing stop codon ('*')
     */
    static String translate(String cds, int phase = 0, boolean translateStartAsMet = true,
                            boolean trimStop = true) {
        if (!cds) return ''
        String seq = cds.toUpperCase()
        int offset = phase > 0 ? phase : 0
        def sb = new StringBuilder((seq.length() - offset).intdiv(3) + 1)
        boolean first = true
        for (int i = offset; i + 3 <= seq.length(); i += 3) {
            String codon = seq.substring(i, i + 3)
            char aa
            if (first && translateStartAsMet && ALT_START_CODONS.contains(codon)) {
                aa = 'M'
            } else {
                aa = CODON_TABLE.getOrDefault(codon, 'X' as char)
            }
            sb.append(aa)
            first = false
        }
        String protein = sb.toString()
        if (trimStop && protein.endsWith('*')) {
            protein = protein.substring(0, protein.length() - 1)
        }
        protein
    }

    /**
     * Translate the coding sequence assembled from one or more CDS segments
     * (multi-exon genes are concatenated in transcription order).
     *
     * @param segments CDS features (all on the same strand, sharing a protein)
     * @param contigSeq the contig nucleotide sequence (1-based GFF coordinates)
     */
    static String translateCds(List<Feature> segments, String contigSeq,
                               boolean translateStartAsMet = true) {
        if (!segments || !contigSeq) return ''
        Strand strand = segments.first().strand
        // Transcription order: ascending start on +, descending on -.
        def ordered = strand == Strand.MINUS ?
            segments.sort(false) { -it.start } :
            segments.sort(false) { it.start }

        def coding = new StringBuilder()
        for (Feature seg : ordered) {
            int s = Math.max(1, seg.start)
            int e = Math.min(contigSeq.length(), seg.end)
            if (e < s) continue
            String piece = contigSeq.substring(s - 1, e)
            if (strand == Strand.MINUS) piece = reverseComplement(piece)
            coding.append(piece)
        }
        int phase = ordered.first().phase
        translate(coding.toString(), phase, translateStartAsMet, true)
    }

    /**
     * Populate {@code genome.contigs[*].proteins} by translating the CDS
     * features already loaded on each contig. Requires each contig to carry
     * its nucleotide {@code sequence} (from the genome FASTA or an embedded
     * {@code ##FASTA} block). CDS segments are grouped per protein by the
     * GFF3 {@code ID} attribute (NCBI multi-exon rows share an ID; prodigal
     * gives each CDS a unique ID), falling back to locus_tag.
     *
     * @return number of proteins created
     */
    static int populateProteins(Genome genome, boolean translateStartAsMet = true) {
        int created = 0
        Set<String> usedIds = new HashSet<>()
        genome.contigs.each { Contig contig ->
            if (!contig.sequence) {
                if (contig.cdsFeatures()) {
                    log.warn("Contig ${contig.id} has CDS features but no sequence; " +
                        "cannot translate. Provide the genome FASTA or an embedded ##FASTA block.")
                }
                return
            }
            // Group CDS segments per protein, preserving first-seen order.
            Map<String, List<Feature>> byProtein = new LinkedHashMap<>()
            contig.cdsFeatures().each { Feature f ->
                String key = proteinKey(f)
                byProtein.computeIfAbsent(key) { [] } << f
            }
            byProtein.each { String key, List<Feature> segments ->
                String aa = translateCds(segments, contig.sequence, translateStartAsMet)
                if (!aa) return
                String pid = uniqueId(key, contig.id, usedIds)
                Feature primary = segments.min { it.start }
                def protein = new Protein(id: pid, sequence: aa, sourceFeature: primary)
                contig.addProtein(protein)
                created++
            }
        }
        created
    }

    /** Grouping/identity key for a CDS feature: protein_id, ID, locus_tag, or synthesized. */
    static String proteinKey(Feature f) {
        f.attributes['protein_id']?.first() ?:
            f.id ?:
            f.locusTag ?:
            "${f.seqId}_${f.start}_${f.end}_${f.strand.symbol}"
    }

    private static String uniqueId(String key, String contigId, Set<String> used) {
        String id = key
        if (used.contains(id)) {
            id = "${contigId}:${key}"
            int n = 2
            while (used.contains(id)) {
                id = "${contigId}:${key}.${n++}"
            }
        }
        used.add(id)
        id
    }
}
