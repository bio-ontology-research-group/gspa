package gspa.predictor.similarity

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * DIAMOND sequence similarity predictor.
 * Runs DIAMOND blastp against a reference database (e.g., UniProt/SwissProt)
 * and transfers GO/EC annotations from best hits.
 */
class DiamondPredictor extends AbstractToolPredictor {

    /** Path to DIAMOND database (.dmnd) */
    String database

    /** E-value threshold */
    double evalue = 1e-5

    /** Maximum target sequences to report */
    int maxTargetSeqs = 10

    /** Minimum query coverage */
    double queryCover = 50.0

    /** Minimum subject coverage */
    double subjectCover = 50.0

    /** Minimum percent identity */
    double minIdentity = 30.0

    /** Number of threads */
    int threads = Runtime.runtime.availableProcessors()

    /**
     * UniProt accession -> GO terms mapping.
     * Loaded from idmapping or GOA file. If null, GO terms are extracted from hit titles.
     */
    Map<String, Set<String>> accessionToGO = null

    /**
     * Load UniProt accession -> GO mapping from a TSV file.
     * Format: accession\tGO:XXXXXXX\tGO:XXXXXXX\t...
     * Or GAF-derived: accession\tGO_term (one per line)
     */
    void loadGoMapping(File mappingFile) {
        accessionToGO = [:].withDefault { [] as Set }
        mappingFile.eachLine { line ->
            if (line.startsWith('#') || line.startsWith('!')) return
            def fields = line.split('\t')
            if (fields.length >= 2) {
                String acc = fields[0].trim()
                for (int i = 1; i < fields.length; i++) {
                    String val = fields[i].trim()
                    if (val.startsWith('GO:')) {
                        accessionToGO[acc] << val
                    }
                }
            }
        }
        log.info("Loaded GO mappings for ${accessionToGO.size()} accessions")
    }

    @Override
    String getName() { 'diamond' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.GO, AnnotationType.EC] as Set
    }

    @Override
    String getExecutable() { 'diamond' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def outputFile = new File(outputDir, 'diamond_results.tsv')
        [
            executablePath ?: executable,
            'blastp',
            '--db', new File(database).absolutePath,
            '--query', inputFasta.absolutePath,
            '--out', outputFile.absolutePath,
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
                         'evalue', 'bitscore', 'stitle',
            '--evalue', evalue.toString(),
            '--max-target-seqs', maxTargetSeqs.toString(),
            '--query-cover', queryCover.toString(),
            '--subject-cover', subjectCover.toString(),
            '--id', minIdentity.toString(),
            '--threads', threads.toString(),
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def outputFile = new File(outputDir, 'diamond_results.tsv')
        if (!outputFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        outputFile.eachLine { line ->
            def fields = line.split('\t')
            if (fields.length < 9) return

            String queryId = fields[0]
            String subjectId = fields[1]
            double pident = fields[2] as double
            double eval = fields[6] as double
            double bitscore = fields[7] as double
            String title = fields[8]

            // Confidence score: normalize by identity
            double score = Math.min(1.0, pident / 100.0)

            // Extract UniProt accession from subject ID (sp|P03707|TERS_LAMBD -> P03707)
            String accession = extractAccession(subjectId)

            // Get GO terms: from mapping file if available, otherwise from title
            def goTerms = (accessionToGO != null && accession) ?
                (accessionToGO[accession] ?: []) :
                extractGoTerms(title)

            goTerms.each { goTerm ->
                results[queryId] << new Annotation(
                    type: AnnotationType.GO,
                    value: goTerm,
                    source: name,
                    score: score,
                    evidence: 'IEA',
                    metadata: [
                        hit: subjectId,
                        hitAccession: accession,
                        pident: pident,
                        evalue: eval,
                        bitscore: bitscore,
                    ]
                )
            }

            // Extract EC numbers if present in title
            def ecNumbers = extractEcNumbers(title)
            ecNumbers.each { ec ->
                results[queryId] << new Annotation(
                    type: AnnotationType.EC,
                    value: ec,
                    source: name,
                    score: score,
                    metadata: [hit: subjectId, pident: pident]
                )
            }

            // If no GO terms found, create a CUSTOM annotation with the product name
            // from the best hit (useful for product naming even without GO)
            if (goTerms.isEmpty() && ecNumbers.isEmpty()) {
                String product = extractProduct(title)
                if (product) {
                    results[queryId] << new Annotation(
                        type: AnnotationType.CUSTOM,
                        value: product,
                        source: name,
                        score: score,
                        evidence: 'ISS',
                        metadata: [
                            hit: subjectId,
                            hitAccession: accession,
                            pident: pident,
                            evalue: eval,
                            annotationType: 'product',
                        ]
                    )
                }
            }
        }

        results
    }

    private static List<String> extractGoTerms(String text) {
        (text =~ /GO:\d{7}/).collect { it as String }
    }

    private static List<String> extractEcNumbers(String text) {
        (text =~ /EC:[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+/).collect { it as String }
    }

    /**
     * Extract protein product name from UniProt FASTA title.
     * E.g., "sp|P03707|TERS_LAMBD Terminase small subunit OS=..." -> "Terminase small subunit"
     */
    private static String extractProduct(String title) {
        // UniProt format: sp|ACC|NAME Description OS=Organism ...
        def matcher = title =~ /(?:sp|tr)\|\S+\|\S+\s+(.+?)\s+OS=/
        if (matcher.find()) {
            return matcher.group(1).trim()
        }
        // Fallback: take everything before OS= or the whole title
        def osIdx = title.indexOf(' OS=')
        if (osIdx > 0) {
            // Strip the leading sp|...|... part
            def spaceIdx = title.indexOf(' ')
            if (spaceIdx > 0 && spaceIdx < osIdx) {
                return title.substring(spaceIdx + 1, osIdx).trim()
            }
        }
        null
    }

    /**
     * Extract UniProt accession from subject ID.
     * Handles formats: sp|P03707|TERS_LAMBD, tr|A0A0...|..., or plain P03707
     */
    private static String extractAccession(String subjectId) {
        def parts = subjectId.split('\\|')
        if (parts.length >= 2) return parts[1]
        return subjectId
    }
}
