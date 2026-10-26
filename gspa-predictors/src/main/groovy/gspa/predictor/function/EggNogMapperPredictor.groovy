package gspa.predictor.function

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * eggNOG-mapper predictor for ortholog group assignment.
 * Provides COG, KEGG, GO annotations via orthology inference.
 */
class EggNogMapperPredictor extends AbstractToolPredictor {

    /** eggNOG database directory */
    String dataDir

    /** Search mode: diamond, mmseqs, hmmer */
    String searchMode = 'diamond'

    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'eggnog-mapper' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.GO, AnnotationType.EC, AnnotationType.COG, AnnotationType.KEGG] as Set
    }

    @Override
    String getExecutable() { 'emapper.py' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def prefix = new File(outputDir, 'eggnog_results')
        def cmd = [
            executablePath ?: executable,
            '-i', inputFasta.absolutePath,
            '-o', prefix.absolutePath,
            '--cpu', threads.toString(),
            '-m', searchMode,
            '--itype', 'proteins',
            '--go_evidence', 'non-electronic',
        ]
        if (dataDir) {
            cmd.addAll(['--data_dir', dataDir])
        }
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def annotFile = new File(outputDir, 'eggnog_results.emapper.annotations')
        // emapper always writes this file (header only when nothing matched), so
        // its absence after a clean exit is a failure, not "zero annotations".
        if (!annotFile.exists()) {
            throw new IllegalStateException("eggnog-mapper produced no ${annotFile.name}")
        }

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        // emapper 2.1.x layout, used only when the file carries no '#query' header:
        // 0 query, 1 seed_ortholog, 2 evalue, 3 score (bit score), 4 eggNOG_OGs,
        // 5 max_annot_lvl, 6 COG_category, 7 Description, 8 Preferred_name, 9 GOs,
        // 10 EC, 11 KEGG_ko ...
        Map<String, Integer> col = [evalue: 2, score: 3, ogs: 4, cog: 6, gos: 9, ec: 10, ko: 11]
        boolean headerSeen = false

        annotFile.eachLine { line ->
            if (line.startsWith('#query')) {
                col = columnsFromHeader(line)
                headerSeen = true
                return
            }
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t', -1)
            int needed = col.values().max() + 1
            if (fields.length < needed || (!headerSeen && fields.length < 21)) return

            // Without a header, validate that the assumed E-value / bit-score
            // columns are numeric before trusting them; otherwise skip the row
            // rather than invent confidences from the wrong column.
            if (!headerSeen && (!isNumeric(fields[col.evalue]) || !isNumeric(fields[col.score]))) return

            String proteinId = fields[0]
            String eggNogOGs = field(fields, col.ogs)    // Orthologous groups
            String goTerms = field(fields, col.gos)      // GO terms
            String ecNumbers = field(fields, col.ec)     // EC numbers
            String keggKOs = field(fields, col.ko)       // KEGG KO
            String cogCats = field(fields, col.cog)      // COG categories
            // Confidence derives from the seed-ortholog E-value (smaller = stronger);
            // an E-value of exactly 0 is the strongest possible hit.
            String evalueRaw = field(fields, col.evalue)
            double confScore
            if (!isNumeric(evalueRaw)) {
                confScore = 0.7
            } else {
                double evalue = evalueRaw as double
                confScore = evalue <= 0 ? 1.0 : Math.max(0.0, Math.min(1.0, -Math.log10(evalue) / 30.0))
            }

            // GO annotations
            if (goTerms && goTerms != '-') {
                goTerms.split(',').each { goTerm ->
                    goTerm = goTerm.trim()
                    if (goTerm.startsWith('GO:')) {
                        results[proteinId] << new Annotation(
                            type: AnnotationType.GO, value: goTerm,
                            source: name, score: confScore, evidence: 'IEA',
                            metadata: [orthologGroup: eggNogOGs]
                        )
                    }
                }
            }

            // EC numbers
            if (ecNumbers && ecNumbers != '-') {
                ecNumbers.split(',').each { ec ->
                    ec = ec.trim()
                    if (ec) {
                        results[proteinId] << new Annotation(
                            type: AnnotationType.EC, value: "EC:${ec}",
                            source: name, score: confScore,
                            metadata: [orthologGroup: eggNogOGs]
                        )
                    }
                }
            }

            // KEGG KO
            if (keggKOs && keggKOs != '-') {
                keggKOs.split(',').each { ko ->
                    ko = ko.trim()
                    if (ko) {
                        results[proteinId] << new Annotation(
                            type: AnnotationType.KEGG, value: ko,
                            source: name, score: confScore,
                        )
                    }
                }
            }

            // COG categories
            if (cogCats && cogCats != '-' && cogCats != 'S') {
                results[proteinId] << new Annotation(
                    type: AnnotationType.COG, value: cogCats,
                    source: name, score: confScore,
                    metadata: [orthologGroup: eggNogOGs]
                )
            }
        }
        results
    }

    /**
     * Resolve column indices from the '#query' header. Older emapper releases
     * name the seed columns seed_evalue / seed_score, newer ones evalue / score.
     */
    private static Map<String, Integer> columnsFromHeader(String headerLine) {
        List<String> names = headerLine.substring(1).split('\t', -1).collect { it.trim() }
        Map<String, List<String>> wanted = [
            evalue: ['evalue', 'seed_evalue'], score: ['score', 'seed_score'],
            ogs: ['eggNOG_OGs'], cog: ['COG_category'], gos: ['GOs'], ec: ['EC'], ko: ['KEGG_ko'],
        ]
        Map<String, Integer> col = [:]
        wanted.each { key, aliases ->
            int idx = names.findIndexOf { it in aliases }
            if (idx < 0) {
                throw new IllegalStateException(
                    "eggnog-mapper annotations header lacks required column ${aliases.join('/')}: ${names}")
            }
            col[key] = idx
        }
        col
    }

    private static String field(String[] fields, int idx) {
        idx < fields.length ? fields[idx] : ''
    }

    private static boolean isNumeric(String s) {
        if (s == null || s.isEmpty() || s == '-') return false
        try { Double.parseDouble(s); return true } catch (NumberFormatException e) { return false }
    }
}
