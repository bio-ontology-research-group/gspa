package gspa.predictor.function

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

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
        // emapper always writes this file (a '#query' header and comment footer even
        // when nothing matched), so absence or emptiness after a clean exit is a
        // failed contract, not "zero annotations".
        def annotFile = requireOutput(outputDir, 'eggnog_results.emapper.annotations',
            'emapper writes the annotations table with a #query header on every run')

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        // emapper 2.1.x layout, used only when the file carries no '#query' header:
        // 0 query, 1 seed_ortholog, 2 evalue, 3 score (bit score), 4 eggNOG_OGs,
        // 5 max_annot_lvl, 6 COG_category, 7 Description, 8 Preferred_name, 9 GOs,
        // 10 EC, 11 KEGG_ko ...
        Map<String, Integer> col = [evalue: 2, score: 3, ogs: 4, cog: 6, gos: 9, ec: 10, ko: 11]
        boolean headerSeen = false
        int dataRows = 0
        int lineNo = 0

        annotFile.eachLine { line ->
            lineNo++
            if (line.startsWith('#query')) {
                col = columnsFromHeader(line)
                headerSeen = true
                return
            }
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t', -1)
            // A data row that is too short or whose E-value is not a finite number >= 0
            // means the table is corrupt or its columns are not what we think they are.
            // Rejecting the file is the only honest outcome: no confidence is invented.
            int needed = headerSeen ? col.values().max() + 1 : 21
            if (fields.length < needed) {
                throw new ToolOutputContractException("eggnog-mapper line ${lineNo}: ${fields.length} columns, " +
                    "expected at least ${needed}${headerSeen ? '' : ' (no #query header found; assumed emapper 2.1 layout)'}")
            }
            dataRows++

            String proteinId = fields[0]
            String eggNogOGs = field(fields, col.ogs)    // Orthologous groups
            String goTerms = field(fields, col.gos)      // GO terms
            String ecNumbers = field(fields, col.ec)     // EC numbers
            String keggKOs = field(fields, col.ko)       // KEGG KO
            String cogCats = field(fields, col.cog)      // COG categories
            // Confidence derives from the seed-ortholog E-value (smaller = stronger);
            // an E-value of exactly 0 is valid and is the strongest possible hit.
            double evalue = parseEvalue(field(fields, col.evalue), lineNo, proteinId, headerSeen)
            double confScore = evalue == 0.0d ? 1.0d : Math.min(1.0d, -Math.log10(evalue) / 30.0d)
            if (confScore < 0.0d) confScore = 0.0d   // E-values above 1 are weak hits, not errors

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
        if (!headerSeen && dataRows == 0) {
            throw new ToolOutputContractException(
                'eggnog-mapper annotations file has neither a #query header nor data rows')
        }
        results
    }

    /**
     * Resolve column indices from the '#query' header. Older emapper releases
     * name the seed columns seed_evalue / seed_score, newer ones evalue / score.
     */
    protected static Map<String, Integer> columnsFromHeader(String headerLine) {
        List<String> names = headerLine.substring(1).split('\t', -1).collect { it.trim() }
        Map<String, List<String>> wanted = [
            evalue: ['evalue', 'seed_evalue'], score: ['score', 'seed_score'],
            ogs: ['eggNOG_OGs'], cog: ['COG_category'], gos: ['GOs'], ec: ['EC'], ko: ['KEGG_ko'],
        ]
        Map<String, Integer> col = [:]
        wanted.each { key, aliases ->
            int idx = names.findIndexOf { it in aliases }
            if (idx < 0) {
                throw new ToolOutputContractException(
                    "eggnog-mapper annotations header lacks required column ${aliases.join('/')}: ${names}")
            }
            col[key] = idx
        }
        col
    }

    protected static String field(String[] fields, int idx) {
        idx < fields.length ? fields[idx] : ''
    }

    /** Strict E-value: a finite number >= 0 (0 is valid). Anything else is a contract violation. */
    protected static double parseEvalue(String raw, int lineNo, String proteinId, boolean headerSeen) {
        double v
        try {
            v = Double.parseDouble(raw?.trim())
        } catch (NumberFormatException | NullPointerException e) {
            throw badEvalue(raw, lineNo, proteinId, headerSeen)
        }
        if (!Double.isFinite(v) || v < 0.0d) throw badEvalue(raw, lineNo, proteinId, headerSeen)
        v
    }

    protected static ToolOutputContractException badEvalue(String raw, int lineNo, String proteinId, boolean headerSeen) {
        new ToolOutputContractException(
            "eggnog-mapper line ${lineNo} (${proteinId}): E-value '${raw}' is not a finite number >= 0" +
            (headerSeen ? ' in the column named evalue/seed_evalue' :
                ' in column 3 (no #query header, so the emapper 2.1 layout was assumed; check the emapper version)'))
    }
}
