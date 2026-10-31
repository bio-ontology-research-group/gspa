package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

/** Protein-only AMRFinderPlus calls, including the explicitly labelled plus categories. */
class AmrFinderPredictor extends AbstractToolPredictor {
    String organism
    String database
    int threads = Runtime.runtime.availableProcessors()
    protected String executedVersion = 'unknown'

    @Override String getName() { 'amrfinder' }
    @Override Set<AnnotationType> getOutputTypes() { [AnnotationType.AMR] as Set }
    @Override String getExecutable() { 'amrfinder' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (threads < 1) throw new IllegalArgumentException('amrfinder threads must be positive')
        executedVersion = getVersion()
        def cmd = [executablePath ?: executable, '-p', inputFasta.absolutePath,
                   '-o', new File(outputDir, 'amrfinder_results.tsv').absolutePath,
                   '--threads', threads.toString(), '--plus']
        if (database) cmd.addAll(['-d', database])
        if (organism) cmd.addAll(['--organism', organism])
        cmd
    }

    @Override
    protected String successNote() {
        'Protein-only search; AMR/STRESS/VIRULENCE categories are tool calls, not resistance phenotypes. ' +
            'Identity/100 is an uncalibrated similarity score; HMM calls without identity remain unscored.'
    }

    // Both pre-v4 and v4 published names. No positional fallback: coordinates and extra
    // columns are optional, and a missing header must never become a fabricated hit.
    private static final Map<String, List<String>> COLUMNS = [
        protein: ['Protein id', 'Protein identifier'], gene: ['Element symbol', 'Gene symbol'],
        description: ['Element name', 'Sequence name'], scope: ['Scope'],
        type: ['Type', 'Element type'], subtype: ['Subtype', 'Element subtype'],
        drugClass: ['Class'], subclass: ['Subclass'], method: ['Method'],
        coverage: ['% Coverage of reference', '% Coverage of reference sequence'],
        identity: ['% Identity to reference', '% Identity to reference sequence'],
        hmmAccession: ['HMM accession', 'HMM id'], hmmDescription: ['HMM description'],
    ]

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def file = requireOutput(outputDir, 'amrfinder_results.tsv', 'AMRFinderPlus always writes a header')
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        List<String> header = null
        Map<String, Integer> columns = [:]
        int lineNumber = 0
        file.eachLine('UTF-8') { line ->
            lineNumber++
            if (!line.trim()) throw invalid(lineNumber, 'blank row')
            def fields = line.split('\t', -1).toList()
            if (header == null) {
                header = fields
                if (header.toSet().size() != header.size() || header.any { !it }) {
                    throw invalid(lineNumber, 'empty or duplicate column name')
                }
                AmrFinderPredictor.COLUMNS.each { key, aliases ->
                    def matches = aliases.findAll { header.contains(it) }
                    if (matches.size() != 1) throw invalid(lineNumber, "missing or ambiguous column ${aliases}")
                    columns[key] = header.indexOf(matches[0])
                }
                return
            }
            if (fields.size() != header.size()) throw invalid(lineNumber, "expected ${header.size()} fields, found ${fields.size()}")
            def val = { String key -> fields[columns[key]] }
            ['protein', 'gene', 'description', 'scope', 'type', 'subtype', 'method'].each { key ->
                if (!val(key)?.trim() || val(key) == 'NA') throw invalid(lineNumber, "missing ${key}")
            }
            if (!(val('scope') in ['core', 'plus'])) throw invalid(lineNumber, 'invalid Scope')
            if (!(val('type') in ['AMR', 'STRESS', 'VIRULENCE'])) throw invalid(lineNumber, 'invalid Type')
            boolean hmm = val('method') == 'HMM'
            Double identity = percentage(val('identity'), hmm, lineNumber, 'identity')
            Double coverage = percentage(val('coverage'), hmm, lineNumber, 'coverage')
            if (hmm && (!val('hmmAccession') || val('hmmAccession') == 'NA')) {
                throw invalid(lineNumber, 'HMM call lacks its HMM accession')
            }
            Map raw = [:]
            header.eachWithIndex { String key, int i -> raw[key] = fields[i] }
            results[val('protein')] << new Annotation(
                type: AnnotationType.AMR, value: "${val('type')}:${val('subtype')}", source: name,
                // Annotation.score is a primitive. Zero is only a storage placeholder when
                // score_available=false; such observations cannot enter functional inference.
                score: identity == null ? 0d : identity / 100d,
                metadata: [geneName: val('gene'), sequenceName: val('description'), scope: val('scope'),
                    drugClass: val('drugClass'), drugSubclass: val('subclass'), method: val('method'),
                    coverage: coverage, identity: identity, hmmAccession: val('hmmAccession'),
                    hmmDescription: val('hmmDescription'), score_available: identity != null,
                    score_model: 'reference_identity_fraction/1', score_calibrated: false,
                    tool_version: executedVersion, rawFields: raw])
        }
        results
    }

    protected static Double percentage(String value, boolean allowMissing, int line, String field) {
        if (value == 'NA' && allowMissing) return null
        double number
        try { number = Double.parseDouble(value) }
        catch (NumberFormatException ignored) { throw invalid(line, "invalid ${field}: ${value}") }
        if (!Double.isFinite(number) || number < 0d || number > 100d) {
            throw invalid(line, "out-of-range ${field}: ${value}")
        }
        number
    }

    protected static ToolOutputContractException invalid(int line, String detail) {
        new ToolOutputContractException("amrfinder row ${line}: ${detail}")
    }
}
