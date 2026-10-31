package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

/**
 * dbCAN predictor for carbohydrate-active enzyme (CAZyme) annotation.
 * Preserves named-method family calls; overview tool counts are not probabilities.
 * The command contract is the legacy run_dbcan CLI, not the dbcan v5 subcommands.
 */
class DbCanPredictor extends AbstractToolPredictor {

    /** dbCAN database directory */
    String dbDir

    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'dbcan' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.CAZYME] as Set }

    @Override
    String getExecutable() { 'run_dbcan' }

    @Override
    boolean isAvailable() {
        // run_dbcan 4.x defines --version as an argument with a value, not a
        // version-reporting action. Its bare --version exits with a usage error.
        try {
            def process = [executablePath ?: executable, '--help'].execute()
            StringBuffer stdout = new StringBuffer(), stderr = new StringBuffer()
            Thread outputReader = process.consumeProcessOutputStream(stdout)
            Thread errorReader = process.consumeProcessErrorStream(stderr)
            process.waitForOrKill(10000)
            outputReader.join(1000)
            errorReader.join(1000)
            process.exitValue() == 0 && ['--db_dir', '--hmm_cpu', '--dia_cpu', '--dbcan_thread'].every {
                stdout.toString().contains(it)
            }
        } catch (Exception ignored) {
            false
        }
    }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (threads < 1) throw new IllegalArgumentException('dbCAN threads must be positive')
        def cmd = [
            executablePath ?: executable,
            inputFasta.absolutePath,
            'protein',
            '--out_dir', outputDir.absolutePath,
            '--hmm_cpu', threads.toString(), '--dia_cpu', threads.toString(),
            '--dbcan_thread', threads.toString(),
        ]
        if (dbDir) cmd.addAll(['--db_dir', dbDir])
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def overviewFile = requireOutput(outputDir, 'overview.txt',
            'dbCAN overview requires a named-column header, including successful no-hit output')
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        List<String> header = null
        Map<String, String> methodColumns = [:]
        Set<String> seen = [] as Set
        int lineNumber = 0
        overviewFile.eachLine { line ->
            lineNumber++
            if (!line.trim()) return
            List<String> fields = line.split('\t', -1).collect { it.trim() }
            if (header == null && line.startsWith('#') && !fields.contains('Gene ID')) return
            if (header == null) {
                header = fields
                if (header.toSet().size() != header.size() ||
                    !header.containsAll(['Gene ID', 'DIAMOND', '#ofTools'])) {
                    throw invalid(lineNumber, 'missing or duplicated named columns')
                }
                [['HMMER', 'dbCAN', 'dbCAN_hmm'], ['dbCAN_sub', 'eCAMI', 'Hotpep']].eachWithIndex { aliases, i ->
                    def present = aliases.findAll { header.contains(it) }
                    if (present.size() != 1) throw invalid(lineNumber, 'missing or ambiguous family-method column')
                    methodColumns[i == 0 ? 'HMMER' : present[0]] = present[0]
                }
                methodColumns['DIAMOND'] = 'DIAMOND'
                return
            }
            if (fields.size() != header.size()) throw invalid(lineNumber, 'row width differs from header')
            Map<String, String> raw = [:]
            header.eachWithIndex { key, i -> raw[key] = fields[i] }
            String proteinId = raw['Gene ID']
            if (!proteinId || !seen.add(proteinId)) throw invalid(lineNumber, 'empty or duplicate protein ID')
            if (!(raw['#ofTools'] ==~ /[0-3]/)) throw invalid(lineNumber, 'invalid number of tools')
            int toolCount = Integer.parseInt(raw['#ofTools'])
            int observedTools = methodColumns.values().count { raw[it] != '-' && !raw[it].isEmpty() }
            if (toolCount != observedTools) throw invalid(lineNumber, 'tool count differs from nonempty method columns')
            Map<String, List<String>> supportingMethods = [:].withDefault { [] }
            methodColumns.each { method, column ->
                String value = raw[column]
                if (!value) throw invalid(lineNumber, 'empty method value; expected explicit dash for no hit')
                if (value == '-') return
                value.split(/\+/, -1).each { token ->
                    // CAZy subfamilies (_4) and dbCAN/eCAMI subfamilies (_e10) are distinct.
                    // Keep the latter in rawFields, and emit their parent CAZy family.
                    def match = token.trim() =~ /^((?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?)(?:_e\d+)?(?:\((\d+)-(\d+)\))?$/
                    if (!match.matches()) throw invalid(lineNumber, "invalid ${method} family token")
                    if (match.group(2) != null &&
                        (new BigInteger(match.group(2)) < 1 || new BigInteger(match.group(3)) < new BigInteger(match.group(2)))) {
                        throw invalid(lineNumber, 'invalid family domain interval')
                    }
                    String family = match.group(1)
                    if (!supportingMethods[family].contains(method)) supportingMethods[family] << method
                }
            }
            supportingMethods.each { family, methods ->
                results[proteinId] << new Annotation(type: AnnotationType.CAZYME,
                    value: family, source: name,
                    // Primitive storage placeholder; score_available=false excludes confidence inference.
                    score: 0d, metadata: [toolsWithHits: toolCount, supportingMethods: methods,
                        score_available: false, score_calibrated: false,
                        score_model: 'unscored_family_calls/1', rawFields: new LinkedHashMap(raw)])
            }
        }
        if (header == null) throw invalid(0, 'missing overview header')
        results
    }

    @Override
    protected String successNote() {
        'CAZyme family calls are unscored observations. Tool counts do not establish family agreement, enzyme specificity or calibrated confidence.'
    }

    private static ToolOutputContractException invalid(int line, String reason) {
        new ToolOutputContractException("dbCAN overview line ${line}: ${reason}")
    }
}
