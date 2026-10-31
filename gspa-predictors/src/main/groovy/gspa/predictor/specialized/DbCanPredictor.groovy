package gspa.predictor.specialized

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.io.FastaReader
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.ToolOutputContractException

/**
 * dbCAN predictor for carbohydrate-active enzyme (CAZyme) annotation.
 * Preserves named-method family calls; overview tool counts are not probabilities.
 * Supports the run_dbcan 4.x positional and 5.2.8 subcommand interfaces.
 */
class DbCanPredictor extends AbstractToolPredictor {

    /** dbCAN database directory */
    String dbDir
    /** Supported command contracts: run_dbcan 4.x positional CLI or 5.2.8 subcommands. */
    String interfaceVersion = '4'

    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'dbcan' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.CAZYME] as Set }

    @Override
    String getExecutable() { 'run_dbcan' }

    @Override
    String getVersion() {
        // Legacy --version selects a database version; it does not report the tool version.
        if (interfaceVersion != '5') return 'unknown'
        try {
            def process = [executablePath ?: executable, 'version'].execute()
            StringBuffer stdout = new StringBuffer(), stderr = new StringBuffer()
            Thread reader = process.consumeProcessOutputStream(stdout)
            Thread errors = process.consumeProcessErrorStream(stderr)
            process.waitForOrKill(10000)
            reader.join(1000)
            errors.join(1000)
            def match = stdout.toString().trim() =~ /^dbCAN version: (\S+)$/
            return process.exitValue() == 0 && match.matches() ? match.group(1) : 'unknown'
        } catch (Exception ignored) {
            return 'unknown'
        }
    }

    @Override
    boolean isAvailable() {
        // run_dbcan 4.x defines --version as an argument with a value, not a
        // version-reporting action. Its bare --version exits with a usage error.
        try {
            validateInterface()
            def command = [executablePath ?: executable]
            if (interfaceVersion == '5') command << 'CAZyme_annotation'
            def process = (command + ['--help']).execute()
            StringBuffer stdout = new StringBuffer(), stderr = new StringBuffer()
            Thread outputReader = process.consumeProcessOutputStream(stdout)
            Thread errorReader = process.consumeProcessErrorStream(stderr)
            process.waitForOrKill(10000)
            outputReader.join(1000)
            errorReader.join(1000)
            List<String> flags = interfaceVersion == '5' ?
                ['--db_dir', '--input_raw_data', '--mode', '--threads'] :
                ['--db_dir', '--hmm_cpu', '--dia_cpu', '--dbcan_thread']
            process.exitValue() == 0 && flags.every {
                stdout.toString().contains(it)
            }
        } catch (Exception ignored) {
            false
        }
    }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        validateInterface()
        if (threads < 1) throw new IllegalArgumentException('dbCAN threads must be positive')
        if (!dbDir) throw new IllegalArgumentException('dbCAN requires an explicit database directory')
        if (interfaceVersion == '5') {
            return [executablePath ?: executable, 'CAZyme_annotation',
                '--input_raw_data', inputFasta.absolutePath, '--mode', 'protein',
                '--output_dir', outputDir.absolutePath, '--db_dir', dbDir,
                '--threads', threads.toString()]
        }
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
        validateInterface()
        def overviewFile = requireOutput(outputDir, interfaceVersion == '5' ? 'overview.tsv' : 'overview.txt',
            'dbCAN overview requires a named-column header, including successful no-hit output')
        File input = new File(outputDir, 'input.faa')
        parseOverview(overviewFile, input.isFile() ? FastaReader.readProteins(input)*.id as Set : null)
    }

    /** Shared parser for native execution and the Nextflow/import route. */
    Map<String, List<Annotation>> parseOverview(File overviewFile, Set<String> allowedProteinIds = null) {
        if (!overviewFile.isFile() || overviewFile.length() == 0L) throw invalid(0, 'missing or empty overview')
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
            if (allowedProteinIds != null && !allowedProteinIds.contains(proteinId)) {
                throw invalid(lineNumber, 'unknown protein ID: ' + proteinId)
            }
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
                    score: 0d, metadata: [proteinToolCount: toolCount, supportingMethods: methods,
                        score_available: false, score_calibrated: false,
                        score_model: 'unscored_family_calls/1', rawFields: new LinkedHashMap(raw),
                        // The shared claims writer stringifies nested objects. Preserve explicit
                        // JSON copies so intervals and subfamily calls remain machine-readable.
                        rawFieldsJson: new ObjectMapper().writeValueAsString(raw),
                        supportingMethodsJson: new ObjectMapper().writeValueAsString(methods)])
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

    private void validateInterface() {
        if (!(interfaceVersion in ['4', '5'])) throw new IllegalArgumentException('dbCAN interfaceVersion must be 4 or 5')
    }
}
