package gspa.predictor.pathway

import gspa.model.*
import gspa.integration.EvidenceType
import groovy.json.JsonOutput
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor
import gspa.predictor.PredictorOutcome
import gspa.predictor.PredictorRunStatus
import gspa.predictor.ToolOutputContractException

/**
 * Protein reaction-candidate search wrapping gapsmith or gapseq.
 * Named-column alignment rows are sequence evidence, not independent pathway support.
 * Weak hits are retained; complex membership is represented separately from EC activity.
 * Original rows remain available for context/coverage reporting. Nucleotide gapseq
 * hits cannot be projected to proteins without a separate coordinate attribution step.
 */
class GapseqPredictor extends AbstractToolPredictor implements GenomePredictor {

    /**
     * External tool to invoke. {@code 'gapsmith'} (default) uses
     * {@code gapsmith find -p all -t Bacteria -o <out-dir>}.
     * {@code 'gapseq'} uses {@code gapseq find -p all -t Bacteria}
     * (legacy R/bash upstream).
     */
    String tool = 'gapsmith'

    /**
     * Reference-data directory for gapsmith ({@code --data-dir}).
     * When null, gapsmith auto-detects (assumes a sibling {@code data/}
     * dir, or {@code GAPSMITH_DATA_DIR}). Ignored for {@code tool=='gapseq'}.
     */
    String dataDir = null

    /** gapseq/gapsmith medium file for gap-filling */
    String medium

    /** Export SBML model */
    boolean exportSbml = false

    /** Taxonomy for pathway prediction */
    String taxonomy = 'Bacteria'

    /**
     * Phase 10: search target for gapseq's reaction-sequence database.
     * <ul>
     *   <li>{@code genome}  — tblastn against genome nucleotide (gapseq only)</li>
     *   <li>{@code proteome} — blastp against the full predicted proteome (default)
     *       (~5× faster but loses small-ORF / frameshift rescue)</li>
     *   <li>{@code reps}    — blastp against cluster representatives only
     *       (requires intra-genome clustering to have run). A reaction is
     *       "found" if any cluster rep has a qualifying hit — equivalent
     *       to querying the whole proteome under the 90% identity
     *       equivalence, and the form Phase 11 will cache across genomes.</li>
     * </ul>
     */
    String target = 'proteome'

    /** Bound external alignment parallelism to the caller's budget. */
    int threads = 1

    /**
     * When {@code target=='reps'}, the list of representative protein IDs
     * (from the intra-genome clusterer). {@code null} with {@code target=='reps'}
     * is an error.
     */
    Set<String> representativeIds = null

    @Override
    String getName() { tool }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.EC, AnnotationType.CUSTOM] as Set
    }

    @Override
    String getExecutable() {
        tool == 'gapseq' ? 'gapseq' : 'gapsmith'
    }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        validateInputConfiguration()
        def exe = executablePath ?: getExecutable()
        if (tool == 'gapseq') {
            return [
                exe,
                'find',
                '-p', 'all',
                '-t', taxonomy,
                '-M', target == 'genome' ? 'nucl' : 'prot',
                '-K', threads.toString(),
                '-f', outputDir.absolutePath,
                inputFasta.absolutePath,
            ]
        }
        // gapsmith: --data-dir is at the top level, out-dir goes to the subcommand
        def cmd = [exe]
        if (dataDir) cmd += ['--data-dir', dataDir]
        cmd += [
            'find',
            '-p', 'all',
            '-t', taxonomy,
            '-K', threads.toString(),
            '-o', outputDir.absolutePath,
            inputFasta.absolutePath,
        ]
        cmd
    }

    /** Original named-column rows, including no-hit and unknown-reference records. */
    List<Map<String, String>> reactionObservations = []
    String reactionTableSha256

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        reactionObservations = []
        reactionTableSha256 = null
        List<File> files = outputDir.listFiles()?.findAll { it.name.endsWith('-Reactions.tbl') && it.isFile() } ?: []
        if (files.size() != 1) throw new ToolOutputContractException("${name}: expected exactly one *-Reactions.tbl, found ${files.size()}")
        File table = files.first()
        if (!table.length()) throw new ToolOutputContractException("${name}: Reactions.tbl is zero bytes")
        List<String> header = null
        List<Map<String, String>> rows = []
        List<String> required = ['pathway','rxn','ec','keyrea','spont','is_complex','qseqid',
                                 'pident','evalue','bitscore','qcovs','stitle','sstart','send','status']
        table.eachLine('UTF-8') { String line, int number ->
            if (!line.trim() || line.startsWith('#')) return
            List<String> fields = line.split('\t', -1).toList()
            if (header == null) {
                header = fields
                if (header.toSet().size() != header.size() || !header.containsAll(required)) {
                    throw new ToolOutputContractException("${name}: unsupported reaction header on line ${number}")
                }
                return
            }
            if (fields.size() != header.size()) throw new ToolOutputContractException("${name}: malformed reaction row ${number}")
            Map<String, String> row = [:]
            header.eachWithIndex { key, i -> row[key] = fields[i] }
            if (!row.pathway.trim() || !row.rxn.trim()) throw new ToolOutputContractException("${name}: empty pathway/reaction identifier")
            for (String key : ['keyrea','spont','is_complex']) flag(row[key])
            if (!(row.status in ['good_blast','bad_blast','no_blast','no_seq_data','spontaneous'])) {
                throw new ToolOutputContractException("${name}: unknown reaction status '${row.status}'")
            }
            if (row.status in ['good_blast','bad_blast']) {
                if (!row.qseqid.trim() || !row.stitle.trim()) throw new ToolOutputContractException("${name}: hit lacks reference or target identity")
                statistic(row, 'pident', 100d); statistic(row, 'qcovs', 100d)
                statistic(row, 'bitscore'); statistic(row, 'evalue')
                for (String key : ['sstart','send']) {
                    if (!(row[key] ==~ /[1-9][0-9]*/)) throw new ToolOutputContractException("${name}: invalid hit coordinate ${key}")
                }
            }
            ecNumbers(row.ec)
            rows << row
        }
        if (header == null) throw new ToolOutputContractException("${name}: missing reaction header")
        Map<List<String>, List<Map<String, String>>> groups = [:]
        rows.findAll { it.status in ['good_blast','bad_blast'] }.each { row ->
            String pid = row.stitle.trim().split(/\s+/, 2)[0]
            List<String> ecs = ecNumbers(row.ec)
            boolean complex = flag(row.is_complex)
            AnnotationType type = complex || !ecs ? AnnotationType.CUSTOM : AnnotationType.EC
            String localReaction = [row.pathway, row.rxn].collect {
                java.net.URLEncoder.encode(it.trim().replaceAll(/^\|+|\|+$/, ''), 'UTF-8')
            }.join(':')
            List<String> values = complex ? ["reaction-subunit:${localReaction}".toString()] :
                (ecs ? ecs.collect { "EC:${it}".toString() } : ["reaction:${localReaction}".toString()])
            values.each { value -> groups.computeIfAbsent([pid,type.name(),value], { [] }) << row }
        }
        Map<String, List<Annotation>> result = [:].withDefault { [] }
        groups.each { key, hits ->
            def best = hits.max { statistic(it, 'pident', 100d) }
            String role = key[2].startsWith('reaction-subunit:') ? 'COMPLEX_SUBUNIT' :
                (key[1] == 'EC' ? 'CATALYTIC' : 'UNSPECIFIED')
            result[key[0]] << new Annotation(type: AnnotationType.valueOf(key[1]), value: key[2],
                source: name, evidenceType: EvidenceType.SEQUENCE_SIMILARITY,
                score: statistic(best, 'pident', 100d) / 100d,
                metadata: [pident:statistic(best,'pident',100d), bitscore:statistic(best,'bitscore'),
                    evalue:statistic(best,'evalue'), qcovs:statistic(best,'qcovs',100d),
                    referenceId:best.qseqid, reaction:best.rxn, candidate_role:role,
                    score_semantics:'alignment_percent_identity_div_100_not_probability',
                    hitsJson:JsonOutput.toJson(hits)])
        }
        reactionObservations = rows
        reactionTableSha256 = java.security.MessageDigest.getInstance('SHA-256').digest(table.bytes).encodeHex().toString()
        result
    }

    protected static boolean flag(String value) {
        if (!(value in ['TRUE','FALSE'])) throw new ToolOutputContractException("invalid reaction boolean '${value}'")
        value == 'TRUE'
    }

    protected static double statistic(Map row, String key, double maximum = Double.MAX_VALUE) {
        double value
        try { value = Double.parseDouble(row[key]) }
        catch (Exception e) { throw new ToolOutputContractException("invalid reaction statistic ${key}: ${row[key]}") }
        if (!Double.isFinite(value) || value < 0d || value > maximum) {
            throw new ToolOutputContractException("invalid reaction statistic ${key}: ${row[key]}")
        }
        value
    }

    protected static List<String> ecNumbers(String value) {
        List<String> numbers = []
        for (String token : value.trim().split('[/,;\\s]+')) {
            String ec = token.replaceFirst(/^EC:/, '')
            if (!ec || ec in ['NA','-']) continue
            // Local reference identifiers (e.g. 1.5.1.M4) are not asserted as ECs.
            // The original token remains in rows/hitsJson, and a reaction candidate
            // is retained when no standard or partial EC is available.
            if (ec ==~ /[0-9]+\.[0-9]+\.[0-9]+\.(?!n[0-9]+$)[A-Za-z][A-Za-z0-9]*/) continue
            if (!(ec ==~ /[0-9]+(?:\.[0-9]+){0,3}/) &&
                    !(ec ==~ /[0-9]+\.(?:[0-9]+|-)\.(?:[0-9]+|-)\.(?:[0-9]+|n[0-9]+|-)/)) {
                throw new ToolOutputContractException("invalid reaction EC '${token}'")
            }
            if (!numbers.contains(ec)) numbers << ec
        }
        numbers
    }

    private void validateInputConfiguration() {
        if (!(tool in ['gapsmith', 'gapseq'])) {
            throw new IllegalArgumentException("Unknown pathway tool: ${tool}; expected gapsmith|gapseq")
        }
        if (!(target in ['genome', 'proteome', 'reps'])) {
            throw new IllegalArgumentException("Unknown pathway target: ${target}; expected genome|proteome|reps")
        }
        if (tool == 'gapsmith' && target == 'genome') {
            throw new IllegalArgumentException('gapsmith find requires protein input; select target=proteome or reps')
        }
        if (threads < 1) throw new IllegalArgumentException('pathway threads must be positive')
    }

    @Override
    Map<String, List<Annotation>> predictBatch(List<Protein> proteins) {
        reactionObservations = []
        reactionTableSha256 = null
        long started = System.currentTimeMillis()
        if (target == 'genome') {
            recordStatus(PredictorRunStatus.failed(name, 'predictBatch receives proteins; target=genome requires predictGenome'), started)
            throw new IllegalArgumentException(lastRunStatus.detail)
        }
        def result = super.predictBatch(proteins)
        if (lastRunStatus?.outcome != PredictorOutcome.FAILED &&
                !(proteins.collect { it.id } as Set).containsAll(result.keySet())) {
            recordStatus(PredictorRunStatus.failed(name,
                'pathway output names proteins absent from the submitted search target'), started)
            throw new ToolOutputContractException(lastRunStatus.detail)
        }
        result
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        lastRunStatus = null
        reactionObservations = []
        reactionTableSha256 = null
        def tmpDir = File.createTempDir("gspa_gapseq_", '')
        long t0 = System.currentTimeMillis()
        try {
            validateInputConfiguration()
            if (target == 'reps' && (representativeIds == null || representativeIds.isEmpty())) {
                throw new IllegalStateException(
                    "GapseqPredictor(target='reps') requires representativeIds to be populated " +
                    "by an earlier intra-genome clustering step")
            }
            File inputFasta
            ContigExport export = null
            List<String> skippedProteins = []
            if (target == 'genome') {
                inputFasta = new File(tmpDir, 'genome.fna')
                export = writeContigFasta(genome, inputFasta)
                if (export.written.isEmpty()) {
                    recordStatus(PredictorRunStatus.failed(name,
                        'target=genome but no contig carries nucleotide sequence'), t0)
                    return [:]
                }
            } else {
                inputFasta = writeTargetFasta(genome, tmpDir)
                skippedProteins = genome.proteins.findAll { p ->
                    (target != 'reps' || representativeIds.contains(p.id)) && !p.sequence
                }.collect { it.id }
                if (inputFasta.length() == 0) {
                    recordStatus(PredictorRunStatus.failed(name,
                        "target=${target} selected no protein sequences to search"), t0)
                    return [:]
                }
            }

            def result = runOrFail(buildCommand(inputFasta, tmpDir), tmpDir, t0)
            if (result == null) return [:]

            def parsed = parseGuarded(t0) {
                def annotations = parseOutput(tmpDir)
                if (target == 'genome' && annotations) {
                    throw new ToolOutputContractException('nucleotide pathway hits require coordinate-to-protein attribution; use target=proteome for protein claims')
                }
                Set<String> searched = genome.proteins.findAll { p -> p.sequence &&
                    (target != 'reps' || representativeIds.contains(p.id)) }.collect { it.id } as Set
                if (target != 'genome' && !searched.containsAll(annotations.keySet())) {
                    throw new ToolOutputContractException('pathway output names proteins absent from the submitted search target')
                }
                annotations
            }
            int n = countAnnotations(parsed)
            if (export?.skipped) {
                recordStatus(PredictorRunStatus.partial(name, n, describeSkipped(export, genome.contigs.size()),
                    export.skipped, result.exitCode), t0)
            } else if (skippedProteins) {
                recordStatus(PredictorRunStatus.partial(name, n,
                    "${skippedProteins.size()} selected proteins lack sequences", skippedProteins, result.exitCode), t0)
            } else {
                recordStatus(PredictorRunStatus.completed(name, n, result.exitCode), t0)
            }
            return parsed
        } catch (Exception e) {
            if (lastRunStatus == null || lastRunStatus.outcome != PredictorOutcome.FAILED) {
                recordStatus(PredictorRunStatus.failed(name, "${e.class.simpleName}: ${e.message}"), t0)
            }
            throw e
        } finally {
            tmpDir.deleteDir()
        }
    }

    /** Emit the FASTA gapseq will search, based on the configured target. */
    private File writeTargetFasta(Genome genome, File dir) {
        switch (target) {
            case 'genome':
                def f = new File(dir, 'genome.fna')
                f.withWriter { w ->
                    genome.contigs.each { contig ->
                        if (contig.sequence) {
                            w.writeLine(">${contig.id}")
                            w.writeLine(contig.sequence)
                        }
                    }
                }
                return f
            case 'proteome':
                def f = new File(dir, 'proteins.faa')
                f.withWriter { w ->
                    genome.proteins.each { p ->
                        if (p.sequence) { w.writeLine(">${p.id}"); w.writeLine(p.sequence) }
                    }
                }
                return f
            case 'reps':
                def f = new File(dir, 'reps.faa')
                f.withWriter { w ->
                    genome.proteins.each { p ->
                        if (p.sequence && representativeIds.contains(p.id)) {
                            w.writeLine(">${p.id}"); w.writeLine(p.sequence)
                        }
                    }
                }
                return f
            default:
                throw new IllegalArgumentException("Unknown gapseq target: ${target}; expected genome|proteome|reps")
        }
    }

    /**
     * Run gapseq in community mode for crossfeeding analysis.
     * Takes multiple genome FASTA files and finds metabolic exchanges.
     */
    Map<String, Object> analyzeCommmunity(List<File> genomeFastas, File outputDir) {
        // gapseq community analysis would be done in a separate step
        // after individual models are built
        log.info("Community analysis with ${genomeFastas.size()} genomes")
        [:]
    }

}
