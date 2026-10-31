package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.model.Protein
import gspa.predictor.AbstractGenomicRegionSidecarPredictor
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor
import gspa.predictor.PredictorOutcome
import gspa.predictor.PredictorRunStatus
import gspa.predictor.ToolOutputContractException
import java.security.MessageDigest

/** antiSMASH regions are unscored genomic observations, not protein enzyme assignments. */
class AntiSmashPredictor extends AbstractToolPredictor implements GenomePredictor {
    int threads = Runtime.runtime.availableProcessors()
    String strictness = 'relaxed'
    String database
    /** Record coverage retained independently of whether any region was found. */
    List<Map<String, Object>> recordObservations = []

    @Override
    String getName() { 'antismash' }
    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.SECONDARY_METABOLITE] as Set }
    @Override
    String getExecutable() { 'antismash' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (!(strictness in ['strict', 'relaxed', 'loose']) || threads < 1) {
            throw new IllegalArgumentException('antismash: invalid strictness or thread count')
        }
        List<String> command = [executablePath ?: executable, inputFasta.absolutePath,
            '--output-dir', outputDir.absolutePath, '--cpus', threads.toString(),
            '--hmmdetection-strictness', strictness, '--minimal', '--genefinding-tool', 'prodigal']
        if (database) command.addAll(['--databases', database])
        command
    }

    protected static String digest(String sequence) {
        MessageDigest.getInstance('SHA-256').digest(sequence.toUpperCase(Locale.ROOT).getBytes('UTF-8')).encodeHex().toString()
    }

    protected static int coordinate(Object value, String label) {
        if (!(value instanceof Number) || !Double.isFinite(value.doubleValue()) ||
                value.doubleValue() != value.longValue() || value.longValue() < 0 || value.longValue() > Integer.MAX_VALUE) {
            throw new ToolOutputContractException("antismash: invalid ${label}: ${value}")
        }
        value.intValue()
    }

    /** Parse native schema-4 areas into contig coordinates; CDS locus tags are never joined. */
    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        recordObservations = []
        List<Map> candidates = []
        def files = outputDir.listFiles()?.findAll { it.isFile() && it.name.endsWith('.json') } ?: []
        files.sort { it.name }.each { file ->
            Object data
            try { data = new groovy.json.JsonSlurper().parse(file) }
            catch (Exception e) { throw new ToolOutputContractException("antismash: cannot parse ${file.name}: ${e.message}") }
            if (data instanceof Map && (data.containsKey('records') || data.containsKey('schema'))) {
                candidates << [file: file, data: data]
            }
        }
        if (candidates.size() != 1) {
            throw new ToolOutputContractException('antismash: expected exactly one native results JSON with records/schema')
        }
        parseNative(candidates[0].data as Map, candidates[0].file as File)
    }

    /** Import one explicitly selected result, without executing antiSMASH. */
    Map<String, List<Annotation>> parseResultFile(File resultFile) {
        recordObservations = []
        Object data = new groovy.json.JsonSlurper().parse(resultFile)
        if (!(data instanceof Map)) throw new ToolOutputContractException('antismash: native result must be an object')
        parseNative(data as Map, resultFile)
    }

    protected Map<String, List<Annotation>> parseNative(Map data, File resultFile) {
        if (data.schema != 4 || !(data.version instanceof String) || !data.version.trim() || !(data.records instanceof List)) {
            throw new ToolOutputContractException('antismash: unsupported results contract; expected schema 4, version and records list')
        }
        String sourceHash = MessageDigest.getInstance('SHA-256').digest(resultFile.bytes).encodeHex().toString()
        Map<String, List<Annotation>> results = [:]
        Set<String> seen = [] as Set
        data.records.each { record ->
            if (!(record instanceof Map) || !(record.id instanceof String) || !record.id ||
                    !(record.seq instanceof Map) || !(record.seq.data instanceof String) || !record.seq.data ||
                    !(record.areas instanceof List) || !(record.modules instanceof Map) || !(record.annotations instanceof Map)) {
                throw new ToolOutputContractException('antismash: record lacks native identity, sequence, areas or modules')
            }
            if (record.containsKey('original_id') && (!(record.original_id instanceof String) || !record.original_id)) {
                throw new ToolOutputContractException('antismash: malformed original record identity')
            }
            String id = record.original_id ?: record.id
            if (!seen.add(id)) throw new ToolOutputContractException("antismash: duplicate record ${id}")
            int length = record.seq.data.length()
            Object detection = record.modules['antismash.detection.hmm_detection']
            boolean analysed = detection instanceof Map && detection.record_id == record.id &&
                detection.rule_results instanceof Map && detection.enabled_types instanceof List
            if (detection != null && !analysed) {
                throw new ToolOutputContractException("antismash: malformed detection results for ${id}")
            }
            recordObservations << [contig_id: id, record_id: record.id, length: length,
                sequence_sha256: digest(record.seq.data), detection_available: analysed,
                region_count: record.areas.size(), source_version: data.version, source_sha256: sourceHash]
            results[id] = []
            record.areas.eachWithIndex { area, index ->
                if (!(area instanceof Map) || !(area.products instanceof List) ||
                        !area.products.every { it instanceof String && it.trim() }) {
                    throw new ToolOutputContractException("antismash: malformed region products for ${id}")
                }
                int start = coordinate(area.start, 'region start'), end = coordinate(area.end, 'region end')
                if (start >= length || end > length || start == end ||
                        (start > end && record.annotations.topology != 'circular')) {
                    throw new ToolOutputContractException("antismash: invalid region bounds for ${id}: ${start}-${end}")
                }
                if (!analysed) throw new ToolOutputContractException("antismash: region without detection results for ${id}")
                List<List<Integer>> segments = start < end ? [[start + 1, end]] : [[start + 1, length]]
                if (start > end && end > 0) segments << [1, end]
                segments.eachWithIndex { segment, part ->
                    results[id] << new Annotation(type: AnnotationType.SECONDARY_METABOLITE,
                        value: area.products ? area.products.join(', ') : 'unclassified region', source: name,
                        contigId: id, genomicStart: segment[0], genomicEnd: segment[1], regionType: 'biosynthetic_region',
                        metadata: [score_available: false, observation_type: 'biosynthetic_region_membership',
                            products: new ArrayList(area.products), region_id: "${id}:region${index + 1}".toString(),
                            segment_index: part, region_segments: segments, source_version: data.version,
                            source_schema: data.schema, source_sha256: sourceHash, record_id: record.id,
                            native_start: start, native_end: end,
                            interpretation: 'overlap with a predicted biosynthetic region; not enzyme activity, expression or a probability'])
                }
            }
        }
        results
    }

    /** The same sequence validation and coordinate projection serve CLI imports and live runs. */
    Map importGenomeResults(Genome genome, File resultFile) {
        genome.metadata.remove("genomic_regions:${name}".toString())
        genome.metadata.remove('antismash_record_coverage')
        projectNative(genome, parseResultFile(resultFile))
    }

    protected Map projectNative(Genome genome, Map regions) {
        Set<String> covered = [] as Set
        recordObservations.each { observation ->
            def contig = genome.findContig(observation.contig_id as String)
            if (contig == null || !contig.sequence || digest(contig.sequence) != observation.sequence_sha256) {
                throw new ToolOutputContractException("antismash: output record does not match input sequence: ${observation.contig_id}")
            }
            if (observation.detection_available) covered.add(contig.id)
        }
        def projection = AbstractGenomicRegionSidecarPredictor.projectRegions(genome, regions)
        List<String> skipped = (genome.contigs.findAll { !covered.contains(it.id) }*.id +
            projection.invalidCds + projection.invalidRegionContigs).unique()
        genome.metadata["genomic_regions:${name}".toString()] = projection.validRegions
        genome.metadata.antismash_record_coverage = recordObservations.collect { new LinkedHashMap(it) }
        [projection: projection, skipped: skipped]
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        lastRunStatus = null
        recordObservations = []
        String key = "genomic_regions:${name}".toString()
        genome.metadata.remove(key)
        genome.metadata.remove('antismash_record_coverage')
        long t0 = System.currentTimeMillis()
        if (workDir) workDir.mkdirs()
        File tmp = workDir ? java.nio.file.Files.createTempDirectory(workDir.toPath(), 'gspa_antismash_').toFile() :
            java.nio.file.Files.createTempDirectory('gspa_antismash_').toFile()
        try {
            File input = new File(tmp, 'genome.fna')
            def exported = writeContigFasta(genome, input)
            if (!exported.written) {
                recordStatus(PredictorRunStatus.failed(name, 'no contig carries nucleotide sequence'), t0)
                return [:]
            }
            def out = new File(tmp, 'output')
            def execution = runOrFail(buildCommand(input, out), tmp, t0)
            if (execution == null) return [:]
            Map regions = parseGuarded(t0) { parseOutput(out) }
            def imported = projectNative(genome, regions)
            def projection = imported.projection
            List<String> skipped = imported.skipped
            int count = projection.validRegions.collect { it.metadata.region_id }.unique().size()
            int overlaps = countAnnotations(projection.byProtein)
            String note = "${count} predicted regions; ${overlaps} protein/region-segment overlaps; membership is unscored and does not assert enzyme function"
            if (skipped) {
                recordStatus(PredictorRunStatus.partial(name, count,
                    note + '; missing detection coverage or invalid CDS coordinates: ' + skipped.take(10).join(', '), skipped, execution.exitCode), t0)
            } else {
                recordStatus(PredictorRunStatus.completed(name, count, execution.exitCode, note), t0)
            }
            projection.byProtein
        } catch (Exception e) {
            genome.metadata.remove(key)
            genome.metadata.remove('antismash_record_coverage')
            if (lastRunStatus?.outcome != PredictorOutcome.FAILED) {
                recordStatus(PredictorRunStatus.failed(name, "${e.class.simpleName}: ${e.message}"), t0)
            }
            throw e
        } finally { tmp.deleteDir() }
    }

    @Override
    Map<String, List<Annotation>> predictBatch(List<Protein> proteins) {
        throw new UnsupportedOperationException('antismash is a genome-level predictor; call predictGenome(Genome)')
    }
}
