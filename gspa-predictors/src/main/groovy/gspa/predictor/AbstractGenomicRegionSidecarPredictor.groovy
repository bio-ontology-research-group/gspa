package gspa.predictor

import gspa.integration.EvidenceType
import gspa.io.FastaWriter
import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein

/**
 * Base class for predictors that delegate to
 * {@code benchmark/neural/run_genomic_predictors.py}.
 *
 * <p>Output is the 6-column genomic-region TSV:
 * {@code contig_id<TAB>region_start<TAB>region_end<TAB>region_type<TAB>score<TAB>attributes}
 * with 1-based inclusive contig coordinates.
 *
 * <p>Each row becomes one {@link Annotation} with
 * {@link Annotation#contigId}, {@link Annotation#genomicStart},
 * {@link Annotation#genomicEnd}, and
 * {@link EvidenceType#GENOMIC_REGION_ML} (correlation group {@code viral}).
 *
 * <p>Subclasses provide:
 * <ul>
 *   <li>{@link #getPredictorName} — the {@code --predictor} flag value
 *   <li>{@link #getOutputAnnotationType} — JVM annotation type
 *       (typically {@code PROPHAGE}, {@code PLASMID}, or
 *       {@code VIRAL_CONTIG})
 *   <li>{@link #extraSidecarArgs} — per-tool flags
 * </ul>
 *
 * <p>Manifest contract differs from the protein sidecars: it carries a
 * {@code genome_fasta} column instead of {@code fasta_path} and an
 * optional {@code gff_path} for tools that need per-CDS features
 * (PhiSpy). The {@link #buildCommand} method writes the manifest
 * accordingly.
 */
abstract class AbstractGenomicRegionSidecarPredictor extends AbstractToolPredictor implements GenomePredictor {

    /** Absolute path to {@code run_genomic_predictors.py}. */
    String sidecarScript

    String pythonExecutable = 'python3'

    /** Drop predictions below this score. */
    double minScore = 0.5

    /** Synthetic tag used in the manifest; identifies the output filename. */
    String tag = 'query'

    /** Optional GFF path for tools that need per-CDS features (PhiSpy). */
    String gffPath

    /** Evidence type reported on emitted annotations. */
    EvidenceType evidenceType = EvidenceType.GENOMIC_REGION_ML

    /** The {@code --predictor} flag value for the sidecar. */
    abstract String getPredictorName()

    /** AnnotationType this predictor primarily emits (e.g. PROPHAGE). */
    abstract AnnotationType getOutputAnnotationType()

    /** Extra sidecar CLI flags (per-tool). */
    protected List<String> extraSidecarArgs() { [] }

    @Override
    String getExecutable() { pythonExecutable }

    /** Missing mandatory configuration (e.g. a model database path), or null when complete. */
    protected String configurationProblem() { null }

    @Override
    boolean isAvailable() {
        if (!sidecarScript || !(new File(sidecarScript).exists())) return false
        if (configurationProblem()) return false
        try {
            def proc = [pythonExecutable, '--version'].execute()
            proc.waitForOrKill(10000)
            return proc.exitValue() == 0
        } catch (Exception ignored) {
            return false
        }
    }

    @Override
    String unavailableReason() {
        if (!sidecarScript) return 'sidecarScript (run_genomic_predictors.py) is not configured'
        if (!new File(sidecarScript).exists()) return "sidecar script not found: ${sidecarScript}"
        configurationProblem() ?: super.unavailableReason()
    }

    /** Region predictors need contig sequences and coordinates: they run on a Genome. */
    @Override
    Map<String, List<Annotation>> predictBatch(List<Protein> proteins) {
        throw new UnsupportedOperationException(
            "${name} is a genome-level predictor; call predictGenome(Genome) with nucleotide contigs")
    }

    /** File name the nucleotide (or tool-specific) input is written under. */
    protected String inputFileName() { 'genome.fna' }

    /**
     * Write the tool input from the genome's contigs. Default: nucleotide FASTA of
     * every contig that carries a sequence (contig ids are the FASTA ids, so the
     * region TSV's contig column joins back to {@link Genome#findContig}).
     */
    protected File writeGenomeInput(Genome genome, File dir) {
        def f = new File(dir, inputFileName())
        f.withWriter { w ->
            genome.contigs.findAll { it.sequence }.each { Contig c ->
                FastaWriter.writeEntry(w, c.id, c.sequence)
            }
        }
        f
    }

    /** Extra metadata stamped on every emitted region (score semantics, coordinate scope). */
    protected Map<String, Object> regionMetadata() { [:] }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        lastRunStatus = null
        genome.metadata.remove("genomic_regions:${name}".toString())
        long t0 = System.currentTimeMillis()
        // The pipeline's placeholder contig holds proteins without a CDS placement; it is not
        // an unscanned sequence, but its proteins can never be attributed to a region.
        List<String> unplaced = genome.contigs.findAll { isPlaceholderContig(it) }.collectMany { it.proteins*.id } as List<String>
        List<String> noSequence = genome.contigs.findAll { !it.sequence && !isPlaceholderContig(it) }*.id
        if (noSequence.size() == genome.contigs.count { !isPlaceholderContig(it) }) {
            recordStatus(PredictorRunStatus.failed(name,
                'no contig carries nucleotide sequence (protein-only input cannot be scanned for genomic regions)'), t0)
            return [:]
        }
        def tmpDir = File.createTempDir("gspa_${name}_", '')
        try {
            File input = writeGenomeInput(genome, tmpDir)
            def result = runOrFail(buildCommand(input, tmpDir), tmpDir, t0)
            if (result == null) return [:]

            // The sidecar exits 0 even when the wrapped tool failed and then writes no TSV;
            // parseOutput turns that into a contract failure. A header-only TSV is the
            // sidecar's "tool ran, no region passed" result.
            Map<String, List<Annotation>> regions = parseGuarded(t0) { parseOutput(tmpDir) }

            def projection = projectRegions(genome, regions)
            // Full contig-level records stay available even where no CDS overlaps.
            genome.metadata["genomic_regions:${name}".toString()] = projection.validRegions
            int n = projection.validRegions.size()

            List<String> problems = []
            List<String> uncovered = []
            if (noSequence) {
                int real = genome.contigs.count { !isPlaceholderContig(it) }
                problems << "${noSequence.size()} of ${real} contigs have no nucleotide sequence and were not " +
                    "analysed: ${noSequence.take(10).join(', ')}${noSequence.size() > 10 ? ', ...' : ''}"
                uncovered.addAll(noSequence)
            }
            if (unplaced) {
                problems << "${unplaced.size()} proteins have no CDS placement (ambiguous or conflicting CDS identity) and " +
                    "cannot be attributed to genomic regions: ${unplaced.take(10).join(', ')}${unplaced.size() > 10 ? ', ...' : ''}"
                uncovered.addAll(unplaced)
            }
            if (projection.invalidRegions) {
                problems << "${projection.invalidRegions.size()} region row(s) rejected for invalid coordinates " +
                    "(outside the contig or end < start): ${projection.invalidRegions.take(5).join('; ')}"
                uncovered.addAll(projection.invalidRegionContigs)
            }
            if (projection.invalidCds) {
                problems << "${projection.invalidCds.size()} CDS have invalid or missing coordinates; region membership " +
                    "could not be determined for: ${projection.invalidCds.take(10).join(', ')}${projection.invalidCds.size() > 10 ? ', ...' : ''}"
                uncovered.addAll(projection.invalidCds)
            }
            if (problems) {
                recordStatus(PredictorRunStatus.partial(name, n, problems.join(' | '), uncovered.unique(), result.exitCode), t0)
            } else {
                recordStatus(PredictorRunStatus.completed(name, n, result.exitCode, successNote()), t0)
            }
            projection.byProtein
        } catch (Exception e) {
            if (lastRunStatus == null || lastRunStatus.outcome != PredictorOutcome.FAILED) {
                recordStatus(PredictorRunStatus.failed(name, "${e.class.simpleName}: ${e.message}"), t0)
            }
            throw e
        } finally {
            tmpDir.deleteDir()
        }
    }

    /** Outcome of validating regions and projecting them onto CDS. */
    static class RegionProjection {
        /** protein id -> region copies for CDS overlapping a valid region */
        Map<String, List<Annotation>> byProtein = [:]
        List<Annotation> validRegions = []
        /** "contig:start-end" descriptions of rejected region rows */
        List<String> invalidRegions = []
        /** contig ids carrying at least one rejected region row */
        Set<String> invalidRegionContigs = [] as Set
        /** ids of proteins whose CDS coordinates are missing/invalid on a contig that has regions */
        List<String> invalidCds = []
    }

    /**
     * Validate region coordinates against the contig (1-based inclusive, start >= 1,
     * start <= end <= contig length when the length is known, contig must exist) and
     * project valid regions onto the proteins whose CDS overlaps them. Proteins with
     * missing or impossible CDS coordinates are not guessed at; they are reported.
     */
    static RegionProjection projectRegions(Genome genome, Map<String, List<Annotation>> regions) {
        def out = new RegionProjection()
        regions.each { String contigId, List<Annotation> anns ->
            Contig contig = genome.findContig(contigId)
            Set<String> badCdsSeen = [] as Set
            anns.each { Annotation region ->
                int rs = region.genomicStart, re = region.genomicEnd
                boolean ok = contig != null && rs >= 1 && re >= rs && (contig.length == 0 || re <= contig.length)
                if (!ok) {
                    out.invalidRegions << "${contigId}:${rs}-${re}".toString()
                    out.invalidRegionContigs << contigId
                    return
                }
                out.validRegions << region
                contig.proteins.each { Protein p ->
                    int s = p.start, e = p.end
                    boolean cdsOk = p.sourceFeature != null && s >= 1 && e >= s && (contig.length == 0 || e <= contig.length)
                    if (!cdsOk) {
                        if (badCdsSeen.add(p.id)) out.invalidCds << p.id
                        return
                    }
                    if (e < rs || s > re) return
                    boolean contained = s >= rs && e <= re
                    Map<String, Object> meta = new LinkedHashMap<>(region.metadata ?: [:])
                    meta.cdsOverlap = contained ? 'contained' : 'partial'
                    (out.byProtein.computeIfAbsent(p.id) { [] } as List<Annotation>) << new Annotation(
                        type: region.type, value: region.value, score: region.score,
                        source: region.source, evidence: region.evidence,
                        evidenceType: region.evidenceType, contigId: region.contigId,
                        genomicStart: rs, genomicEnd: re,
                        regionType: region.regionType, metadata: meta)
                }
            }
        }
        out
    }

    @Override
    Set<AnnotationType> getOutputTypes() { [outputAnnotationType] as Set }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        if (!sidecarScript) {
            throw new IllegalStateException(
                "${name}: sidecarScript is unset. Point it at run_genomic_predictors.py.")
        }
        def manifest = new File(outputDir, 'manifest.tsv')
        def gff = gffPath ?: '-'
        manifest.text = "tag\tgenome_fasta\tgff_path\toutput_dir\n" +
                "${tag}\t${inputFasta.absolutePath}\t${gff}\t${outputDir.absolutePath}\n"
        def cmd = [
            pythonExecutable,
            new File(sidecarScript).absolutePath,
            '--predictor', predictorName,
            '--manifest', manifest.absolutePath,
            '--min-score', minScore.toString(),
        ] as List<String>
        cmd.addAll(extraSidecarArgs())
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // Sidecar contract (run_genomic_predictors.py): the 6-column TSV is opened, header first,
        // only after the wrapped tool produced results; on tool failure the sidecar logs a warning,
        // writes NOTHING and still exits 0. So: no TSV = failed run; header-only = no region passed.
        // The bundled geNomad runner requires both validated native summaries before publishing.
        // PhiSpy still treats a missing native coordinate file as zero (see successNote()).
        def tsv = requireOutput(outputDir, "${tag}.${predictorName}.genomic.tsv",
            'run_genomic_predictors.py writes a header row whenever the wrapped tool produced output')
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        boolean firstLine = true
        int lineNo = 0
        tsv.eachLine { line ->
            lineNo++
            if (firstLine) {
                firstLine = false
                if (line != 'contig_id\tregion_start\tregion_end\tregion_type\tscore\tattributes') {
                    throw new ToolOutputContractException("${name} ${tsv.name}: invalid genomic TSV header")
                }
                return
            }
            if (!line || line.startsWith('#')) return
            def parts = line.split('\t', -1)
            if (parts.size() != 6) {
                throw new ToolOutputContractException("${name} ${tsv.name} line ${lineNo}: ${parts.size()} columns, expected 6")
            }
            try {
                String contig = parts[0]
                int start = parts[1] as int
                int end = parts[2] as int
                String regionType = parts[3]
                double score = parts[4] as double
                String attributesRaw = parts[5]
                if (!Double.isFinite(score) || score < 0.0d || score > 1.0d) {
                    throw new ToolOutputContractException("${name} ${tsv.name} line ${lineNo}: nonfinite or out-of-range score for ${contig}")
                }
                if (!(regionType in ['prophage', 'plasmid', 'viral_contig'])) {
                    throw new ToolOutputContractException("${name} ${tsv.name} line ${lineNo}: unknown region type '${regionType}'")
                }
                // Resolve annotation type from region_type when the predictor
                // emits multiple kinds (e.g. geNomad: prophage / plasmid /
                // viral_contig).
                AnnotationType t = outputAnnotationType
                switch (regionType) {
                    case 'prophage':
                        t = AnnotationType.PROPHAGE
                        break
                    case 'plasmid':
                        t = AnnotationType.PLASMID
                        break
                    case 'viral_contig':
                        t = AnnotationType.VIRAL_CONTIG
                        break
                }
                Map<String, Object> meta = [:]
                if (attributesRaw) {
                    attributesRaw.split('\\|', -1).each { kv ->
                        def i = kv.indexOf('=')
                        if (i <= 0 || meta.containsKey(kv.substring(0, i))) {
                            throw new ToolOutputContractException("${name} ${tsv.name} line ${lineNo}: malformed or duplicate attribute '${kv}'")
                        }
                        meta[kv.substring(0, i)] = kv.substring(i + 1)
                    }
                }
                // Validate complete rows even when they will be filtered out.
                if (score < minScore) return
                meta.putAll(regionMetadata())
                // Genomic-region annotations are keyed by contig (not protein);
                // we use the contig id as the map key so downstream consumers
                // can group / flatten by chromosomal location.
                results[contig] << new Annotation(
                    type: t,
                    value: regionType,
                    source: name,
                    score: score,
                    evidence: 'IEA',
                    evidenceType: evidenceType,
                    contigId: contig,
                    genomicStart: start,
                    genomicEnd: end,
                    regionType: regionType,
                    metadata: meta,
                )
            } catch (NumberFormatException e) {
                throw new ToolOutputContractException("${name} ${tsv.name} line ${lineNo}: non-numeric coordinate/score (${e.message})")
            }
        }
        results
    }
}
