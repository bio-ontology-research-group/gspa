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

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        long t0 = System.currentTimeMillis()
        if (!genome.contigs.any { it.sequence }) {
            recordStatus(PredictorRunStatus.failed(name,
                'no contig carries nucleotide sequence (protein-only input cannot be scanned for genomic regions)'), t0)
            return [:]
        }
        def tmpDir = File.createTempDir("gspa_${name}_", '')
        try {
            File input = writeGenomeInput(genome, tmpDir)
            def command = buildCommand(input, tmpDir)
            log.info("Running ${name}: ${command.join(' ')}")
            def result = execute(command, tmpDir)
            if (result.exitCode != 0) {
                log.error("${name} failed with exit code ${result.exitCode}: ${result.stderr}")
                recordStatus(PredictorRunStatus.failed(name,
                    "exit code ${result.exitCode}: ${result.stderr?.trim()?.take(500)}", result.exitCode), t0)
                return [:]
            }
            // The sidecar exits 0 even when the wrapped tool failed; it then writes
            // no TSV. A header-only TSV is the genuine "no regions" result.
            def tsv = new File(tmpDir, "${tag}.${predictorName}.genomic.tsv")
            if (!tsv.exists()) {
                recordStatus(PredictorRunStatus.failed(name,
                    "sidecar wrote no ${tsv.name}; wrapped tool failed (see sidecar log)", result.exitCode), t0)
                return [:]
            }
            Map<String, List<Annotation>> regions = parseOutput(tmpDir)
            // Full contig-level records stay available even where no CDS overlaps.
            genome.metadata["genomic_regions:${name}".toString()] = regions.values().flatten()
            Map<String, List<Annotation>> byProtein = mapRegionsToProteins(genome, regions)
            recordStatus(PredictorRunStatus.completed(name, countAnnotations(regions), result.exitCode), t0)
            byProtein
        } catch (Exception e) {
            recordStatus(PredictorRunStatus.failed(name, "${e.class.simpleName}: ${e.message}"), t0)
            throw e
        } finally {
            tmpDir.deleteDir()
        }
    }

    /**
     * Project contig-keyed regions onto the proteins whose CDS overlaps them
     * (same contig, 1-based inclusive coordinates). Each protein receives a copy
     * of the region annotation tagged with how much of the CDS lies in the region.
     */
    static Map<String, List<Annotation>> mapRegionsToProteins(Genome genome, Map<String, List<Annotation>> regions) {
        Map<String, List<Annotation>> out = [:]
        regions.each { String contigId, List<Annotation> anns ->
            Contig contig = genome.findContig(contigId)
            if (contig == null) return
            anns.each { Annotation region ->
                contig.proteins.each { Protein p ->
                    int s = p.start, e = p.end
                    if (e == 0 || e < region.genomicStart || s > region.genomicEnd) return
                    boolean contained = s >= region.genomicStart && e <= region.genomicEnd
                    Map<String, Object> meta = new LinkedHashMap<>(region.metadata ?: [:])
                    meta.cdsOverlap = contained ? 'contained' : 'partial'
                    (out.computeIfAbsent(p.id) { [] } as List<Annotation>) << new Annotation(
                        type: region.type, value: region.value, score: region.score,
                        source: region.source, evidence: region.evidence,
                        evidenceType: region.evidenceType, contigId: region.contigId,
                        genomicStart: region.genomicStart, genomicEnd: region.genomicEnd,
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
        def tsv = new File(outputDir, "${tag}.${predictorName}.genomic.tsv")
        if (!tsv.exists()) return [:]
        Map<String, List<Annotation>> results = [:].withDefault { [] }
        boolean firstLine = true
        tsv.eachLine { line ->
            if (firstLine) { firstLine = false; return }
            if (!line || line.startsWith('#')) return
            def parts = line.split('\t')
            if (parts.size() < 5) return
            try {
                String contig = parts[0]
                int start = parts[1] as int
                int end = parts[2] as int
                String regionType = parts[3]
                double score = parts[4] as double
                String attributesRaw = parts.size() >= 6 ? parts[5] : ''
                if (Double.isNaN(score) || score < minScore) return
                // Resolve annotation type from region_type when the predictor
                // emits multiple kinds (e.g. geNomad: prophage / plasmid /
                // viral_contig); fall back to the predictor's primary type.
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
                    attributesRaw.tokenize('|').each { kv ->
                        def i = kv.indexOf('=')
                        if (i > 0) meta[kv.substring(0, i)] = kv.substring(i + 1)
                    }
                }
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
            } catch (NumberFormatException ignored) {
                // skip malformed row
            }
        }
        results
    }
}
