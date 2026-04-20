# GSPA — a walkthrough

*2026-04-20T06:52:55Z by Showboat 0.6.1*
<!-- showboat-id: 56033b5b-1f96-415d-9490-70fcb88592de -->

GSPA (Genome-Scale Protein Annotation) is a multi-module Groovy/Java pipeline that turns raw protein FASTA plus predictor outputs into calibrated per-protein GO/EC posteriors, with quality metrics and suggestions for unannotated 'dark matter' proteins. This walkthrough traces the code path of `gspa integrate`, the benchmark-friendly entry point that consumes pre-parsed claims and produces a TSV of posteriors. We'll start at the module layout, descend to the data model, and climb back up through the integrator, the priors, and the dark-matter / reaction-local suggesters.

```bash
cat settings.gradle.kts && echo --- && ls gspa-core/src/main/groovy/gspa/
```

```output
rootProject.name = "gspa"

include("gspa-core")
include("gspa-predictors")
include("gspa-cli")
---
config
integration
io
metrics
model
ontology
```

Three Gradle modules. `gspa-core` holds the data model (proteins, annotations), the ontology bridge (OWL API + ELK), the integration engine (claims → posteriors), and quality metrics. `gspa-predictors` wraps external tools (DIAMOND, InterProScan, FoldSeek, eggNOG, gapseq, operon callers) behind a common `Predictor` interface. `gspa-cli` is a thin picocli driver that wires config files into one of five subcommands: `annotate`, `evaluate`, `compare`, `report`, `integrate`.

Every `java -jar gspa.jar <subcommand>` goes through `GspaMain`. It's a stock picocli root class that registers its five subcommands declaratively. The top of the file shows the wiring:

```bash
sed -n '26,50p' gspa-cli/src/main/groovy/gspa/cli/GspaMain.groovy
```

```output
@Command(
    name = 'gspa',
    description = 'Genome-Scale Protein Annotation pipeline',
    version = 'gspa 0.1.0-SNAPSHOT',
    mixinStandardHelpOptions = true,
    subcommands = [
        AnnotateCommand,
        EvaluateCommand,
        CompareCommand,
        ReportCommand,
        IntegrateCommand,
    ]
)
class GspaMain implements Runnable {

    @Override
    void run() {
        CommandLine.usage(this, System.out)
    }

    static void main(String[] args) {
        int exitCode = new CommandLine(new GspaMain()).execute(args)
        System.exit(exitCode)
    }
}
```

The integrate subcommand is defined in `IntegrateCommand.groovy`. It's the runtime target of this walkthrough: read a claims JSONL + a theta.json of hyperparameters, run the Phase 7 iterative integrator, and write one posterior per (protein, function) to TSV. Its options reveal the moving parts of the system: claims, theta, priors, dark-matter, outer loop, RLGC.

```bash
grep -n '@Option' gspa-cli/src/main/groovy/gspa/cli/IntegrateCommand.groovy | head -25
```

```output
54:    @Option(names = ['--claims'], required = true,
58:    @Option(names = ['--theta'],
62:    @Option(names = ['--out'], required = true,
66:    @Option(names = ['--provenance'],
70:    @Option(names = ['--go-owl'],
74:    @Option(names = ['--essential-functions'],
78:    @Option(names = ['--essential-profile'],
82:    @Option(names = ['--ec2go'], description = 'EC → GO mapping file.')
85:    @Option(names = ['--pathways'], description = 'Pathway definitions TSV.')
88:    @Option(names = ['--taxonomy'],
92:    @Option(names = ['--taxon-constraints'],
96:    @Option(names = ['--operons'],
100:    @Option(names = ['--gaps'],
104:    @Option(names = ['--orthogroups'],
108:    @Option(names = ['--cluster-consensus'],
112:    @Option(names = ['--enable-priors'],
116:    @Option(names = ['--lite'], description = 'Skip ELK initialization (no process coherence).')
119:    @Option(names = ['--reasoner-cache'],
123:    @Option(names = ['--dark-matter'],
127:    @Option(names = ['--suggestions-out'],
133:    @Option(names = ['--iterate-gapseq'],
137:    @Option(names = ['--max-gapseq-iter'],
141:    @Option(names = ['--gapseq-tau-cover'],
145:    @Option(names = ['--gapseq-q-base'],
149:    @Option(names = ['--gapseq-q-step'],
```

The `run()` method of IntegrateCommand is the linear narrative of the whole pipeline. It loads claims, loads theta, builds the combiner, wires reference data (GO ontology, pathway DB, operons, gaps, reaction graph) into an `IntegrationState`, builds a prior engine from the enabled list, runs an `IterativeRefiner`, optionally runs the dark-matter or RLGC suggester, and writes a TSV. Each of those steps is a class we'll open next.

A predictor produces `Annotation` objects: (type, value, score, source, evidence, goAspect). Multiple predictors can emit annotations for the same protein and same GO term. Before integration, `ClaimExtractor` lifts each Annotation into an `EvidenceClaim` — the same information plus a resolved `EvidenceType` and a calibrated probability. Claims are the currency of the integration layer.

```bash
sed -n '13,44p' gspa-core/src/main/groovy/gspa/model/Annotation.groovy
```

```output
@Canonical
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class Annotation {

    /** The annotation type */
    AnnotationType type

    /** The annotation value (e.g., GO:0006412, PF00001, EC:2.7.1.1) */
    String value

    /** Confidence score from the predictor (0.0 - 1.0) */
    double score = 0.0

    /** Which predictor produced this annotation */
    String source

    /** Evidence code (e.g., IEA, ISS, ISO for GO annotations) */
    String evidence

    /** Free-form metadata */
    Map<String, Object> metadata = [:]

    /** For GO annotations: the GO aspect (MF, BP, CC) */
    String goAspect

    /**
     * Optional evidence type classification used by the Phase 7 integration
     * layer. Predictors may set this directly; otherwise {@code ClaimExtractor}
     * falls back to a source→type lookup table.
     */
    EvidenceType evidenceType

```

```bash
sed -n '24,69p' gspa-core/src/main/groovy/gspa/integration/EvidenceClaim.groovy
```

```output
class EvidenceClaim {

    /** The protein this claim is about. */
    String proteinId

    /** Functional annotation type (GO, EC, KEGG, AMR, CAZYME, etc.). */
    AnnotationType functionType

    /** Functional identifier (e.g., GO:0006412, EC:2.7.1.1, blaCTX-M-15). */
    String functionId

    /** GO aspect (MF, BP, CC) if applicable, null otherwise. */
    String goAspect

    /** Evidence type (homology, structure, context, LM, domain-specific, ...). */
    EvidenceType evidenceType

    /** Producing predictor name, for provenance. */
    String source

    /** Raw predictor score in [0, 1] (already normalized by the predictor). */
    double rawScore

    /**
     * Calibrated probability in [0, 1]. Produced by
     * {@link CalibrationTable} from the raw score + source combination.
     * This is what enters the combiner as p = P(F | this evidence).
     */
    double calibratedProb

    /** Free-form metadata: target accession, e-value, database, etc. */
    Map<String, Object> metadata = [:]

    /** Back-reference to the original {@link Annotation} for auditability. */
    Annotation origin

    /**
     * Optional context key (e.g., operon id, pathway id) that the prior
     * engine can use to group claims.
     */
    String contextKey

    /** Key used to group claims referring to the same (protein, function). */
    String functionKey() {
        "${proteinId}|${functionType}|${functionId}".toString()
    }
```

The killer field in `EvidenceClaim` is `evidenceType`. It determines the claim's *correlation group* — which is how GSPA stops DIAMOND and eggNOG (both homology tools drawing on overlapping DBs) from being double-counted as independent evidence. Look at the enum:

```bash
sed -n '12,34p' gspa-core/src/main/groovy/gspa/integration/EvidenceType.groovy
```

```output
enum EvidenceType {

    SEQUENCE_SIMILARITY,        // DIAMOND, MMseqs2
    SEQUENCE_DOMAIN,            // InterProScan, Pfam/HMMER, TIGRFAM
    SEQUENCE_MOTIF,             // PROSITE, ELM
    SEQUENCE_DEEPLEARNING,      // DeepGO, ESM-based (Phase 9)
    STRUCTURE_SIMILARITY,       // FoldSeek
    STRUCTURE_DEEPLEARNING,     // DeepFRI, GraphGOSeq (Phase 9)
    PROTEIN_LM_EMBEDDING,       // SaProt, GOPredSim (Phase 9)
    ORTHOLOGY,                  // eggNOG-mapper, OMA
    GENOMIC_CONTEXT,            // operon co-occurrence
    METABOLIC_CONTEXT,          // gapseq pathway / gap-fill
    GENOMIC_LANGUAGE_MODEL,     // nucleotide LM over operon / regulon (Phase 9)
    LOCALIZATION,               // SignalP, DeepTMHMM
    DOMAIN_SPECIFIC_AMR,        // AMRFinder
    DOMAIN_SPECIFIC_CAZY,       // dbCAN
    DOMAIN_SPECIFIC_BGC,        // antiSMASH
    DOMAIN_SPECIFIC_VF,         // VFDB
    DARK_MATTER,                // claims promoted by DarkMatterSuggester (Phase 10)
    REACTION_LOCAL_CONTEXT,     // claims from ReactionLocalContextSuggester (Phase 12)
    CROSS_GENOME_TRANSFER,      // conditional-LR-based cross-genome transfer (Phase 12)
    ML_RANKER                   // learned ranker output (Phase 12 M3+)

```

```bash
sed -n '40,66p' gspa-core/src/main/groovy/gspa/integration/EvidenceType.groovy
```

```output
    String correlationGroup() {
        switch (this) {
            case SEQUENCE_SIMILARITY:
            case SEQUENCE_DOMAIN:
            case ORTHOLOGY:
                return 'homology'
            case SEQUENCE_DEEPLEARNING:
            case PROTEIN_LM_EMBEDDING:
                return 'ml_protein_seq'
            case STRUCTURE_SIMILARITY:
            case STRUCTURE_DEEPLEARNING:
                return 'structure'
            case GENOMIC_CONTEXT:
            case METABOLIC_CONTEXT:
                return 'context'
            case GENOMIC_LANGUAGE_MODEL:
                return 'ml_genomic'
            case DARK_MATTER:
            case REACTION_LOCAL_CONTEXT:
            case CROSS_GENOME_TRANSFER:
            case ML_RANKER:
                // Isolated: context-inferred claims share this group so
                // Phase 10 DM and Phase 12 RLGC/cross-genome/ML ranker
                // alternatives collapse correctly (we never emit more
                // than one class of context inference for the same
                // gap in production).
                return 'inferred_context'
```

So DIAMOND, HMMER, and eggNOG all live in group `homology`. FoldSeek + DeepFRI share `structure`. gapseq + operon evidence share `context`. When the combiner sees multiple claims in one group it keeps only the strongest — that's the correlation collapse.

Predictors implement one tiny interface. Everything else — tool installation checks, command building, output parsing — is in `AbstractToolPredictor`. The contract:

```bash
sed -n '12,35p' gspa-predictors/src/main/groovy/gspa/predictor/Predictor.groovy
```

```output
interface Predictor {

    /** Unique predictor name (e.g., "diamond", "interproscan", "foldseek") */
    String getName()

    /** Version string of the wrapped tool */
    String getVersion()

    /** What annotation types this predictor produces */
    Set<AnnotationType> getOutputTypes()

    /** Check if the external tool is available (installed and runnable) */
    boolean isAvailable()

    /** Predict annotations for a single protein */
    List<Annotation> predict(Protein protein)

    /** Predict annotations for a batch of proteins (default: iterate) */
    default Map<String, List<Annotation>> predictBatch(List<Protein> proteins) {
        proteins.collectEntries { p ->
            [(p.id): predict(p)]
        }
    }
}
```

Concrete wrappers set a database path and a few thresholds, then delegate command-building to `AbstractToolPredictor`. The DIAMOND wrapper is typical: a handful of fields (database, evalue, maxTargetSeqs, minIdentity), then getName()/getOutputTypes() and a buildCommand() that lays out the CLI invocation.

```bash
sed -n '12,35p' gspa-predictors/src/main/groovy/gspa/predictor/similarity/DiamondPredictor.groovy
```

```output
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
```

`ClaimExtractor` is the bridge from the predictor layer to the integration layer. It (a) resolves the evidence type using a static `source → EvidenceType` lookup table (unless the predictor already set the field), and (b) calibrates the raw score through a `CalibrationTable` (learned by the Phase 7.4 benchmark from held-out genomes). The lookup table makes it clear which predictor is mapped to which correlation group:

```bash
sed -n '22,53p' gspa-core/src/main/groovy/gspa/integration/ClaimExtractor.groovy
```

```output
    /** Source name → evidence type lookup. */
    static final Map<String, EvidenceType> SOURCE_TO_TYPE = [
        'diamond'        : EvidenceType.SEQUENCE_SIMILARITY,
        'mmseqs2'        : EvidenceType.SEQUENCE_SIMILARITY,
        'hmmer'          : EvidenceType.SEQUENCE_DOMAIN,
        'pfam'           : EvidenceType.SEQUENCE_DOMAIN,
        'interproscan'   : EvidenceType.SEQUENCE_DOMAIN,
        'foldseek'       : EvidenceType.STRUCTURE_SIMILARITY,
        'eggnog-mapper'  : EvidenceType.ORTHOLOGY,
        'eggnog'         : EvidenceType.ORTHOLOGY,
        'operon'         : EvidenceType.GENOMIC_CONTEXT,
        'gapseq'         : EvidenceType.METABOLIC_CONTEXT,
        'crossfeeding'   : EvidenceType.METABOLIC_CONTEXT,
        'signalp'        : EvidenceType.LOCALIZATION,
        'deeptmhmm'      : EvidenceType.LOCALIZATION,
        'amrfinder'      : EvidenceType.DOMAIN_SPECIFIC_AMR,
        'dbcan'          : EvidenceType.DOMAIN_SPECIFIC_CAZY,
        'antismash'      : EvidenceType.DOMAIN_SPECIFIC_BGC,
        'vfdb'           : EvidenceType.DOMAIN_SPECIFIC_VF,
        // Phase 9 predictors — reserved
        'deepgo-plus'    : EvidenceType.SEQUENCE_DEEPLEARNING,
        'deepgo'         : EvidenceType.SEQUENCE_DEEPLEARNING,
        'esm2-go'        : EvidenceType.SEQUENCE_DEEPLEARNING,
        'saprot'         : EvidenceType.PROTEIN_LM_EMBEDDING,
        'proteinbert'    : EvidenceType.PROTEIN_LM_EMBEDDING,
        'deepfri'        : EvidenceType.STRUCTURE_DEEPLEARNING,
        'graphgoseq'     : EvidenceType.STRUCTURE_DEEPLEARNING,
        'nucleotide-transformer' : EvidenceType.GENOMIC_LANGUAGE_MODEL,
        'dnabert2'       : EvidenceType.GENOMIC_LANGUAGE_MODEL,
        'evo'            : EvidenceType.GENOMIC_LANGUAGE_MODEL,
        'genslm'         : EvidenceType.GENOMIC_LANGUAGE_MODEL,
    ] as Map<String, EvidenceType>
```

```bash
sed -n '141,177p' gspa-core/src/main/groovy/gspa/integration/ClaimExtractor.groovy
```

```output
    List<EvidenceClaim> readClaimsJsonl(File file) {
        ObjectMapper mapper = new ObjectMapper()
        List<EvidenceClaim> out = []
        file.withReader { reader ->
            reader.eachLine { line ->
                line = line.trim()
                if (line.isEmpty() || line.startsWith('#')) return
                Map rec = mapper.readValue(line, Map)
                String source = rec.source as String
                double raw = (rec.raw_score as Number)?.doubleValue() ?: 0.0d
                EvidenceType type = null
                if (rec.evidence_type) {
                    try { type = EvidenceType.valueOf(rec.evidence_type as String) } catch (ignored) {}
                }
                if (type == null && source) {
                    type = SOURCE_TO_TYPE[(source as String).toLowerCase(Locale.ROOT)]
                }
                if (type == null) return    // unresolved claim; skip
                double calibrated = rec.calibrated_prob != null
                    ? ((Number) rec.calibrated_prob).doubleValue()
                    : calibration.calibrate(source ?: '', raw)
                out << new EvidenceClaim(
                    proteinId: rec.protein_id as String,
                    functionType: AnnotationType.valueOf((rec.function_type as String).toUpperCase(Locale.ROOT)),
                    functionId: rec.function_id as String,
                    goAspect: rec.go_aspect as String,
                    evidenceType: type,
                    source: source,
                    rawScore: raw,
                    calibratedProb: calibrated,
                    metadata: (rec.metadata ?: [:]) as Map,
                )
            }
        }
        out
    }
}
```

`readClaimsJsonl` is what the `integrate` subcommand calls with `--claims`. Each JSON line becomes an `EvidenceClaim` with the evidence type resolved, raw score preserved, and calibrated probability filled in. Claims without a resolvable type are silently skipped — the integrator only sees typed claims.

Now the mathematical core. Given a bundle of claims for the same (protein, function), how do we get one posterior probability? GSPA uses a Noisy-OR log-odds combination with per-type reliability weights. The reliability table is learned, but defaults are hand-picked:

```bash
sed -n '51,71p' gspa-core/src/main/groovy/gspa/integration/EvidenceCombiner.groovy
```

```output
    static Map<EvidenceType, Double> defaultReliability() {
        Map<EvidenceType, Double> r = new EnumMap<>(EvidenceType)
        EvidenceType.values().each { r[it] = 0.6d }
        r[EvidenceType.SEQUENCE_SIMILARITY]      = 0.70d
        r[EvidenceType.SEQUENCE_DOMAIN]          = 0.75d
        r[EvidenceType.STRUCTURE_SIMILARITY]     = 0.80d
        r[EvidenceType.STRUCTURE_DEEPLEARNING]   = 0.75d
        r[EvidenceType.ORTHOLOGY]                = 0.70d
        r[EvidenceType.SEQUENCE_DEEPLEARNING]    = 0.65d
        r[EvidenceType.PROTEIN_LM_EMBEDDING]     = 0.65d
        r[EvidenceType.GENOMIC_CONTEXT]          = 0.45d
        r[EvidenceType.METABOLIC_CONTEXT]        = 0.55d
        r[EvidenceType.GENOMIC_LANGUAGE_MODEL]   = 0.55d
        r[EvidenceType.LOCALIZATION]             = 0.50d
        r[EvidenceType.DOMAIN_SPECIFIC_AMR]      = 0.85d
        r[EvidenceType.DOMAIN_SPECIFIC_CAZY]     = 0.80d
        r[EvidenceType.DOMAIN_SPECIFIC_BGC]      = 0.75d
        r[EvidenceType.DOMAIN_SPECIFIC_VF]       = 0.75d
        r[EvidenceType.SEQUENCE_MOTIF]           = 0.50d
        r
    }
```

```bash
sed -n '81,130p' gspa-core/src/main/groovy/gspa/integration/EvidenceCombiner.groovy
```

```output
    double combineLikelihood(List<EvidenceClaim> claimsForSameFunction) {
        if (claimsForSameFunction == null || claimsForSameFunction.isEmpty()) {
            return lMin
        }

        // Step 1 + 2: group by correlation group, pick the strongest in each group.
        Map<String, List<EvidenceClaim>> byGroup = [:]
        for (EvidenceClaim c : claimsForSameFunction) {
            String g = c.evidenceType.correlationGroup()
            byGroup.computeIfAbsent(g, { k -> new ArrayList<EvidenceClaim>() }) << c
        }

        double oneMinusProd = 1.0d
        for (Map.Entry<String, List<EvidenceClaim>> entry : byGroup.entrySet()) {
            EvidenceClaim best = null
            double bestScore = -1.0d
            EvidenceClaim secondBest = null
            double secondBestScore = -1.0d
            for (EvidenceClaim c : entry.value) {
                double w = reliability.getOrDefault(c.evidenceType, 0.5d)
                double s = w * c.calibratedProb
                if (s > bestScore) {
                    secondBest = best
                    secondBestScore = bestScore
                    best = c
                    bestScore = s
                } else if (s > secondBestScore) {
                    secondBest = c
                    secondBestScore = s
                }
            }

            double p = bestScore
            // Second-opinion bonus: only if the two claims are from distinct databases.
            if (secondOpinionBonus > 0.0 && secondBest != null
                    && distinctDatabases(best, secondBest)) {
                p += secondOpinionBonus * secondBestScore
            }
            p = Math.min(0.999d, Math.max(0.001d, p))
            oneMinusProd *= (1.0d - p)
        }

        // Step 3: Noisy-OR combination across groups.
        double pPost = 1.0d - oneMinusProd
        pPost = Math.min(0.999999d, Math.max(1.0e-6d, pPost))

        // Step 4: convert to log-odds + clip.
        double logOdds = Math.log(pPost / (1.0d - pPost))
        Math.min(lMax, Math.max(lMin, logOdds))
    }
```

Read the code top-down: (1) group claims by correlation group; (2) within each group keep the strongest weighted-by-reliability claim, optionally adding a small 'second-opinion bonus' when a distinct DB corroborates; (3) apply Noisy-OR across groups — 1 − ∏(1 − p_g); (4) convert to log-odds and clip to ±12. The result is a single scalar 'likelihood log-odds' per (protein, function) — the starting point for iterative refinement.

Likelihood alone isn't enough — GSPA also injects contextual priors (does this annotation complete a pathway? is it taxonomically valid? is the protein in a triggered operon?). `IterativeRefiner` runs a Jacobi-style fixed-point loop: posterior = likelihood + Σ λ_k · prior_k(boost), damping each step so non-monotone priors don't oscillate.

```bash
sed -n '29,45p' gspa-core/src/main/groovy/gspa/integration/IterativeRefiner.groovy
```

```output
class IterativeRefiner {

    private static final Logger log = LoggerFactory.getLogger(IterativeRefiner)

    EvidenceCombiner combiner
    PriorEngine priorEngine = new PriorEngine()

    int maxIter = 6
    double epsilon = 0.005
    double damping = 0.5

    IterativeRefiner(EvidenceCombiner combiner) {
        this.combiner = combiner
    }

    /**
     * Refine a set of claims into posterior annotations.
```

```bash
sed -n '71,100p' gspa-core/src/main/groovy/gspa/integration/IterativeRefiner.groovy
```

```output
        for (int iter = 0; iter < maxIterations; iter++) {
            priorEngine.beginIteration(state)

            Map<String, Double> newLogOdds = new LinkedHashMap<>()
            for (Map.Entry<String, List<EvidenceClaim>> entry : byKey.entrySet()) {
                String key = entry.key
                double lLik = likelihood[key]
                double lPri = priorEngine.totalBoost(entry.value.first().proteinId, key, state)
                double lNew = clip(lLik + lPri, key, state)

                double lOld = posteriorLogOdds.getOrDefault(key, lLik)
                // Jacobi-style under-relaxation: new = (1-d) * old + d * computed.
                double lDamped = (1.0d - damping) * lOld + damping * lNew
                // Re-apply the pin floor after damping. Damping can drag a
                // floored lNew back toward an unfloored lOld (seed = raw
                // likelihood on iteration 1), which would silently undo the
                // pin. Clamp to [floor, lMax] if a floor exists.
                Double floor = state.pinnedFloors != null ? state.pinnedFloors[key] : null
                if (floor != null && lDamped < floor) {
                    lDamped = Math.min(combiner.lMax, (double) floor)
                }
                newLogOdds[key] = lDamped
            }

            double delta = meanAbsDelta(newLogOdds, posteriorLogOdds)
            log.debug("iter ${iter}: mean |Δp| = ${String.format(Locale.ROOT, '%.5f', delta)}")

            posteriorLogOdds = newLogOdds
            state.updatePosteriors(posteriorLogOdds)
            iterationsRun = iter + 1
```

Default: up to 6 iterations, damping 0.5, convergence ε = 0.005 on mean |Δp|. At each iteration the engine calls `beginIteration` on every prior (so expensive priors like Consistency and Coherence cache their state), then recomputes `l_new = l_likelihood + Σ λ_k · prior_k.boost`, damps it, and reapplies any 'pinned floor' from the Phase 10 outer loop. Divergence detection rolls back to the last stable state.

Each prior is a tiny plugin conforming to this interface:

```bash
sed -n '10,32p' gspa-core/src/main/groovy/gspa/integration/Prior.groovy
```

```output
interface Prior {
    String name()

    /**
     * Called once at the start of each refinement iteration.
     *
     * Expensive priors (Consistency, Coherence) should recompute their
     * per-iteration state here (e.g. run the SAT solver, compute the
     * currently-annotated set) so that {@link #logOddsBoost} becomes a
     * cheap lookup. Default: no-op.
     */
    default void beginIteration(IntegrationState state) {}

    /**
     * Log-odds contribution for (proteinId, functionKey). Zero means no
     * effect; negative values downweight the claim; positive values boost.
     */
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state)

    /** Declares which parts of state this prior reads, for re-run triggers. */
    Set<String> inputs()
}
```

GSPA ships six concrete priors: `Essentiality`, `Coherence`, `Consistency`, `GapFilling`, `GenomicContext`, and `HomologyTransfer` (for cross-genome). Each contributes a signed log-odds boost per (protein, function). The `PriorEngine` just sums them with per-prior λ weights from theta.json.

The `CoherencePrior` is ontology-driven. In each iteration it asks the ELK reasoner: which (whole, part) pairs in GO's has_part hierarchy have a 'whole' annotated but the required 'part' missing? Each such unsatisfied F-term earns a boost on any candidate claim (protein, F) — scaled by how close the pathway already is to being complete. 'Closed stays closed', but 'nearly closed' gets promoted.

```bash
sed -n '54,82p' gspa-core/src/main/groovy/gspa/integration/prior/CoherencePrior.groovy
```

```output
    void beginIteration(IntegrationState state) {
        processMissingTerms = new LinkedHashMap<>()
        pathwayMissingTerms = new LinkedHashMap<>()

        Set<String> annotated = state.goReasoner != null
            ? state.currentlyAnnotatedGoTermsPropagated()
            : state.currentlyAnnotatedGoTerms()

        // --- Process coherence: needs GoReasoner for has_part pairs ---
        if (state.goReasoner != null) try {
            def coherence = new Coherence(state.goOntology, state.goReasoner)
            ProcessCoherenceResult result = coherence.evaluateProcessCoherence(annotated)
            int triggered = result.triggered
            int satisfied = result.satisfied
            double fracAnnotated = triggered == 0 ? 1.0 : (satisfied / (double) triggered)
            double weight = Math.max(0.0d, 1.0d - fracAnnotated)

            for (Map.Entry<String, String> unsat : result.unsatisfied) {
                // Each entry: C present but F missing. Boost F.
                String fTerm = unsat.value
                if (fTerm != null) {
                    processMissingTerms.merge(fTerm, weight, { a, b -> Math.max(a as double, b as double) })
                }
            }
        } catch (Exception ex) {
            log.warn("CoherencePrior: process-coherence check failed: ${ex.message}")
        }

        // --- Pathway coherence: for each pathway, find missing required terms ---
```

`ConsistencyPrior` is where the codebase's signature design decision shows up. Taxon constraints ('never_in_taxon', 'only_in_taxon') are negative statements — 'no photosynthesis in Archaea' — which ELK cannot handle because it's a tractable OWL 2 EL reasoner with no disjointness + negation. So GSPA encodes each iteration's currently-annotated set as a SAT instance over propositional taxon variables, runs SAT4J, and uses the UNSAT core as the set of 'conflicting GO terms'. Any claim about a conflicting term gets a negative log-odds penalty.

```bash
sed -n '48,82p' gspa-core/src/main/groovy/gspa/integration/prior/ConsistencyPrior.groovy
```

```output
    void beginIteration(IntegrationState state) {
        conflictingTerms = Collections.emptySet()
        SatConsistencyChecker checker = state.satConsistencyChecker
        if (checker == null) {
            log.debug('ConsistencyPrior disabled: no SatConsistencyChecker in state')
            return
        }

        Set<String> annotated = state.currentlyAnnotatedGoTerms()
        if (annotated.isEmpty()) return

        try {
            ConsistencyResult result = checker.check(annotated)
            if (result.consistent) return
            Set<String> terms = new LinkedHashSet<>()
            for (ConsistencyViolation v : result.violations ?: []) {
                if (v.involvedGoTerms != null) terms.addAll(v.involvedGoTerms)
            }
            conflictingTerms = terms
            log.debug("ConsistencyPrior: ${terms.size()} conflicting GO terms from ${result.violations?.size() ?: 0} violations")
        } catch (Exception ex) {
            log.warn("ConsistencyPrior: SAT check failed: ${ex.message}")
        }
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (conflictingTerms.isEmpty()) return 0.0d
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null || parts[1] != 'GO') return 0.0d
        if (!conflictingTerms.contains(parts[2])) return 0.0d
        double penalty = hardFilter ? -1000.0d : -alphaCons
        return penalty
    }
}
```

Default penalty is soft: `-alphaCons` log-odds (~3.0). Strong likelihoods can survive one soft penalty — the right call for HGT and contamination cases. `hardFilter = true` swaps in -1000 and effectively removes the claim.

Priors can reweight existing claims — but what about proteins with no annotation at all? Phase 8's `DarkMatterSuggester` handles that. Given a metabolic gap (Pathway P, Reaction R missing the gene for function f_R) and the genome's operon structure, it assigns the missing function to the most likely protein in the most likely operon using a four-layer Bayesian score:

```bash
sed -n '12,40p' gspa-core/src/main/groovy/gspa/integration/suggester/DarkMatterSuggester.groovy
```

```output
/**
 * Phase 8 "Dark Matter / Contextual Gap" Suggester.
 *
 * <p>Takes a metabolic gap {@code (P, R)} and the genomic context
 * (operons + integrated posteriors) and assigns the missing function to
 * the most likely protein in the most likely operon. Implements the
 * four-layer Bayesian algorithm from plan §A.1:</p>
 *
 * <ol>
 *   <li><b>Layer 1</b> — Bayes factor {@code BF(O, P)} that operon O
 *       participates in pathway P, computed as a soft-weighted product
 *       over operon members' annotations (no "present / absent"
 *       threshold).</li>
 *   <li><b>Layer 2</b> — per protein
 *       {@code L_R(p) = L_lik + L_op + L_lm}, where {@code L_op =
 *       log(BF/(1+BF))} is shared across the operon and {@code L_lik}
 *       is the integrator's current posterior for {@code (p, f_R)} or a
 *       low "absent" default.</li>
 *   <li><b>Layer 3</b> — softmax over operon members:
 *       {@code q(p) = π_R(p) / Σ π_R(p')}.</li>
 *   <li><b>Layer 4</b> — singleton if the top q &gt; 0.5, otherwise a
 *       disjunctive suggestion over the smallest top-k whose cumulative
 *       q exceeds {@code coverageThreshold}.</li>
 * </ol>
 *
 * <p>Phase 8 is boost/suggest-only — suggestions never re-enter the
 * refinement loop. They are written to
 * {@link IntegratedAnnotationSet#suggestions} as a separate channel.</p>
 */
```

Suggestions are emitted as either a `SingletonSuggestion` (top q > 0.5, one confident protein) or a `DisjunctiveSuggestion` (smallest top-k whose cumulative q reaches the coverage threshold — 'one of these three genes probably does it'). Suggestions are a separate channel from posteriors; they don't feed back into the refiner in Phase 8. Phase 10's `OuterIterativeRefiner` promotes high-q singleton suggestions into new DARK_MATTER claims and re-runs the refiner with them pinned.

Phase 12's `ReactionLocalContextSuggester` (RLGC) is the newer, operon-free alternative. Instead of requiring a pre-called operon structure, it uses a continuous Gaussian kernel density over genome coordinates plus BFS on the gapsmith reaction graph. This lets it work on genomes with poorly-predicted operons and partial synteny.

```bash
sed -n '14,42p' gspa-core/src/main/groovy/gspa/integration/suggester/ReactionLocalContextSuggester.groovy
```

```output
/**
 * Phase 12 Reaction-Local Genomic Context (RLGC) Suggester.
 *
 * <p>Replaces {@link DarkMatterSuggester}'s pathway × operon × BF-gate
 * machinery with a data-driven triple:</p>
 * <ol>
 *   <li><b>Reaction-graph locality</b> — score by evidence for R's
 *       reaction-graph neighbours (distance-weighted via {@code alpha^d})
 *       rather than pathway membership.</li>
 *   <li><b>Continuous genomic density</b> — kernel-smoothed field over
 *       anchor (neighbour-catalysing) proteins; no operon call required.</li>
 *   <li><b>Diversity / commitment / direction priors</b> — Noisy-OR over
 *       distinct neighbour reactions, commitment penalty from current
 *       posteriors, strand-consistency, and intergenic-gap penalty.</li>
 * </ol>
 *
 * <p>For each gap {@code (P, R, f_R)}, the algorithm:</p>
 * <ol>
 *   <li>Builds {@code N_k(R)} via {@link ReactionGraph#bfs}.</li>
 *   <li>Collects anchor proteins: those with posterior &gt; τ for the
 *       GO term of any neighbour r' ∈ N_k(R).</li>
 *   <li>Constructs a {@link GenomicDensityField} of per-anchor weights.</li>
 *   <li>Ranks candidate genes on {@code log D + β·Div − γ·commitment −
 *       δ·self + ε·strand − ζ·log(intergenic_gap)}.</li>
 *   <li>Softmaxes the top-Q candidates; emits singleton if p > 0.5,
 *       otherwise disjunctive over top-k covering {@code coverageThreshold}.</li>
 * </ol>
 */
class ReactionLocalContextSuggester {
```

In the CLI this is gated by `--rlc-suggester` and needs two extra inputs: `--reaction-graph` (gapsmith seed_reactions TSV) and `--genome-layout` (a TSV of protein_id / contig / start / end / strand). The output format is the same as DarkMatter — Singleton or Disjunctive — so both can flow into the same downstream pipeline, and the M2 cross-genome re-scorer can layer on top of either.

```bash
ls gspa-core/src/main/groovy/gspa/metrics/
```

```output
Coherence.groovy
Completeness.groovy
Consistency.groovy
HtmlReportWriter.groovy
InformationContent.groovy
MagAdjuster.groovy
QualityPipeline.groovy
QualityReportWriter.groovy
QualityScorer.groovy
```

Separate from integration is GAEF — the Genome Annotation Evaluation Framework — in `gspa/metrics/`. Three axes: Completeness (essential-function coverage), Coherence (process / pathway / complex closure via has_part), and Consistency (taxon constraints via SAT4J, the same technique as the ConsistencyPrior). `QualityScorer` rolls them into a composite score; `QualityPipeline` orchestrates; `HtmlReportWriter` renders the dashboards. The `gspa evaluate` subcommand is the user-facing driver:

```bash
sed -n '155,165p' README.md
```

````output
# Convert integrated TSV to synthetic GFF + GAF for gspa evaluate
python3 benchmark/make_synth_gff_gaf.py \
  --fasta bsubtilis.faa --tsv integrated.tsv \
  --gff-out synth.gff --gaf-out synth.gaf

java -jar gspa.jar evaluate \
  -i synth.gff -a synth.gaf --go-owl go.owl \
  --ec2go ec2go --pathways kegg_pathways.tsv \
  --reasoner-cache ./reasoner-cache \
  -k bacteria -o quality.json
```
````

Three pillars, in short. (1) **Evidence integration** collapses correlated predictors into independent groups, then combines them with a reliability-weighted Noisy-OR. (2) **Priors-as-log-odds** adds pathway-coherence boosts and taxon-consistency penalties on each iteration of a damped fixed-point loop, with ELK for subsumption and SAT4J for anything involving negation. (3) **Dark-matter suggestions** reach beyond the claimed proteins: DarkMatterSuggester for operon-aware genomes, ReactionLocalContextSuggester for operon-free reaction-graph + kernel-density inference. The outer loop can feed high-confidence singleton suggestions back as new claims, iterating until the gap set stabilizes. The product is a TSV of (protein, function) posteriors with per-prior provenance, plus a separate suggestions channel for downstream metabolic modeling.

