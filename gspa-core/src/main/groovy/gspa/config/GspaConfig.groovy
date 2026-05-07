package gspa.config

import gspa.model.OrganismDomain
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * Central configuration for the GSPA pipeline.
 * Populated from YAML config files with hierarchical merging:
 * built-in defaults -> kingdom preset -> user config -> CLI flags.
 */
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class GspaConfig {

    // --- Input ---
    InputConfig input = new InputConfig()

    // --- Predictors ---
    PredictorConfig predictors = new PredictorConfig()

    // --- Quality assessment ---
    QualityConfig quality = new QualityConfig()

    // --- Output ---
    OutputConfig output = new OutputConfig()

    // --- Evidence integration (Phase 7+) ---
    IntegrationConfig integration = new IntegrationConfig()

    static class InputConfig {
        /** Input format: fasta, genbank, gff */
        String format = 'fasta'
        /** Input type: genome, mag, metagenome */
        String type = 'genome'
        /** Organism kingdom: bacteria, archaea, eukaryote, virus, auto */
        String kingdom = 'auto'
        /** Path to protein FASTA (if pre-called) */
        String proteinFasta
    }

    static class PredictorConfig {
        SimilarityConfig similarity = new SimilarityConfig()
        StructureConfig structure = new StructureConfig()
        DomainConfig domains = new DomainConfig()
        OperonConfig operons = new OperonConfig()
        PathwayConfig pathway = new PathwayConfig()
        LocalizationConfig localization = new LocalizationConfig()
        DisorderConfig disorder = new DisorderConfig()
        NeuralConfig neural = new NeuralConfig()
        FossPredictorConfig foss = new FossPredictorConfig()
        ViralConfig viral = new ViralConfig()

        /** Additional predictors to enable (by name) */
        Set<String> enable = []
        /** Predictors to explicitly disable */
        Set<String> disable = []
    }

    static class SimilarityConfig {
        boolean enabled = true
        /** Tool: diamond, mmseqs2 */
        String tool = 'diamond'
        /** Path to reference database */
        String database
        double evalue = 1e-5
        int maxTargetSeqs = 10
    }

    static class StructureConfig {
        boolean enabled = false
        /** Tool: foldseek */
        String tool = 'foldseek'
        /** Path to structure database (e.g., AlphaFold DB) */
        String database
        /** Whether to predict structures first (via ESMFold) */
        boolean predictStructures = false
        /** Path to ProstT5 model directory for sequence→structure search */
        String prostt5Model
        /**
         * Path to a FoldSeek function-centroid database (built by
         * {@code benchmark/neural/build_foldseek_centroids.py}). When set
         * together with a {@link #centroidMode} other than {@code "none"},
         * FoldSeek runs in centroid mode against this DB instead of doing
         * homology transfer against {@link #database}.
         */
        String centroidDb
        /** Centroid mode: none | go | ec | both (default none). */
        String centroidMode = 'none'
    }

    static class DomainConfig {
        boolean enabled = true
        /** InterProScan applications to run */
        List<String> applications = ['Pfam', 'TIGRFAM', 'CDD', 'SUPERFAMILY']
        /** Use InterProScan or direct HMMER */
        String tool = 'interproscan'
    }

    static class OperonConfig {
        boolean enabled = true
        /** Maximum intergenic distance (bp) to consider genes as co-operonic */
        int maxIntergenicDistance = 300
        /** Require genes to be on the same strand */
        boolean requireSameStrand = true
    }

    static class PathwayConfig {
        boolean enabled = true
        /** Tool: gapseq */
        String tool = 'gapseq'
        /** Enable crossfeeding analysis for MAGs/metagenomes */
        boolean crossfeeding = false
        /** Export SBML models */
        boolean exportSbml = false
        /** Run FBA via COBRApy (requires Python) */
        boolean runFba = false
    }

    static class LocalizationConfig {
        boolean enabled = false
        boolean signalP = true   // licensed; FOSS replacement: deepSig
        boolean tmhmm = true     // licensed; FOSS replacement: tmbed
        boolean psort = false
        boolean deepSig = false  // FOSS Sec/Tat signal peptide
        boolean tmbed = false    // FOSS TM helices via ProtT5
        boolean tppred3 = false  // FOSS N-terminal targeting peptide
        boolean psortb = false   // FOSS bacterial localization
        String deepSigKingdom = 'gramn'   // gramp | gramn | euk
        String tppred3Kingdom = 'nonplant' // plant | nonplant
        String psortbGram = 'negative'    // positive | negative | archaea
    }

    /** Intrinsic-disorder prediction (Metapredict). */
    static class DisorderConfig {
        boolean enabled = false
        /** Minimum residue length of an IDR to emit. */
        int minRegionLen = 10
        /** Minimum mean disorder score to emit. */
        double minScore = 0.5
    }

    /**
     * Neural / deep-learning predictors running via the Python sidecar at
     * {@code benchmark/neural/run_neural_predictors.py}.
     */
    static class NeuralConfig {
        /** Absolute path to {@code run_neural_predictors.py}. Required for
         *  any of the neural predictors. */
        String sidecarScript
        /** Python executable used to invoke the sidecar. */
        String pythonExecutable = 'python3'

        Esm2DeepGoPlusConfig esm2DeepGoPlus = new Esm2DeepGoPlusConfig()
        ProteInferConfig proteinfer = new ProteInferConfig()
        CleanNeuralConfig clean = new CleanNeuralConfig()
        Esm2CentroidConfig esm2Centroid = new Esm2CentroidConfig()
    }

    static class Esm2DeepGoPlusConfig {
        boolean enabled = false
        /** Trained ESM2Head checkpoint (.pt). Required when enabled. */
        String checkpoint
        /** Vocabulary file: one GO term per line, column index = FC output idx. */
        String terms
        /** ESM2 variant matching the checkpoint. */
        String model = 'esm2_t33_650M_UR50D'
        int batchSize = 16
        double minScore = 0.1
    }

    static class ProteInferConfig {
        boolean enabled = false
        /** Directory with a ProteInfer model release. */
        String modelDir
        int batchSize = 16
        double minScore = 0.1
    }

    static class CleanNeuralConfig {
        boolean enabled = false
        /** Directory with a CLEAN checkpoint (contains CLEAN_infer_fasta.py). */
        String modelDir
        int batchSize = 16
        double minScore = 0.1
    }

    static class Esm2CentroidConfig {
        boolean enabled = false
        /** NPZ file with centroids / terms / annotation_types arrays. */
        String db
        /** ESM2 variant used to embed queries at inference time. */
        String model = 'esm2_t12_35M_UR50D'
        /** Top-k centroid neighbours to keep per protein. */
        int topK = 5
        int batchSize = 16
        double minScore = 0.2
    }

    /**
     * FOSS region/term/site predictors that delegate to the
     * {@code run_region_predictors.py}, {@code run_term_predictors.py},
     * and {@code run_site_predictors.py} sidecars under {@code benchmark/neural/}.
     */
    static class FossPredictorConfig {
        /** Absolute paths to the three sidecar scripts. */
        String regionSidecar
        String termSidecar
        String siteSidecar
        String pythonExecutable = 'python3'

        DeepFriConfig    deepFri    = new DeepFriConfig()
        DeepEcConfig     deepEc     = new DeepEcConfig()
        DeepArgConfig    deepArg    = new DeepArgConfig()
        MusiteDeepConfig musiteDeep = new MusiteDeepConfig()
        ScanNetConfig    scanNet    = new ScanNetConfig()
        MdFConfig        mdf        = new MdFConfig()
    }

    static class DeepFriConfig {
        boolean enabled = false
        /** Path to DeepFRI repo clone (predict.py + bundled weights). */
        String modelDir
        /** Mode: 'seq' (default) or 'struct'. */
        String mode = 'seq'
        double minScore = 0.5
    }

    /**
     * metagenomic-deepFRI (Bezshapkin et al. 2026, BSD-3-Clause).
     * Successor to {@link DeepFriConfig} with FoldComp-database
     * structure retrieval + ONNX inference. Wired through the term
     * sidecar {@code run_term_predictors.py --predictor mdf}.
     */
    static class MdFConfig {
        boolean enabled = false
        /** Path to mDeepFRI weights directory (from
         *  {@code mDeepFRI get-models --version 1.0 --output ...}). */
        String weightsDir
        /** Path to the {@code mDeepFRI} binary; null = first on PATH. */
        String mdfExecutable
        /** Sequence-only path (default). False to use structure path; requires {@code foldcompDb}. */
        boolean skipPdb = true
        /** FoldComp-format structure DB (AFDBv4, ESM Atlas, etc.). */
        String foldcompDb
        int threads = 4
        double minScore = 0.1
    }

    static class DeepEcConfig {
        boolean enabled = false
        /** Path to DeepEC repo clone. */
        String modelDir
        double minScore = 0.5
        /** AGPL-3.0 acknowledgement. Set true to acknowledge network-clause. */
        boolean acknowledgeAgpl = false
    }

    static class DeepArgConfig {
        boolean enabled = false
        /** Path to DeepARG database directory. */
        String modelDir
        /** prot | nucl */
        String type = 'prot'
        double minScore = 0.5
    }

    static class MusiteDeepConfig {
        boolean enabled = false
        /** Path to MusiteDeep_web repo clone. */
        String modelDir
        /** Underscore-joined PTM types matching MusiteDeep --residue-types. */
        String residueTypes = 'Phosphoserine_Phosphothreonine'
        double minScore = 0.5
    }

    static class ScanNetConfig {
        boolean enabled = false
        String modelDir
        /** Directory of {@code <tag>/*.pdb} structure files. */
        String structureDir
        double minScore = 0.5
    }

    /**
     * v1.3 phage / prophage / viral predictors. All FOSS; consume
     * {@code genome_fasta} (not protein FASTA).
     */
    static class ViralConfig {
        /** Absolute path to {@code run_genomic_predictors.py}. */
        String genomicSidecar
        String pythonExecutable = 'python3'

        GenomadConfig genomad = new GenomadConfig()
        CheckVConfig  checkv  = new CheckVConfig()
        PhiSpyConfig  phispy  = new PhiSpyConfig()
    }

    static class GenomadConfig {
        boolean enabled = false
        String dbPath          // path to geNomad model DB
        String sif             // optional Singularity image
        double minScore = 0.7
    }

    static class CheckVConfig {
        boolean enabled = false
        String dbPath
        String sif
        int threads = 4
        double minScore = 0.5
    }

    static class PhiSpyConfig {
        boolean enabled = false
        String sif
        String trainset        // optional path to a custom trainset
        double minScore = 0.5
    }

    static class QualityConfig {
        CompletenessConfig completeness = new CompletenessConfig()
        CoherenceConfig coherence = new CoherenceConfig()
        ConsistencyConfig consistency = new ConsistencyConfig()
        /** Path to GO OWL file */
        String goOwlFile
    }

    static class CompletenessConfig {
        /** Essential function profile name or path to custom TSV */
        String profile = 'auto'
        /** Additional GO terms to add to the profile for this run */
        List<String> addTerms = []
        /** GO terms to remove from the profile for this run */
        List<String> removeTerms = []
        /** Path to custom essential functions file (replaces profile) */
        String customFile
    }

    static class CoherenceConfig {
        boolean process = true
        boolean pathway = true
        boolean complex = true
    }

    static class ConsistencyConfig {
        boolean taxonConstraints = true
        /** Path to taxonomy hierarchy file */
        String taxonomyFile
    }

    static class OutputConfig {
        /** Output formats to generate */
        List<String> formats = ['gff3', 'tsv']
        /** Generate HTML quality report */
        boolean htmlReport = true
        /** Output directory */
        String outputDir = 'gspa_output'
    }

    /**
     * Phase 7 evidence-integration configuration. Disabled by default so
     * existing pipelines are unaffected.
     */
    static class IntegrationConfig {
        /** Master switch for the integration layer. */
        boolean enabled = false

        /** Max refinement iterations (Phase 7.3+). */
        int maxIter = 6

        /** Convergence threshold on mean |Δp|. */
        double epsilon = 0.005

        /** Damping / under-relaxation factor. */
        double damping = 0.5

        /** Path to a learned theta.json with reliability + prior weights. */
        String thetaFile = null

        /** Path to a YAML override file with integration defaults. */
        String defaultsFile = null

        /** Dark-matter suggester (Phase 8). */
        DarkMatterConfig darkMatter = new DarkMatterConfig()

        /** Phase 10 intra-genome clustering (cluster proteome, predict reps only, propagate). */
        IntragenomeClusterConfig intragenomeCluster = new IntragenomeClusterConfig()

        /** Phase 10 outer iterative gapseq+DarkMatter fixed-point loop. */
        OuterLoopConfig outerLoop = new OuterLoopConfig()
    }

    /** Phase 10 intra-genome clustering config. Disabled by default. */
    static class IntragenomeClusterConfig {
        boolean enabled = false
        /** Minimum sequence identity (0.9 = 90%). */
        double identity = 0.9d
        /** Minimum alignment coverage (0.8 = 80%). */
        double coverage = 0.8d
    }

    /** Phase 10 outer-loop config (iterative gapseq + DarkMatter). */
    static class OuterLoopConfig {
        boolean enabled = false
        int maxIter = 5
        double qBase = 0.5d
        double qStep = 0.05d
        double qCap = 0.75d
        double tauCover = 0.5d
        boolean pinPromotions = true
    }

    /** Phase 8 dark-matter suggester config. Disabled by default. */
    static class DarkMatterConfig {
        boolean enabled = false

        /** Output mode: separate | top-q | distributed. */
        String mode = 'separate'

        /** Bayes factor threshold for candidate operons. */
        double bfMin = 10.0

        /** Likelihood elevation factor for "in pathway" functions. */
        double gammaInP = 50.0

        /** Credible-set coverage for disjunctive output. */
        double coverageThreshold = 0.9

        /** Enable SAT-based disambiguation (Phase 8.3). */
        boolean satDisambiguation = false

        /** Use Phase 9 genomic-LM context in suggester. */
        boolean useLmContext = false
    }

    /**
     * Resolve the organism domain from config or auto-detect.
     */
    OrganismDomain resolveOrganismDomain() {
        if (input.kingdom == 'auto') return OrganismDomain.UNKNOWN
        return OrganismDomain.fromName(input.kingdom)
    }

    /**
     * Resolve the essential function profile name based on kingdom.
     */
    String resolveCompletenessProfile() {
        if (quality.completeness.customFile) return 'custom'
        if (quality.completeness.profile != 'auto') return quality.completeness.profile
        switch (resolveOrganismDomain()) {
            case OrganismDomain.BACTERIA: return 'bacteria'
            case OrganismDomain.ARCHAEA: return 'archaea'
            case OrganismDomain.EUKARYA: return 'eukaryote'
            default: return 'bacteria'
        }
    }
}
