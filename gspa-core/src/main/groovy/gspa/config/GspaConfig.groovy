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
        boolean signalP = true
        boolean tmhmm = true
        boolean psort = false
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
