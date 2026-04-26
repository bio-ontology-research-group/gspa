package gspa.integration

/**
 * The type of evidence that supports an {@link EvidenceClaim}.
 *
 * Types are grouped by correlation: claims whose evidence types share a
 * {@link #correlationGroup} are not combined as independent by
 * {@code EvidenceCombiner} — only the strongest claim per group enters the
 * Noisy-OR combination. This addresses the non-independence of tools like
 * DIAMOND and eggNOG that draw from overlapping homology databases.
 */
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
    ML_RANKER,                  // learned ranker output (Phase 12 M3+)
    SEQUENCE_REGION_ML          // region-level ML predictors (Metapredict, SignalP region, TMHMM helix)

    /**
     * Correlation group. Claims in the same group are collapsed to the single
     * strongest claim before Noisy-OR combination, since their evidence is
     * drawn from overlapping sources and is not truly independent.
     */
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
            case SEQUENCE_MOTIF:
                return 'motif'
            case LOCALIZATION:
                return 'localization'
            case SEQUENCE_REGION_ML:
                // Region-level ML predictors (disorder, signal peptide
                // cleavage site, TM helix) predict positional features
                // that do not compete with whole-protein homology or
                // structure evidence, so they get their own group.
                return 'region_features'
            case DOMAIN_SPECIFIC_AMR:
            case DOMAIN_SPECIFIC_CAZY:
            case DOMAIN_SPECIFIC_BGC:
            case DOMAIN_SPECIFIC_VF:
                // Each domain-specific type is its own group; they speak
                // about different function namespaces (AMR gene IDs vs
                // CAZyme families etc.) so correlation isn't an issue.
                return "domain_${name().toLowerCase(Locale.ROOT)}".toString()
            default:
                return name().toLowerCase(Locale.ROOT)
        }
    }
}
