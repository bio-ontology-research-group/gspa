package gspa.integration

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.model.Protein

/**
 * Converts {@link Annotation} instances attached to a {@link Genome}'s
 * proteins into {@link EvidenceClaim}s for the Phase 7 integration layer.
 *
 * Each predictor's {@code source} string is mapped to an {@link EvidenceType}
 * via a static lookup table. Predictors that wish to override the default
 * can populate {@link Annotation#evidenceType} directly — the extractor
 * prefers the annotation's own value when set.
 *
 * Raw scores are calibrated through {@link CalibrationTable}.
 */
class ClaimExtractor {

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

    final CalibrationTable calibration

    ClaimExtractor(CalibrationTable calibration = new CalibrationTable()) {
        this.calibration = calibration
    }

    /**
     * Resolve the evidence type for an annotation. Prefers the annotation's
     * explicit {@code evidenceType}, then falls back to the source lookup
     * table. Returns {@code null} if neither is set; caller must decide
     * whether to skip or default.
     */
    static EvidenceType resolveType(Annotation a) {
        if (a.evidenceType != null) return a.evidenceType
        if (a.source != null) return SOURCE_TO_TYPE[a.source.toLowerCase(Locale.ROOT)]
        return null
    }

    /**
     * Convert a single annotation on a protein into a claim.
     * Returns {@code null} if the evidence type cannot be resolved — the
     * caller should log a warning and drop the claim.
     */
    EvidenceClaim fromAnnotation(String proteinId, Annotation ann) {
        EvidenceType type = resolveType(ann)
        if (type == null) return null
        new EvidenceClaim(
            proteinId: proteinId,
            functionType: ann.type,
            functionId: ann.value,
            goAspect: ann.goAspect,
            evidenceType: type,
            source: ann.source,
            rawScore: ann.score,
            calibratedProb: calibration.calibrate(ann.source ?: '', ann.score),
            metadata: ann.metadata != null ? new LinkedHashMap<>(ann.metadata) : [:],
            origin: ann,
        )
    }

    /**
     * Extract all claims from a protein's annotations. Annotations whose
     * evidence type cannot be resolved are skipped.
     */
    List<EvidenceClaim> fromProtein(Protein protein) {
        if (protein.annotations == null) return []
        List<EvidenceClaim> out = []
        for (Annotation ann : protein.annotations.annotations) {
            EvidenceClaim claim = fromAnnotation(protein.id, ann)
            if (claim != null) out << claim
        }
        out
    }

    /** Extract all claims from all proteins in a genome. */
    List<EvidenceClaim> fromGenome(Genome genome) {
        List<EvidenceClaim> out = []
        for (Protein p : genome.proteins) {
            out.addAll(fromProtein(p))
        }
        out
    }

    /**
     * Group claims by {@link EvidenceClaim#functionKey()}. This is the
     * grouping used by {@link EvidenceCombiner} to combine multiple claims
     * about the same (protein, function) assignment.
     */
    static Map<String, List<EvidenceClaim>> groupByFunctionKey(List<EvidenceClaim> claims) {
        Map<String, List<EvidenceClaim>> out = new LinkedHashMap<>()
        for (EvidenceClaim c : claims) {
            out.computeIfAbsent(c.functionKey(), { k -> new ArrayList<EvidenceClaim>() }) << c
        }
        out
    }
}
