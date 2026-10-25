package gspa.config

/** Biological workflow defaults shared by annotate and integrate. Never resource paths. */
class WorkflowPresets {
    static final Map<String, String> DESCRIPTIONS = [
        prokaryote: 'One prokaryotic genome, domain unspecified; use bacteria/archaea for essentiality.',
        bacteria: 'One bacterial genome: bacterial essentials and taxon constraints.',
        archaea: 'One archaeal genome: archaeal essentials and taxon constraints.',
        eukaryote: 'One eukaryotic proteome or existing gene models; no operon assumptions.',
        mag: 'One incomplete prokaryotic bin: meta gene calling, no completeness-driven repair.',
        metagenome: 'Mixed community: meta gene calling, no pooled single-organism constraints.'
    ].asImmutable()

    private static final Map<String, String> ALIASES = [
        prokaryotic: 'prokaryote', bacterial: 'bacteria', archaeal: 'archaea',
        eukaryotic: 'eukaryote', metagenomics: 'metagenome'
    ].asImmutable()

    static String canonicalName(String name) {
        String key = name?.trim()?.toLowerCase(Locale.ROOT)
        key = ALIASES[key] ?: key
        if (!DESCRIPTIONS.containsKey(key)) {
            throw new IllegalArgumentException("Unknown preset '${name}'. Choose: ${DESCRIPTIONS.keySet().join(', ')}. Run 'gspa presets' for guidance.")
        }
        key
    }

    static Map annotationDefaults(String name, String kingdom = null) {
        String preset = canonicalName(name)
        String domain = kingdom ?: (preset in ['bacteria', 'archaea', 'eukaryote'] ? preset : 'auto')
        boolean mixed = preset == 'metagenome'
        boolean incomplete = preset == 'mag'
        boolean known = domain in ['bacteria', 'archaea', 'eukaryote']
        boolean eukaryote = domain == 'eukaryote'
        [
            input: [kingdom: domain, type: mixed ? 'metagenome' : incomplete ? 'mag' : 'genome',
                    requireExistingGenes: eukaryote],
            predictors: [
                // Start with domain annotation and a configured neural backbone. Optional
                // database searches/localization/reconstruction must be requested explicitly.
                similarity: [enabled: false],
                domains: [enabled: true],
                operons: [enabled: !eukaryote && !mixed],
                pathway: [enabled: false, crossfeeding: false],
                localization: [enabled: false]
            ],
            quality: [
                enabled: !mixed, scope: 'genome',
                completeness: [profile: mixed || !known ? 'skip' : domain, enforce: false],
                coherence: [enforce: false],
                consistency: [taxonConstraints: !mixed, enforce: known && !mixed,
                              organismTaxon: known && !mixed ? domain : null]
            ]
        ]
    }

    /** Keys are CLI option names; the CLI applies only options the user did not supply. */
    static Map<String, Object> integrationDefaults(String name) {
        String preset = canonicalName(name)
        boolean mixed = preset == 'metagenome'
        boolean incomplete = preset == 'mag'
        boolean known = preset in ['bacteria', 'archaea', 'eukaryote']
        List<String> priors = []
        if (!mixed) {
            if (known) priors << 'essentiality'
            if (!incomplete) priors << 'coherence'
            priors.addAll(['consistency', 'anatomical'])
            if (!incomplete && preset != 'eukaryote') priors << 'gap_filling'
            if (preset != 'eukaryote') priors << 'genomic_context'
        }
        [
            '--essential-profile': known ? preset : null,
            '--taxon': known ? preset : null,
            '--enable-priors': priors.join(','),
            '--enforce-stage2': !mixed && !incomplete,
            '--enforce-consistency-mode': 'relax',
            '--lite': true
        ]
    }
}
