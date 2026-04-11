package gspa.config

import gspa.model.OrganismDomain
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Manages essential function profiles for completeness evaluation.
 * Each profile defines a set of GO terms that are expected to be present
 * in a genome of a given organism type.
 *
 * Profiles can be loaded from built-in presets or custom TSV files,
 * and modified at runtime (add/remove terms).
 */
class EssentialFunctions {

    private static final Logger log = LoggerFactory.getLogger(EssentialFunctions)

    /** Profile name */
    String profileName

    /** GO term -> category (e.g., "Core", "Glucose Metabolism") */
    Map<String, String> functions = [:]

    /** GO term -> human-readable description */
    Map<String, String> descriptions = [:]

    /**
     * Load from a TSV file.
     * Format: GO_term\tcategory\tdescription
     */
    static EssentialFunctions loadFromTsv(File tsvFile, String profileName = null) {
        def ef = new EssentialFunctions(profileName: profileName ?: tsvFile.name.replaceAll(/\.tsv$/, ''))
        tsvFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t')
            if (fields.length >= 2) {
                String goTerm = fields[0].trim()
                String category = fields[1].trim()
                String desc = fields.length >= 3 ? fields[2].trim() : ''
                ef.functions[goTerm] = category
                ef.descriptions[goTerm] = desc
            }
        }
        log.info("Loaded ${ef.functions.size()} essential functions from ${tsvFile} (profile: ${ef.profileName})")
        ef
    }

    /**
     * Load from a built-in resource.
     */
    static EssentialFunctions loadPreset(String profileName) {
        def stream = EssentialFunctions.classLoader.getResourceAsStream(
            "essential-functions/${profileName}.tsv")
        if (stream == null) {
            log.warn("No built-in essential function profile for: ${profileName}")
            return getDefault(OrganismDomain.fromName(profileName))
        }
        def tmpFile = File.createTempFile("ef_${profileName}_", '.tsv')
        try {
            tmpFile.text = stream.text
            return loadFromTsv(tmpFile, profileName)
        } finally {
            tmpFile.delete()
        }
    }

    /**
     * Get the default essential function set for a domain.
     * Tries to load from built-in TSV resource first, falls back to hardcoded.
     */
    static EssentialFunctions getDefault(OrganismDomain domain) {
        String profileName = switch (domain) {
            case OrganismDomain.BACTERIA -> 'bacteria'
            case OrganismDomain.ARCHAEA -> 'archaea'
            case OrganismDomain.EUKARYA -> 'eukaryote'
            default -> 'bacteria'
        }

        // Try loading from resource TSV
        def stream = EssentialFunctions.classLoader.getResourceAsStream(
            "essential-functions/${profileName}.tsv")
        if (stream != null) {
            def tmpFile = File.createTempFile("ef_${profileName}_", '.tsv')
            try {
                tmpFile.text = stream.text
                return loadFromTsv(tmpFile, profileName)
            } finally {
                tmpFile.delete()
            }
        }

        // Fallback to hardcoded
        def ef = new EssentialFunctions()
        ef.profileName = profileName
        switch (domain) {
            case OrganismDomain.BACTERIA:
                addBacteriaEssentials(ef)
                break
            case OrganismDomain.ARCHAEA:
                addArchaeaEssentials(ef)
                break
            case OrganismDomain.EUKARYA:
                addEukaryoteEssentials(ef)
                break
            default:
                addBacteriaEssentials(ef)
        }
        log.info("Using hardcoded essential functions for ${ef.profileName}: ${ef.functions.size()} terms")
        ef
    }

    /**
     * Apply runtime modifications (add/remove terms) from config.
     */
    EssentialFunctions withModifications(List<String> addTerms, List<String> removeTerms) {
        def modified = new EssentialFunctions(
            profileName: this.profileName,
            functions: new LinkedHashMap(this.functions),
            descriptions: new LinkedHashMap(this.descriptions)
        )
        addTerms?.each { term ->
            modified.functions[term] = 'Custom'
            modified.descriptions[term] = 'User-added essential function'
        }
        removeTerms?.each { term ->
            modified.functions.remove(term)
            modified.descriptions.remove(term)
        }
        modified
    }

    /** Get all GO term IDs in this profile */
    Set<String> getGoTerms() {
        functions.keySet()
    }

    /** Get GO terms in a specific category */
    Set<String> getTermsByCategory(String category) {
        functions.findAll { it.value == category }.keySet()
    }

    /** Get all categories */
    Set<String> getCategories() {
        functions.values() as Set
    }

    // --- Hardcoded essential functions from the paper ---
    // Based on Table S1: Syn1.0 essential genes mapped to GO terms

    private static void addBacteriaEssentials(EssentialFunctions ef) {
        // Core functions (from paper's Table 1 / Syn1.0)
        def core = [
            'GO:0006259': 'DNA metabolic process',
            'GO:0006260': 'DNA replication',
            'GO:0006281': 'DNA repair',
            'GO:0006351': 'DNA-templated transcription',
            'GO:0006412': 'translation',
            'GO:0006418': 'tRNA aminoacylation for protein translation',
            'GO:0006457': 'protein folding',
            'GO:0051301': 'cell division',
            'GO:0009058': 'biosynthetic process',
            'GO:0006629': 'lipid metabolic process',
            'GO:0009059': 'macromolecule biosynthetic process',
            'GO:0006520': 'amino acid metabolic process',
            'GO:0009117': 'nucleotide metabolic process',
            'GO:0006091': 'generation of precursor metabolites and energy',
            'GO:0055085': 'transmembrane transport',
            'GO:0006810': 'transport',
            'GO:0005975': 'carbohydrate metabolic process',
            'GO:0016192': 'vesicle-mediated transport',
            'GO:0006399': 'tRNA metabolic process',
            'GO:0006396': 'RNA processing',
            'GO:0009306': 'protein secretion',
            'GO:0006886': 'intracellular protein transport',
        ]
        core.each { goId, desc ->
            ef.functions[goId] = 'Core'
            ef.descriptions[goId] = desc
        }
    }

    private static void addArchaeaEssentials(EssentialFunctions ef) {
        // Start with core functions shared with bacteria
        addBacteriaEssentials(ef)
        ef.profileName = 'archaea'

        // Add archaeal-specific
        ef.functions['GO:0015948'] = 'Archaeal-specific'
        ef.descriptions['GO:0015948'] = 'methanogenesis'
        ef.functions['GO:0043571'] = 'Archaeal-specific'
        ef.descriptions['GO:0043571'] = 'maintenance of CRISPR repeat elements'
    }

    private static void addEukaryoteEssentials(EssentialFunctions ef) {
        // Core functions shared across all life
        addBacteriaEssentials(ef)
        ef.profileName = 'eukaryote'

        // Eukaryote-specific
        def eukSpecific = [
            'GO:0006606': 'protein import into nucleus',
            'GO:0006913': 'nucleocytoplasmic transport',
            'GO:0016192': 'vesicle-mediated transport',
            'GO:0007049': 'cell cycle',
            'GO:0000278': 'mitotic cell cycle',
            'GO:0006914': 'autophagy',
            'GO:0006888': 'endoplasmic reticulum to Golgi vesicle-mediated transport',
            'GO:0007005': 'mitochondrion organization',
            'GO:0006119': 'oxidative phosphorylation',
            'GO:0008380': 'RNA splicing',
        ]
        eukSpecific.each { goId, desc ->
            ef.functions[goId] = 'Eukaryote-specific'
            ef.descriptions[goId] = desc
        }

        // Remove prokaryote-specific terms that don't apply
        ef.functions.remove('GO:0009306') // protein secretion (prokaryotic sense)
    }
}
