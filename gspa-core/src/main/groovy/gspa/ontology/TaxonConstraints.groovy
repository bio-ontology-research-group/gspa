package gspa.ontology

import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Parses and stores taxon constraints for GO terms.
 * Sources: go-computed-taxon-constraints.obo or GO-Plus OWL axioms.
 *
 * Taxon constraints define:
 * - only_in_taxon (RO:0002160): GO term can only occur in this taxon
 * - never_in_taxon (RO:0002161): GO term never occurs in this taxon
 */
class TaxonConstraints {

    private static final Logger log = LoggerFactory.getLogger(TaxonConstraints)

    /** GO term -> set of taxa where the term can only occur */
    Map<String, Set<String>> onlyInTaxon = [:].withDefault { [] as Set }

    /** GO term -> set of taxa where the term never occurs */
    Map<String, Set<String>> neverInTaxon = [:].withDefault { [] as Set }

    /**
     * Load taxon constraints from a GO-Plus OWL ontology.
     * Extracts only_in_taxon and never_in_taxon axioms.
     */
    void loadFromGoOntology(GoOntology go) {
        log.info("Extracting taxon constraints from GO ontology...")

        // Extract only_in_taxon relationships
        def onlyInAxioms = go.getExistentialAxioms(GoOntology.ONLY_IN_TAXON)
        onlyInAxioms.each { axiom ->
            onlyInTaxon[axiom.subject] << axiom.object
        }

        // Extract never_in_taxon relationships
        def neverInAxioms = go.getExistentialAxioms(GoOntology.NEVER_IN_TAXON)
        neverInAxioms.each { axiom ->
            neverInTaxon[axiom.subject] << axiom.object
        }

        log.info("Loaded ${onlyInTaxon.size()} only_in_taxon and ${neverInTaxon.size()} never_in_taxon constraints")
    }

    /**
     * Load from the GO-computed taxon constraints OBO file
     * (go-computed-taxon-constraints.obo).
     *
     * Format: blocks of [Term] with id: GO:NNNNNNN and
     * relationship: RO:0002162 NCBITaxon:NNN  (only_in_taxon)
     * property_value: RO:0002161 NCBITaxon:NNN (never_in_taxon)
     */
    void loadFromObo(File oboFile) {
        log.info("Loading taxon constraints from OBO: ${oboFile}")
        String currentId = null
        oboFile.eachLine { line ->
            line = line.trim()
            if (line == '[Term]') {
                currentId = null
            } else if (line.startsWith('id: GO:')) {
                currentId = line.substring(4).trim()
            } else if (currentId != null && line.startsWith('relationship: RO:0002162 ')) {
                // only_in_taxon (RO:0002162 = in_taxon used as only_in via relationship)
                String taxon = line.substring('relationship: RO:0002162 '.length()).trim()
                onlyInTaxon[currentId] << taxon
            } else if (currentId != null && line.startsWith('property_value: RO:0002161 ')) {
                // never_in_taxon
                String taxon = line.substring('property_value: RO:0002161 '.length()).trim()
                // OBO property_value may have a trailing xsd:string — strip it
                if (taxon.contains(' ')) taxon = taxon.split('\\s+')[0]
                neverInTaxon[currentId] << taxon
            }
        }
        log.info("Loaded from OBO: ${onlyInTaxon.size()} only_in_taxon, ${neverInTaxon.size()} never_in_taxon constraints")
    }

    /**
     * Load taxon constraints from a TSV file.
     * Expected format: GO_term\trelation\tNCBITaxon_ID
     */
    void loadFromTsv(File tsvFile) {
        log.info("Loading taxon constraints from: ${tsvFile}")

        tsvFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t')
            if (fields.length < 3) return

            String goTerm = fields[0].trim()
            String relation = fields[1].trim()
            String taxon = fields[2].trim()

            switch (relation) {
                case ['only_in_taxon', 'RO:0002160']:
                    onlyInTaxon[goTerm] << taxon
                    break
                case ['never_in_taxon', 'RO:0002161']:
                    neverInTaxon[goTerm] << taxon
                    break
            }
        }

        log.info("Loaded ${onlyInTaxon.size()} only_in_taxon and ${neverInTaxon.size()} never_in_taxon constraints")
    }

    /**
     * Get all taxon constraints triggered by a set of GO annotations.
     * Returns the aggregated positive and negative taxon requirements.
     */
    TaxonRequirements getRequirements(Set<String> goTerms) {
        Set<String> positiveTaxa = []
        Set<String> negativeTaxa = []
        Map<String, Set<String>> positiveSource = [:].withDefault { [] as Set }
        Map<String, Set<String>> negativeSource = [:].withDefault { [] as Set }

        goTerms.each { goTerm ->
            onlyInTaxon[goTerm]?.each { taxon ->
                positiveTaxa << taxon
                positiveSource[taxon] << goTerm
            }
            neverInTaxon[goTerm]?.each { taxon ->
                negativeTaxa << taxon
                negativeSource[taxon] << goTerm
            }
        }

        new TaxonRequirements(
            positiveTaxa: positiveTaxa,
            negativeTaxa: negativeTaxa,
            positiveSource: positiveSource,
            negativeSource: negativeSource
        )
    }

    /** Check if a GO term has any taxon constraints */
    boolean hasConstraints(String goTerm) {
        onlyInTaxon.containsKey(goTerm) || neverInTaxon.containsKey(goTerm)
    }

    /** Total number of constrained GO terms */
    int constrainedTermCount() {
        (onlyInTaxon.keySet() + neverInTaxon.keySet()).size()
    }
}

/**
 * Aggregated taxon requirements from a set of GO annotations.
 */
class TaxonRequirements {
    /** Taxa the organism must belong to */
    Set<String> positiveTaxa = []
    /** Taxa the organism must NOT belong to */
    Set<String> negativeTaxa = []
    /** Which GO terms required each positive taxon */
    Map<String, Set<String>> positiveSource = [:]
    /** Which GO terms required each negative taxon */
    Map<String, Set<String>> negativeSource = [:]
}
