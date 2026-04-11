package gspa.model

import groovy.transform.Canonical
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * Taxonomy information for a genome, supporting both NCBI and GTDB lineages.
 */
@Canonical
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class TaxonomyInfo {

    /** NCBI Taxonomy ID */
    String ncbiTaxId

    /** GTDB lineage string (e.g., "d__Bacteria;p__Proteobacteria;...") */
    String gtdbLineage

    /** Scientific name */
    String scientificName

    /** Organism domain */
    OrganismDomain domain

    /** Lineage ranks: domain, phylum, class, order, family, genus, species */
    Map<String, String> lineage = [:]

    /** Parse GTDB lineage string into ranks */
    static TaxonomyInfo fromGtdb(String gtdbLineage) {
        def info = new TaxonomyInfo(gtdbLineage: gtdbLineage)
        def prefixToRank = [
            'd__': 'domain', 'p__': 'phylum', 'c__': 'class',
            'o__': 'order', 'f__': 'family', 'g__': 'genus', 's__': 'species'
        ]
        gtdbLineage.split(';').each { part ->
            def prefix = part.take(3)
            def rank = prefixToRank[prefix]
            if (rank) {
                info.lineage[rank] = part.substring(3)
            }
        }
        def domainName = info.lineage['domain']
        if (domainName) {
            info.domain = OrganismDomain.fromName(domainName)
        }
        info
    }
}
