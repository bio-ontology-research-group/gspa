package gspa.model

/**
 * Types of functional annotations that GSPA can handle.
 */
enum AnnotationType {
    GO,         // Gene Ontology
    EC,         // Enzyme Commission number
    KEGG,       // KEGG Orthology
    COG,        // Clusters of Orthologous Groups
    PFAM,       // Pfam domain
    TIGRFAM,    // TIGRFAM domain
    INTERPRO,   // InterPro entry
    CDD,        // Conserved Domain Database
    SUPERFAMILY,// SUPERFAMILY domain
    AMR,        // Antimicrobial resistance
    CAZYME,     // Carbohydrate-active enzyme
    VFDB,       // Virulence factor
    CRISPR,     // CRISPR-Cas system
    SIGNAL_PEPTIDE,
    TRANSMEMBRANE,
    SUBCELLULAR_LOCALIZATION,
    SECONDARY_METABOLITE,
    DISORDER,           // Intrinsically disordered region
    TARGETING_PEPTIDE,  // N-terminal mito/chloro/Sec/Tat targeting
    PTM_SITE,           // Post-translational modification site (1-residue)
    PPI_INTERFACE,      // Protein-protein interaction interface residue
    PROPHAGE,           // Genomic region: integrated bacteriophage
    PLASMID,            // Genomic region: extrachromosomal plasmid
    VIRAL_CONTIG,       // Whole-contig viral classification
    PHAGE_FUNCTION,     // Per-protein phage gene category (capsid, integrase, ...)
    CUSTOM              // User-defined annotation type
}
