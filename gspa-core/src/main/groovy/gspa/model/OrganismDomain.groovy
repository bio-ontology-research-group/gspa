package gspa.model

/**
 * Biological domain/kingdom of the organism.
 * Used for selecting appropriate predictors and essential function profiles.
 */
enum OrganismDomain {
    BACTERIA,
    ARCHAEA,
    EUKARYA,
    VIRUS,
    UNKNOWN

    static OrganismDomain fromName(String name) {
        switch (name?.toLowerCase()) {
            case 'bacteria': return BACTERIA
            case 'archaea': return ARCHAEA
            case ['eukarya', 'eukaryota', 'eukaryote']: return EUKARYA
            case ['virus', 'viruses']: return VIRUS
            default: return UNKNOWN
        }
    }

    boolean isProkaryote() {
        this == BACTERIA || this == ARCHAEA
    }

    boolean isEukaryote() {
        this == EUKARYA
    }
}
