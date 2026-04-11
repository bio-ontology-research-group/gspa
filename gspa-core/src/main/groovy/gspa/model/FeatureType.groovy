package gspa.model

enum FeatureType {
    CDS('CDS'),
    GENE('gene'),
    MRNA('mRNA'),
    TRNA('tRNA'),
    RRNA('rRNA'),
    REPEAT_REGION('repeat_region'),
    CRISPR('CRISPR'),
    SIGNAL_PEPTIDE('signal_peptide'),
    TRANSMEMBRANE('transmembrane_region'),
    REGION('region'),
    MISC_FEATURE('misc_feature')

    final String gff3Type

    FeatureType(String gff3Type) {
        this.gff3Type = gff3Type
    }

    static FeatureType fromGff3(String type) {
        values().find { it.gff3Type.equalsIgnoreCase(type) } ?: MISC_FEATURE
    }
}
