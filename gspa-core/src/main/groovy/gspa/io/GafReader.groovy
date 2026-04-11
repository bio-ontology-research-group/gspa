package gspa.io

import gspa.model.Annotation
import gspa.model.AnnotationType

/**
 * Reads Gene Association File (GAF) format.
 * GAF 2.2 is the standard format for GO annotations.
 * See: http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
 */
class GafReader {

    /**
     * Parse a GAF file and return a map of protein ID -> list of GO annotations.
     * Uses DB_Object_ID (column 2) as the protein identifier by default.
     */
    static Map<String, List<Annotation>> readGaf(File gafFile, String source = 'GOA') {
        Map<String, List<Annotation>> result = [:].withDefault { [] }

        def reader = gafFile.name.endsWith('.gz') ?
            new BufferedReader(new InputStreamReader(new java.util.zip.GZIPInputStream(new FileInputStream(gafFile)))) :
            gafFile.newReader()

        try {
            reader.eachLine { line ->
                if (line.startsWith('!') || line.trim().isEmpty()) return

                def fields = line.split('\t')
                if (fields.length < 15) return

                String proteinId = fields[1]    // DB_Object_ID
                String qualifier = fields[3]    // Qualifier (NOT, contributes_to, etc.)
                String goId = fields[4]         // GO ID
                String evidenceCode = fields[6] // Evidence code
                String aspect = fields[8]       // Aspect (F=MF, P=BP, C=CC)

                // Skip NOT qualifiers
                if (qualifier.toUpperCase().contains('NOT')) return

                String goAspect = switch (aspect) {
                    case 'F' -> 'MF'
                    case 'P' -> 'BP'
                    case 'C' -> 'CC'
                    default -> aspect
                }

                def annotation = new Annotation(
                    type: AnnotationType.GO,
                    value: goId,
                    source: source,
                    evidence: evidenceCode,
                    goAspect: goAspect,
                    score: evidenceCodeScore(evidenceCode)
                )

                result[proteinId] << annotation
            }
        } finally {
            reader.close()
        }
        result
    }

    /**
     * Rough confidence score based on GO evidence code.
     * Experimental evidence gets higher scores.
     */
    static double evidenceCodeScore(String code) {
        switch (code) {
            // Experimental
            case ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP']:
                return 1.0
            // Phylogenetic
            case ['IBA', 'IBD', 'IKR', 'IRD']:
                return 0.8
            // Computational
            case ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA']:
                return 0.6
            // Author/curator
            case ['TAS', 'NAS']:
                return 0.7
            // Electronic
            case 'IEA':
                return 0.4
            default:
                return 0.5
        }
    }
}
