package gspa.predictor.function

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * eggNOG-mapper predictor for ortholog group assignment.
 * Provides COG, KEGG, GO annotations via orthology inference.
 */
class EggNogMapperPredictor extends AbstractToolPredictor {

    /** eggNOG database directory */
    String dataDir

    /** Search mode: diamond, mmseqs, hmmer */
    String searchMode = 'diamond'

    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'eggnog-mapper' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.GO, AnnotationType.EC, AnnotationType.COG, AnnotationType.KEGG] as Set
    }

    @Override
    String getExecutable() { 'emapper.py' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def prefix = new File(outputDir, 'eggnog_results')
        def cmd = [
            executablePath ?: executable,
            '-i', inputFasta.absolutePath,
            '-o', prefix.absolutePath,
            '--cpu', threads.toString(),
            '-m', searchMode,
            '--itype', 'proteins',
            '--go_evidence', 'non-electronic',
        ]
        if (dataDir) {
            cmd.addAll(['--data_dir', dataDir])
        }
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def annotFile = new File(outputDir, 'eggnog_results.emapper.annotations')
        if (!annotFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        annotFile.eachLine { line ->
            if (line.startsWith('#')) return
            def fields = line.split('\t')
            if (fields.length < 21) return

            String proteinId = fields[0]
            String eggNogOGs = fields[4]       // Orthologous groups
            String goTerms = fields[9]         // GO terms
            String ecNumbers = fields[10]      // EC numbers
            String keggKOs = fields[11]        // KEGG KO
            String cogCats = fields[6]         // COG categories
            double score = safeDouble(fields[3]) // seed evalue

            double confScore = score > 0 ? Math.min(1.0, -Math.log10(score) / 30.0) : 0.7

            // GO annotations
            if (goTerms && goTerms != '-') {
                goTerms.split(',').each { goTerm ->
                    goTerm = goTerm.trim()
                    if (goTerm.startsWith('GO:')) {
                        results[proteinId] << new Annotation(
                            type: AnnotationType.GO, value: goTerm,
                            source: name, score: confScore, evidence: 'IEA',
                            metadata: [orthologGroup: eggNogOGs]
                        )
                    }
                }
            }

            // EC numbers
            if (ecNumbers && ecNumbers != '-') {
                ecNumbers.split(',').each { ec ->
                    ec = ec.trim()
                    if (ec) {
                        results[proteinId] << new Annotation(
                            type: AnnotationType.EC, value: "EC:${ec}",
                            source: name, score: confScore,
                            metadata: [orthologGroup: eggNogOGs]
                        )
                    }
                }
            }

            // KEGG KO
            if (keggKOs && keggKOs != '-') {
                keggKOs.split(',').each { ko ->
                    ko = ko.trim()
                    if (ko) {
                        results[proteinId] << new Annotation(
                            type: AnnotationType.KEGG, value: ko,
                            source: name, score: confScore,
                        )
                    }
                }
            }

            // COG categories
            if (cogCats && cogCats != '-' && cogCats != 'S') {
                results[proteinId] << new Annotation(
                    type: AnnotationType.COG, value: cogCats,
                    source: name, score: confScore,
                    metadata: [orthologGroup: eggNogOGs]
                )
            }
        }
        results
    }

    private static double safeDouble(String s) {
        try { return s as double } catch (Exception e) { return 0.0 }
    }
}
