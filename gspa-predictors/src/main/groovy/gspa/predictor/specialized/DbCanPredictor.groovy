package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.predictor.AbstractToolPredictor

/**
 * dbCAN predictor for carbohydrate-active enzyme (CAZyme) annotation.
 * Uses HMMER, DIAMOND, and eCAMI for CAZyme family identification.
 */
class DbCanPredictor extends AbstractToolPredictor {

    /** dbCAN database directory */
    String dbDir

    int threads = Runtime.runtime.availableProcessors()

    @Override
    String getName() { 'dbcan' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.CAZYME] as Set }

    @Override
    String getExecutable() { 'run_dbcan' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def cmd = [
            executablePath ?: executable,
            inputFasta.absolutePath,
            'protein',
            '--out_dir', outputDir.absolutePath,
        ]
        if (dbDir) cmd.addAll(['--db_dir', dbDir])
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def overviewFile = new File(outputDir, 'overview.txt')
        if (!overviewFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        overviewFile.eachLine { line ->
            if (line.startsWith('Gene ID') || line.startsWith('#')) return
            def fields = line.split('\t')
            if (fields.length < 6) return

            String proteinId = fields[0]
            String hmmer = fields[1]       // HMMER CAZy family
            String diamond = fields[3]     // DIAMOND CAZy family
            int toolCount = fields[5] as int  // number of tools agreeing

            // Use consensus
            String cazyFamily = hmmer != '-' ? hmmer : (diamond != '-' ? diamond : null)
            if (!cazyFamily) return

            double score = toolCount >= 3 ? 0.95 : (toolCount >= 2 ? 0.8 : 0.6)

            // Parse multiple families (can be e.g., "GH1+GH2")
            cazyFamily.split(/\+/).each { family ->
                family = family.replaceAll(/\(.*\)/, '').trim()
                if (family && family != '-') {
                    results[proteinId] << new Annotation(
                        type: AnnotationType.CAZYME,
                        value: family,
                        source: name,
                        score: score,
                        metadata: [toolsAgreeing: toolCount]
                    )
                }
            }
        }
        results
    }
}
