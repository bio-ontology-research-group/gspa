package gspa.cli

import picocli.CommandLine.Command
import picocli.CommandLine.Option

import java.nio.file.Files
import java.nio.file.StandardCopyOption

/**
 * {@code gspa visualize} — emit a self-contained HTML browser for a GSPA
 * workspace.
 *
 * Reads the conventional layout produced by {@code gspa annotate} +
 * {@code gspa integrate}:
 *
 * <pre>
 *   workdir/
 *     prokka_out/{genome}.{tsv,gff}                  per-CDS metadata + coords
 *     gspa_out/                                       integrator outputs
 *       integrated.tsv                                posteriors
 *       provenance.json                               supporting predictors
 *       operons.tsv                                   GSPA OperonPredictor
 *       quality_gspa.json                             GAEF report (with detail)
 *     mdf_out/{genome}.mdf.tsv                        mDeepFRI sidecar (optional)
 *     proteinfer_out/{genome}.proteinfer.tsv          ProteInfer sidecar (optional)
 *     clean_out/{genome}.clean.tsv                    CLEAN sidecar (optional)
 *     amrfinder_out/{genome}.amr.tsv                  AmrFinderPlus (optional)
 *     antismash_out/regions.js                        antiSMASH BGCs (optional)
 *     ipr_out/{genome}.faa.tsv                        InterProScan (optional)
 *     input/{genome}_assembly.fa                      reference for igv.js (optional)
 * </pre>
 *
 * Implementation: extracts a Python templating script bundled as a JAR
 * resource, then invokes {@code python3} with the workspace path passed via
 * env vars. Python is required (any version 3.8+); the script has no
 * external dependencies beyond the standard library.
 *
 * The Python sibling lives at {@code visualize/make_viz.py}. Keeping it as
 * a Python script (vs reimplementing in Groovy) lets the same tool be
 * invoked from {@code gspa-nf} pipelines too.
 */
@Command(
    name = 'visualize',
    description = 'Emit a single self-contained HTML browser for a GSPA workspace.',
    mixinStandardHelpOptions = true
)
class VisualizeCommand implements Runnable {

    @Option(names = ['-w', '--workdir'], required = true,
            description = 'GSPA workspace directory (contains prokka_out/, gspa_out/, etc.)')
    File workdir

    @Option(names = ['-o', '--out'],
            description = 'HTML output path. Default: <workdir>/<genome>_browser.html')
    File out

    @Option(names = ['-r', '--run-dir'],
            description = 'Subdir of <workdir> with integrator output (integrated.tsv, ' +
                          'provenance.json, operons.tsv, quality_gspa.json). Default: gspa_out.')
    String runDir = 'gspa_out'

    @Option(names = ['-g', '--genome-id'],
            description = 'Genome id used to resolve per-tool sidecar filenames ' +
                          '(<id>.mdf.tsv, <id>.proteinfer.tsv, etc.). Default: workdir basename.')
    String genomeId

    @Option(names = ['--go-obo'],
            description = 'Path to go.obo for GO term name lookup. Default: ' +
                          '$GSPA_REF/go.obo or /data/hohndor/gspa/reference/go.obo.')
    File goObo

    @Option(names = ['--ec2go'],
            description = 'Path to ec2go.txt for EC name lookup. Default: ' +
                          '$GSPA_REF/ec2go.txt or /data/hohndor/gspa/reference/ec2go.txt.')
    File ec2goFile

    @Option(names = ['--fasta'],
            description = 'FASTA reference (required by igv.js for the genome browser tab; ' +
                          'embedded into the HTML as a base64 data URL). ' +
                          'Default: <workdir>/input/<genome>_assembly.fa.')
    File fasta

    @Option(names = ['--python'],
            description = 'Python interpreter to invoke. Default: python3.')
    String python = 'python3'

    @Override
    void run() {
        if (!workdir.isDirectory()) {
            throw new IllegalArgumentException("Workdir does not exist: ${workdir}")
        }
        String gid = genomeId ?: workdir.name
        File htmlOut = out ?: new File(new File(workdir, runDir), "${gid}_browser.html")

        println "GSPA visualize"
        println "  Workdir:   ${workdir}"
        println "  Run dir:   ${runDir}"
        println "  Genome id: ${gid}"
        println "  Output:    ${htmlOut}"

        // Extract the bundled Python script to a temp file.
        File scriptDir = Files.createTempDirectory('gspa-visualize-').toFile()
        scriptDir.deleteOnExit()
        File makeViz = extractResource('visualize/make_viz.py', new File(scriptDir, 'make_viz.py'))
        // Bundle predict_operons.py too (handy for callers but not invoked here).
        extractResource('visualize/predict_operons.py', new File(scriptDir, 'predict_operons.py'))

        ProcessBuilder pb = new ProcessBuilder(python, makeViz.absolutePath)
        pb.redirectErrorStream(true)
        pb.inheritIO()
        Map<String, String> env = pb.environment()
        env['GSPA_WORKDIR']  = workdir.absolutePath
        env['GSPA_RUN_DIR']  = runDir
        env['GSPA_GENOME_ID']= gid
        env['GSPA_OUT']      = htmlOut.absolutePath
        if (goObo != null)    env['GSPA_GO_OBO'] = goObo.absolutePath
        if (ec2goFile != null)env['GSPA_EC2GO']  = ec2goFile.absolutePath
        if (fasta != null)    env['GSPA_FASTA']  = fasta.absolutePath

        Process p = pb.start()
        int code = p.waitFor()
        if (code != 0) {
            throw new RuntimeException("make_viz.py exited with code ${code}")
        }
        if (htmlOut.exists()) {
            println "  OK: ${htmlOut} (${htmlOut.length() / 1_000_000} MB)"
        }
    }

    /**
     * Copy a bundled JAR resource to a temp file and return the file handle.
     * Used to make the bundled Python templater visible on disk so an
     * external interpreter can read it.
     */
    private File extractResource(String resourcePath, File dest) {
        InputStream is = getClass().classLoader.getResourceAsStream(resourcePath)
        if (is == null) {
            throw new IllegalStateException("Bundled resource missing: ${resourcePath}")
        }
        try {
            Files.copy(is, dest.toPath(), StandardCopyOption.REPLACE_EXISTING)
        } finally {
            is.close()
        }
        dest.setReadable(true)
        return dest
    }
}
