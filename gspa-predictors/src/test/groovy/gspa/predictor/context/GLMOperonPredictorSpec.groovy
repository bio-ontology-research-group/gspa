package gspa.predictor.context

import gspa.model.*
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Files
import java.nio.file.Path

/**
 * Unit tests for {@link GLMOperonPredictor}. We don't run the real Python
 * sidecar in CI — we substitute a tiny shell stub that writes canned
 * TSV/NPZ output and assert the wrapper parses it and routes BP transfer
 * through the inherited {@link OperonPredictor#transferAnnotations}.
 */
class GLMOperonPredictorSpec extends Specification {

    @TempDir
    Path tmp

    def "wrapper parses operons.tsv, drops sub-threshold operons, and transfers BP terms"() {
        given:
        // Build a 4-protein genome where the gLM sidecar (mocked here)
        // returns two operons; the second is below confidence threshold
        // and must be dropped before BP transfer.
        Genome genome = buildGenome([
            [id: 'p1', start: 1,    end: 300,  strand: Strand.PLUS],
            [id: 'p2', start: 320,  end: 600,  strand: Strand.PLUS],
            [id: 'p3', start: 700,  end: 1000, strand: Strand.PLUS],
            [id: 'p4', start: 1010, end: 1300, strand: Strand.PLUS],
        ])
        // p1 has a BP term that should propagate inside operon (p1, p2).
        genome.findProtein('p1').annotations.add(
            new Annotation(type: AnnotationType.GO, value: 'GO:0006412',
                source: 'diamond', goAspect: 'BP'),
        )

        File runDir = tmp.toFile()
        File stub = writeShellStub(runDir, """\
#!/usr/bin/env bash
# Mock gLM sidecar — emits canned outputs based on the --operons-out path.
set -e
operons_out=""
conf_out=""
cents_out=""
embs_out=""
while [[ \$# -gt 0 ]]; do
  case "\$1" in
    --operons-out) operons_out="\$2"; shift 2;;
    --confidence-out) conf_out="\$2"; shift 2;;
    --centroids-out) cents_out="\$2"; shift 2;;
    --protein-embeddings-out) embs_out="\$2"; shift 2;;
    *) shift;;
  esac
done
printf 'p1\\tp2\\n' > "\$operons_out"
printf 'p3\\tp4\\n' >> "\$operons_out"
printf 'operon_idx\\tsize\\tconfidence\\n' > "\$conf_out"
printf 'op0\\t2\\t0.90\\n' >> "\$conf_out"
printf 'op1\\t2\\t0.20\\n' >> "\$conf_out"
# Empty NPZ files; wrapper does not parse them in step 1.
: > "\$cents_out"
: > "\$embs_out"
""")

        GLMOperonPredictor predictor = new GLMOperonPredictor(
            sidecarPath: stub.absolutePath,
            pythonExecutable: 'bash',
            mode: 'mock',
            workDir: runDir,
            minOperonConfidence: 0.5d,
        )

        when:
        Map<String, List<Annotation>> transfers = predictor.predictGenome(genome)

        then: 'high-confidence operon (p1, p2) survives'
        transfers['p2']?.any { it.value == 'GO:0006412' && it.goAspect == 'BP' }
        transfers['p2'][0].source == 'glm-operon'      // overridden getName() flows through
        transfers['p2'][0].evidence == 'IGC'
        transfers['p2'][0].score == 0.4

        and: 'low-confidence operon (p3, p4) is dropped — no transfers there'
        !(transfers['p3']?.any { it.value == 'GO:0006412' })
        !(transfers['p4']?.any { it.value == 'GO:0006412' })
    }

    def "isAvailable rejects missing sidecar"() {
        expect:
        !new GLMOperonPredictor(sidecarPath: null).isAvailable()
        !new GLMOperonPredictor(sidecarPath: '/no/such/file').isAvailable()
    }

    def "isAvailable rejects real-mode without weights dir"() {
        given:
        File stub = writeShellStub(tmp.toFile(), '#!/usr/bin/env bash\nexit 0\n')

        expect:
        !new GLMOperonPredictor(
            sidecarPath: stub.absolutePath,
            mode: 'real',
            weightsDir: null,
            pythonExecutable: 'bash',
        ).isAvailable()
    }

    // ------------------------------------------------------------ helpers

    private File writeShellStub(File dir, String body) {
        File f = new File(dir, "stub_${System.nanoTime()}.sh")
        f.text = body
        f.setExecutable(true)
        f
    }

    private Genome buildGenome(List<Map> geneSpecs) {
        Genome genome = new Genome(id: 'test')
        Contig contig = new Contig(id: 'c1', sequence: 'A' * 2000)
        geneSpecs.each { Map spec ->
            Feature feature = new Feature(
                seqId: 'c1', type: FeatureType.CDS,
                start: spec.start, end: spec.end,
                strand: spec.strand,
            )
            feature.id = spec.id
            contig.addFeature(feature)
            Protein protein = new Protein(
                id: spec.id, sequence: 'M' * ((spec.end - spec.start) / 3 as int),
                sourceFeature: feature,
            )
            contig.addProtein(protein)
        }
        contig.sortFeatures()
        genome.addContig(contig)
        genome
    }
}
