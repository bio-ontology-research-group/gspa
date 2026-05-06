package gspa.predictor.cluster

import gspa.model.Genome
import gspa.model.Protein
import gspa.model.ProteinCluster
import gspa.model.ProteinClusterSet
import gspa.model.ProteinRef
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Single-genome implementation of {@link ProteinClusterer} using
 * MMseqs2 {@code linclust}.
 *
 * <p>Part 1 of Phase 10 uses this to cluster a genome's predicted
 * proteome at ~90% identity / 80% coverage so per-protein predictors
 * (DIAMOND, InterProScan, etc.) can run on representatives only; the
 * annotations are then propagated to members by
 * {@link ClusterAnnotationPropagator}.</p>
 *
 * <p>Emits singleton clusters for unique proteins. The representative
 * of each cluster is the one chosen by linclust (longest by default).</p>
 *
 * <p>Phase 11 will provide a cross-genome clusterer using the same
 * interface; this class intentionally accepts {@code List<Genome>} so
 * the wiring stays identical.</p>
 */
class IntragenomeClusterer implements ProteinClusterer {

    private static final Logger log = LoggerFactory.getLogger(IntragenomeClusterer)

    /** Absolute path to the mmseqs executable. Defaults to PATH lookup. */
    String mmseqsExecutable = 'mmseqs'

    /** Directory used for temporary FASTA / mmseqs DBs. Defaults to a tempdir per call. */
    File workDir

    /** Override the mmseqs invocation for testing (receives the command list, returns exit code). */
    Closure<Integer> commandRunner = { List<String> cmd ->
        def pb = new ProcessBuilder(cmd).redirectErrorStream(true)
        def p = pb.start()
        def stdout = p.inputStream.text
        int code = p.waitFor()
        if (code != 0) log.warn("mmseqs failed (exit ${code}):\n${stdout}")
        code
    }

    @Override
    ProteinClusterSet cluster(List<Genome> genomes, double identity, double coverage) {
        if (genomes == null || genomes.isEmpty()) {
            return new ProteinClusterSet(clusters: [], identityThreshold: identity, coverageThreshold: coverage)
        }
        if (genomes.size() > 1) {
            // Intragenome-only: defer multi-genome to Phase 11's CrossGenomeClusterer.
            throw new IllegalArgumentException(
                "IntragenomeClusterer only accepts one genome; use CrossGenomeClusterer for ${genomes.size()} genomes")
        }
        Genome g = genomes[0]
        List<Protein> proteins = g.proteins ?: []
        if (proteins.isEmpty()) {
            return new ProteinClusterSet(clusters: [], identityThreshold: identity, coverageThreshold: coverage)
        }

        File tmp = workDir ?: File.createTempDir('intragenome-cluster-')
        if (!tmp.exists()) tmp.mkdirs()

        try {
            File inputFasta = new File(tmp, 'input.faa')
            writeFasta(proteins, inputFasta)

            File resultPrefix = new File(tmp, 'result')
            File tmpDir = new File(tmp, 'mmseqs-tmp')
            tmpDir.mkdirs()

            List<String> cmd = [
                mmseqsExecutable, 'easy-linclust',
                inputFasta.absolutePath,
                resultPrefix.absolutePath,
                tmpDir.absolutePath,
                '--min-seq-id', String.valueOf(identity),
                '-c', String.valueOf(coverage),
                '--cov-mode', '0',
                '--threads', String.valueOf(Runtime.runtime.availableProcessors()),
            ]
            int code = commandRunner.call(cmd)
            if (code != 0) {
                log.warn("MMseqs2 linclust failed (exit ${code}); returning singleton clusters")
                return singletonClusters(g, proteins, identity, coverage)
            }

            File clusterTsv = new File(resultPrefix.absolutePath + '_cluster.tsv')
            if (!clusterTsv.exists()) {
                log.warn("linclust output ${clusterTsv} missing; returning singletons")
                return singletonClusters(g, proteins, identity, coverage)
            }
            parseClusterTsv(clusterTsv, g, identity, coverage)
        } finally {
            if (workDir == null) tmp.deleteDir()
        }
    }

    /** Parse MMseqs2 _cluster.tsv (rep_id TAB member_id, one row per member). */
    static ProteinClusterSet parseClusterTsv(File tsv, Genome g, double identity, double coverage) {
        Map<String, List<String>> repToMembers = new LinkedHashMap<>()
        tsv.eachLine { String line ->
            if (line.isEmpty()) return
            String[] parts = line.split('\t', 2)
            if (parts.length < 2) return
            String rep = parts[0]
            String member = parts[1]
            repToMembers.computeIfAbsent(rep, { new ArrayList<String>() }) << member
        }
        List<ProteinCluster> clusters = []
        int idx = 0
        repToMembers.each { rep, members ->
            ProteinRef repRef = new ProteinRef(genomeId: g.id, proteinId: rep)
            List<ProteinRef> memberRefs = members.collect { new ProteinRef(genomeId: g.id, proteinId: it) }
            clusters << new ProteinCluster(
                clusterId: "c${idx++}".toString(),
                representative: repRef,
                members: memberRefs,
            )
        }
        new ProteinClusterSet(clusters: clusters, identityThreshold: identity, coverageThreshold: coverage)
    }

    private static ProteinClusterSet singletonClusters(Genome g, List<Protein> proteins, double identity, double coverage) {
        List<ProteinCluster> clusters = []
        proteins.eachWithIndex { Protein p, int i ->
            ProteinRef ref = new ProteinRef(genomeId: g.id, proteinId: p.id)
            clusters << new ProteinCluster(clusterId: "c${i}".toString(), representative: ref, members: [ref])
        }
        new ProteinClusterSet(clusters: clusters, identityThreshold: identity, coverageThreshold: coverage)
    }

    private static void writeFasta(List<Protein> proteins, File out) {
        out.withWriter('UTF-8') { w ->
            for (Protein p : proteins) {
                w.writeLine('>' + p.id)
                w.writeLine(p.sequence ?: '')
            }
        }
    }
}
