package gspa.predictor.cluster

import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import gspa.model.ProteinRef
import spock.lang.Specification
import spock.lang.TempDir

class IntragenomeClustererSpec extends Specification {

    @TempDir
    File tmp

    def "parses mmseqs cluster TSV into a ProteinClusterSet"() {
        given:
        def g = new Genome(id: 'g1')
        // Simulate mmseqs linclust output: rep TAB member, one row per member.
        def tsv = new File(tmp, 'cluster.tsv')
        tsv.text = "p1\tp1\np1\tp2\np1\tp3\np4\tp4\n"

        when:
        def set = IntragenomeClusterer.parseClusterTsv(tsv, g, 0.9d, 0.8d)

        then:
        set.clusterCount() == 2
        set.proteinCount() == 4
        set.identityThreshold == 0.9d
        set.coverageThreshold == 0.8d
        def c1 = set.clusterFor(new ProteinRef(genomeId: 'g1', proteinId: 'p2'))
        c1 != null
        c1.representative == new ProteinRef(genomeId: 'g1', proteinId: 'p1')
        c1.members.size() == 3
        set.isRepresentative(new ProteinRef(genomeId: 'g1', proteinId: 'p1'))
        !set.isRepresentative(new ProteinRef(genomeId: 'g1', proteinId: 'p2'))
    }

    def "empty genome yields empty cluster set without invoking mmseqs"() {
        given:
        def clusterer = new IntragenomeClusterer()
        clusterer.commandRunner = { List<String> cmd -> throw new IllegalStateException("should not run") }
        def g = new Genome(id: 'g1')

        when:
        def set = clusterer.cluster([g], 0.9d, 0.8d)

        then:
        set.clusterCount() == 0
    }

    def "more than one genome is rejected — use CrossGenomeClusterer"() {
        given:
        def clusterer = new IntragenomeClusterer()

        when:
        clusterer.cluster([new Genome(id: 'a'), new Genome(id: 'b')], 0.9d, 0.8d)

        then:
        thrown(IllegalArgumentException)
    }

    def "falls back to singletons when mmseqs exits non-zero"() {
        given:
        def contig = new Contig(id: 'c1', proteins: [new Protein(id: 'p1', sequence: 'MAAA'), new Protein(id: 'p2', sequence: 'MBBB')])
        def g = new Genome(id: 'g1', contigs: [contig])
        def clusterer = new IntragenomeClusterer()
        clusterer.workDir = new File(tmp, 'work')
        clusterer.commandRunner = { List<String> cmd -> 1 }   // fail

        when:
        def set = clusterer.cluster([g], 0.9d, 0.8d)

        then:
        set.clusterCount() == 2
        set.clusters.every { it.members.size() == 1 }
    }
}
