package gspa.model

import spock.lang.Specification

class ProteinClusterSetSpec extends Specification {

    def "resolves protein → cluster → representative and flags reps"() {
        given:
        def rep = new ProteinRef('g1', 'p1')
        def m2 = new ProteinRef('g1', 'p2')
        def m3 = new ProteinRef('g1', 'p3')
        def solo = new ProteinRef('g1', 'p4')
        def c1 = new ProteinCluster(clusterId: 'c1', representative: rep, members: [rep, m2, m3])
        def c2 = new ProteinCluster(clusterId: 'c2', representative: solo, members: [solo])
        def set = new ProteinClusterSet(clusters: [c1, c2], identityThreshold: 0.9d, coverageThreshold: 0.8d)

        expect:
        set.clusterFor(m2) == c1
        set.clusterFor(rep) == c1
        set.representativeOf(m3) == rep
        set.isRepresentative(rep)
        !set.isRepresentative(m2)
        set.isRepresentative(solo)
        c2.singleton
        set.representatives() == [rep, solo]
        set.clusterCount() == 2
        set.proteinCount() == 4
    }

    def "representativeOf returns the protein itself when unclustered"() {
        given:
        def set = new ProteinClusterSet(clusters: [], identityThreshold: 0.9d, coverageThreshold: 0.8d)
        def p = new ProteinRef('g1', 'orphan')

        expect:
        set.clusterFor(p) == null
        set.representativeOf(p) == p
    }
}
