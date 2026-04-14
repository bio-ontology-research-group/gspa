package gspa.predictor.cluster

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import gspa.model.ProteinCluster
import gspa.model.ProteinClusterSet
import gspa.model.ProteinRef
import spock.lang.Specification

class ClusterAnnotationPropagatorSpec extends Specification {

    def "copies rep's annotations to members with propagation metadata"() {
        given:
        def rep = new Protein(id: 'p1')
        def m2 = new Protein(id: 'p2')
        def m3 = new Protein(id: 'p3')
        rep.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:0001', source: 'diamond'))
        rep.annotations.add(new Annotation(type: AnnotationType.PFAM, value: 'PF001', source: 'hmmer'))

        def g = new Genome(id: 'g1', contigs: [new Contig(id: 'c1', proteins: [rep, m2, m3])])
        def clusters = new ProteinClusterSet(clusters: [
            new ProteinCluster(
                clusterId: 'c0',
                representative: new ProteinRef(genomeId: 'g1', proteinId: 'p1'),
                members: [
                    new ProteinRef(genomeId: 'g1', proteinId: 'p1'),
                    new ProteinRef(genomeId: 'g1', proteinId: 'p2'),
                    new ProteinRef(genomeId: 'g1', proteinId: 'p3'),
                ]
            )
        ], identityThreshold: 0.9, coverageThreshold: 0.8)

        when:
        int copied = new ClusterAnnotationPropagator().propagate(clusters, [g])

        then:
        copied == 4     // 2 annotations × 2 non-rep members
        m2.annotations.size() == 2
        m3.annotations.size() == 2
        rep.annotations.size() == 2  // rep unchanged
        // Propagated annotations carry provenance.
        m2.annotations.annotations.every { it.metadata.propagatedFrom == 'p1' }
        m2.annotations.annotations.every { it.metadata.clusterId == 'c0' }
        // Rep's own annotations are NOT stamped.
        rep.annotations.annotations.every { (it.metadata ?: [:]).propagatedFrom == null }
    }

    def "filters by sourcePredictor when specified"() {
        given:
        def rep = new Protein(id: 'p1')
        def m2 = new Protein(id: 'p2')
        rep.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:1', source: 'diamond'))
        rep.annotations.add(new Annotation(type: AnnotationType.PFAM, value: 'PF1', source: 'hmmer'))

        def g = new Genome(id: 'g1', contigs: [new Contig(id: 'c1', proteins: [rep, m2])])
        def clusters = new ProteinClusterSet(clusters: [
            new ProteinCluster(
                clusterId: 'c0',
                representative: new ProteinRef(genomeId: 'g1', proteinId: 'p1'),
                members: [
                    new ProteinRef(genomeId: 'g1', proteinId: 'p1'),
                    new ProteinRef(genomeId: 'g1', proteinId: 'p2'),
                ]
            )
        ])

        when: "propagate only diamond annotations"
        int copied = new ClusterAnnotationPropagator().propagate(clusters, [g], 'diamond')

        then:
        copied == 1
        m2.annotations.size() == 1
        m2.annotations.annotations[0].value == 'GO:1'
    }

    def "singleton clusters produce no copies"() {
        given:
        def rep = new Protein(id: 'p1')
        rep.annotations.add(new Annotation(type: AnnotationType.GO, value: 'GO:1', source: 'diamond'))
        def g = new Genome(id: 'g1', contigs: [new Contig(id: 'c1', proteins: [rep])])
        def clusters = new ProteinClusterSet(clusters: [
            new ProteinCluster(
                clusterId: 'c0',
                representative: new ProteinRef(genomeId: 'g1', proteinId: 'p1'),
                members: [new ProteinRef(genomeId: 'g1', proteinId: 'p1')],
            )
        ])

        expect:
        new ClusterAnnotationPropagator().propagate(clusters, [g]) == 0
        rep.annotations.size() == 1   // rep unchanged
    }
}
