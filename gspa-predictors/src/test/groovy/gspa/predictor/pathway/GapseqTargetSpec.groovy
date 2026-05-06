package gspa.predictor.pathway

import gspa.model.Contig
import gspa.model.Genome
import gspa.model.Protein
import spock.lang.Specification
import spock.lang.TempDir

class GapseqTargetSpec extends Specification {

    @TempDir
    File tmp

    private Genome tinyGenome() {
        new Genome(id: 'g1', contigs: [
            new Contig(id: 'c1', sequence: 'ATGAAA', proteins: [
                new Protein(id: 'p1', sequence: 'MKS'),
                new Protein(id: 'p2', sequence: 'MQA'),
                new Protein(id: 'p3', sequence: 'MQA'),
            ])
        ])
    }

    def "target=genome writes contig sequence, not proteins"() {
        given:
        def gp = new GapseqPredictor(target: 'genome')

        when: "reflectively call the private writer"
        def m = gp.class.getDeclaredMethod('writeTargetFasta', Genome, File); m.accessible = true
        def out = (File) m.invoke(gp, tinyGenome(), tmp)

        then:
        out.name == 'genome.fna'
        out.text.contains('>c1')
        !out.text.contains('>p1')
    }

    def "target=proteome writes every protein"() {
        given:
        def gp = new GapseqPredictor(target: 'proteome')

        when:
        def m = gp.class.getDeclaredMethod('writeTargetFasta', Genome, File); m.accessible = true
        def out = (File) m.invoke(gp, tinyGenome(), tmp)

        then:
        out.name == 'proteins.faa'
        out.text.contains('>p1')
        out.text.contains('>p2')
        out.text.contains('>p3')
    }

    def "target=reps writes only representative proteins"() {
        given:
        def gp = new GapseqPredictor(target: 'reps', representativeIds: ['p1', 'p3'] as Set)

        when:
        def m = gp.class.getDeclaredMethod('writeTargetFasta', Genome, File); m.accessible = true
        def out = (File) m.invoke(gp, tinyGenome(), tmp)

        then:
        out.name == 'reps.faa'
        out.text.contains('>p1')
        out.text.contains('>p3')
        !out.text.contains('>p2')
    }

    def "target=reps without representativeIds throws at predictGenome"() {
        given:
        def gp = new GapseqPredictor(target: 'reps', representativeIds: null)

        when:
        gp.predictGenome(tinyGenome())

        then:
        thrown(IllegalStateException)
    }

    def "unknown target is rejected"() {
        given:
        def gp = new GapseqPredictor(target: 'bogus')

        when:
        def m = gp.class.getDeclaredMethod('writeTargetFasta', Genome, File); m.accessible = true
        m.invoke(gp, tinyGenome(), tmp)

        then:
        def e = thrown(Exception)
        e.cause instanceof IllegalArgumentException || e instanceof IllegalArgumentException
    }
}
