package gspa.integration.workflow

import gspa.integration.contextbridge.PathwayGapBuilder
import spock.lang.Specification
import spock.lang.TempDir

class PathwayGapBuilderSpec extends Specification {
    @TempDir File dir
    static Map row(Map extra = [:]) {
        [pathway:'|P|', rxn:'R', ec:'1.1.1.1', status:'no_blast', spont:'FALSE', is_complex:'FALSE'] + extra
    }
    static Map observation(List rows = [row()], Map extra = [:]) {
        [schema:PathwayGapBuilder.OBSERVATION_SCHEMA, genome_id:'g', source:'gapsmith', target:'proteome',
         outcome:'SUCCESS', reaction_table_sha256:'a'*64, pathway_table_sha256:'b'*64,
         pathway_rows:[[ID:'P', Prediction:'TRUE']], reaction_rows:rows] + extra
    }

    def 'all evidence for a reaction is considered before gap eligibility'() {
        when:
        def result = PathwayGapBuilder.build(observation([row(), row(veto)]), 'g')

        then:
        result.available
        result.gaps.isEmpty()
        result.report.excluded_reactions[reason] == 1

        where:
        veto                        | reason
        [status:'good_blast']        | 'strong_hit_present'
        [status:'no_seq_data']       | 'reference_unknown'
        [status:'spontaneous']       | 'spontaneous'
        [spont:'TRUE']              | 'spontaneous'
        [is_complex:'TRUE']         | 'complex_requires_subunit_inference'
    }

    def 'weak-hit flag and EC GO alternatives are order independent with no duplicates'() {
        given:
        def rows = [row(), row(ec:'EC:1.1.1.1/2.2.2.2', status:'bad_blast'), row()]
        def mapping = ['1.1.1.1':['GO:0000001', 'GO:0000002', 'GO:0000001']]

        when:
        def forward = PathwayGapBuilder.build(observation(rows), 'g', mapping)
        def reverse = PathwayGapBuilder.build(observation(rows.reverse()), 'g', mapping)

        then:
        forward.gaps == reverse.gaps
        forward.gaps.size() == 3
        forward.gaps.every { it.gapseqGuessed && it.pathwayId == 'gapsmith:P' }
        forward.gaps.find { it.ecNumber == '2.2.2.2' }.goTerm == null
        forward.report.eligible_reactions == 1
    }

    def 'unpredicted pathway is not inferred from reaction rows'() {
        expect:
        PathwayGapBuilder.build(observation([row()], [pathway_rows:paths]), 'g').gaps.isEmpty()

        where:
        paths << [[], [[ID:'P',Prediction:'FALSE']], [[ID:'Q',Prediction:'TRUE']]]
    }

    def 'incomplete target or missing pathway observations stay unavailable'() {
        when:
        def result = PathwayGapBuilder.build(observation([row()], overrides), 'g')

        then:
        !result.available
        result.gaps.isEmpty()
        result.report.status == 'unavailable'

        where:
        overrides << [[outcome:'PARTIAL'], [outcome:'FAILED'], [outcome:'UNAVAILABLE'],
                      [target:'reps'], [target:'genome'], [pathway_rows:null], [pathway_table_sha256:null]]
    }

    def 'empty completed tables are observed zero, not missing'() {
        when:
        def result = PathwayGapBuilder.build(observation([], [pathway_rows:[], outcome:'SUCCESS_ZERO']), 'g')

        then:
        result.available
        result.gaps.isEmpty()
        result.report.candidate_gap_records == 0
    }

    def 'bad observations fail explicitly'() {
        when:
        PathwayGapBuilder.build(observation(rows, overrides), 'g')

        then:
        thrown(IllegalArgumentException)

        where:
        rows                         | overrides
        [row()]                      | [genome_id:'another-genome']
        [row()]                      | [reaction_table_sha256:'not-a-hash']
        [row(status:'unknown')]       | [:]
        [row(spont:'maybe')]          | [:]
        [row(ec:'garbage')]           | [:]
        [row()]                      | [pathway_rows:[[ID:'P',Prediction:'TRUE'], [ID:'|P|',Prediction:'FALSE']]]
    }

    def 'partial ECs survive and reference-local tokens remain reaction-only hypotheses'() {
        when:
        def result = PathwayGapBuilder.build(observation([row(ec:ec)]), 'g')

        then:
        result.gaps*.ecNumber == expected

        where:
        ec             | expected
        '1.1'          | ['1.1']
        '1.1.1.n5'     | ['1.1.1.n5']
        '1.5.1.M4'     | [null]
        ''             | [null]
    }

    def 'EC mapping reader retains multiple authoritative GO mappings'() {
        given:
        def file = new File(dir, 'ec2go')
        file.text = '! comment\nEC:1.1.1.1 > GO:one ; GO:0000001\nEC:1.1.1.1 > GO:two ; GO:0000002\nEC:1.1.1.1\tGO:0000001\n'

        expect:
        PathwayGapBuilder.ec2go(file)['1.1.1.1'] == ['GO:0000001', 'GO:0000002']
    }
}
