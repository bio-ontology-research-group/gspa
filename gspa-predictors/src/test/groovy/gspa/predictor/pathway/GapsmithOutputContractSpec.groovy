package gspa.predictor.pathway

import groovy.json.JsonSlurper
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceType
import gspa.model.AnnotationType
import gspa.predictor.ToolOutputContractException
import spock.lang.Specification
import spock.lang.TempDir

class GapsmithOutputContractSpec extends Specification {
    @TempDir File temp
    String table() { getClass().getResource('/observed-output/gapsmith/training-subset-Reactions.tbl').text }
    File write(String body) { def f = new File(temp,'real-Reactions.tbl'); f.text=body; f }

    def 'observed target identifiers and EC values survive the real parser and claim extraction'() {
        given:
        write(table())
        def p = new GapseqPredictor()
        when:
        def result = p.parseOutput(temp)
        def annotation = result['ABB58505.1']?.find { it.type==AnnotationType.EC }
        then:
        annotation != null
        annotation.value == 'EC:4.3.2.1'
        Math.abs(annotation.score - 0.473d) < 1e-12
        annotation.metadata.pident == 47.3d
        annotation.metadata.bitscore == 388d
        annotation.metadata.referenceId == 'sp|P11447|ARLY_ECOLI'
        !result.containsKey('sp|P11447|ARLY_ECOLI')
        new ClaimExtractor().fromAnnotation('ABB58505.1',annotation).evidenceType == EvidenceType.SEQUENCE_SIMILARITY
    }

    def 'weak observed hits remain candidates and complexes are not standalone EC assignments'() {
        given:
        write(table())
        when:
        def result = new GapseqPredictor().parseOutput(temp)
        then:
        result['ABB58544.1']?.any { it.value=='EC:2.1.3.11' && Math.abs(it.score-0.277d)<1e-12 }
        result['ABB56978.1']?.any { it.type==AnnotationType.CUSTOM && it.value=='reaction-subunit:2FLHYD:2FLABCTRANSPORT' }
        result['ABB58152.1']?.every { it.type!=AnnotationType.EC }
        result['ABB58152.1']?.any { it.value=='reaction-subunit:ARGSYN-SUCC-PWY:CARBPSYN-RXN' && it.metadata.candidate_role=='COMPLEX_SUBUNIT' }
        result['ABB57067.1']?.any { it.type==AnnotationType.CUSTOM && it.value=='reaction:AGMA-PROD:rxn1' }
        result.size() == 5
    }

    def 'all rows including missing search evidence remain available for a separate context report'() {
        given:
        write(table())
        def p = new GapseqPredictor()
        when:
        p.parseOutput(temp)
        then:
        p.reactionObservations.size() == 8
        p.reactionObservations*.status.toSet() == ['good_blast','bad_blast','no_blast','no_seq_data','spontaneous'] as Set
        p.reactionTableSha256 == java.security.MessageDigest.getInstance('SHA-256').digest(table().bytes).encodeHex().toString()
    }

    def 'duplicate pathway evidence collapses one candidate while retaining every original row'() {
        given:
        def lines=table().readLines()
        write(lines.join('\n')+'\n'+lines[4]+'\n')
        when:
        def result = new GapseqPredictor().parseOutput(temp)
        then:
        result['ABB58505.1'].size() == 1
        new JsonSlurper().parseText(result['ABB58505.1'][0].metadata.hitsJson).size() == 2
    }

    def 'column order is named and missing required columns are rejected'() {
        given:
        def lines=table().readLines().collect { it.split('\t',-1).toList().reverse().join('\t') }
        write(lines.join('\n')+'\n')
        when:
        def result = new GapseqPredictor().parseOutput(temp)
        then:
        result['ABB58505.1']*.value == ['EC:4.3.2.1']
        when:
        write(table().replaceFirst('pathway\t','wrong\t'))
        new GapseqPredictor().parseOutput(temp)
        then:
        thrown(ToolOutputContractException)
    }

    def 'malformed observed statistics and statuses are rejected: #old to #replacement'() {
        given:
        write(table().replace(old,replacement))
        when:
        new GapseqPredictor().parseOutput(temp)
        then:
        thrown(ToolOutputContractException)
        where:
        old                  | replacement
        '\t47.3\t'           | '\tNaN\t'
        '\t47.3\t'           | '\t101\t'
        '\t388\t'            | '\t-1\t'
        '\t9.850e-133\t'      | '\t-1\t'
        '\t99.1\t'           | '\tInfinity\t'
        '\tgood_blast\t'     | '\tunknown_status\t'
        '\tFALSE\t'          | '\tmaybe\t'
        '\t4.3.2.1\t'        | '\tgarbage\t'
    }

    def 'pathway-local reaction IDs never merge unrelated custom reactions'() {
        given:
        def lines=table().readLines()
        write(table()+lines[3].replace('AGMA-PROD', 'OTHER-PATHWAY')+'\n')
        when:
        def result = new GapseqPredictor().parseOutput(temp)
        then:
        result['ABB57067.1']*.value.toSet() == ['reaction:AGMA-PROD:rxn1','reaction:OTHER-PATHWAY:rxn1'] as Set
    }

    def 'reference-local IDs are retained as reaction candidates and partial ECs stay partial: #ec'() {
        given:
        write(table().replace('4.3.2.1',ec))
        when:
        def result = new GapseqPredictor().parseOutput(temp)
        def annotation = result['ABB58505.1'][0]
        then:
        annotation.type.name()==kind
        annotation.value==value
        new JsonSlurper().parseText(annotation.metadata.hitsJson)[0].ec==ec
        where:
        ec          | kind     | value
        '1.5.1.M4'  | 'CUSTOM' | 'reaction:ARGSYN-SUCC-PWY:ARGSUCCINLYA-RXN'
        '2.5.1.bl'  | 'CUSTOM' | 'reaction:ARGSYN-SUCC-PWY:ARGSUCCINLYA-RXN'
        '2.4.1'     | 'EC'     | 'EC:2.4.1'
        '1.1'       | 'EC'     | 'EC:1.1'
        '1.1.1.n2'  | 'EC'     | 'EC:1.1.1.n2'
    }

    def 'truncated rows and multiple result tables are rejected'() {
        given:
        write(table().readLines()[0]+'\nP\tR\n')
        when:
        new GapseqPredictor().parseOutput(temp)
        then:
        thrown(ToolOutputContractException)
        when:
        write(table())
        new File(temp,'other-Reactions.tbl').text=table()
        new GapseqPredictor().parseOutput(temp)
        then:
        thrown(ToolOutputContractException)
    }
}
