package gspa.integration.learned

import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceType
import gspa.integration.IntegrationState
import gspa.model.AnnotationType
import gspa.model.Genome
import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Unroll

import static gspa.integration.learned.LearnedTestSupport.claim

/** Adapter roles must survive raw-claim export into the shared learned candidate path. */
class ClaimRoleBridgeSpec extends Specification {
    @TempDir File scratch
    def extractor = new ContextFeatureExtractor()
    def state = new IntegrationState(new Genome(id: 'g'))

    private static def observation(Object role, double score = 0.01d) {
        claim('p', 'reaction-subunit:path:rxn', score, 'gapsmith',
            EvidenceType.SEQUENCE_SIMILARITY, AnnotationType.CUSTOM, [candidate_role: role])
    }

    @Unroll
    def "supplied #role role survives candidate construction even at a weak score"() {
        when:
        def candidates = extractor.extractAll([observation(role)], state)

        then:
        candidates*.role == [CandidateRole.valueOf(role)]
        candidates[0].get('role.' + role.toLowerCase(Locale.ROOT)) == 1d
        candidates[0].get('seq.group.homology.max_raw') == 0.01d
        candidates[0].provenance.any { it.contains('candidate_role') && it.contains('gapsmith') && it.contains(role) }

        where:
        role << CandidateRole.values()*.name()
    }

    def "competing supplied roles are deterministic alternatives and duplicates add no candidates"() {
        given:
        def claims = [observation('REGULATOR'), observation('COMPLEX_SUBUNIT'), observation('REGULATOR')]

        expect:
        extractor.extractAll(claims, state)*.role == [CandidateRole.REGULATOR, CandidateRole.COMPLEX_SUBUNIT]
        extractor.extractAll(claims.reverse(), state)*.candidateKey() == extractor.extractAll(claims, state)*.candidateKey()
    }

    def "explicit structured role remains authoritative and raw role remains auditable"() {
        given:
        def c = observation('COMPLEX_SUBUNIT')
        def ctx = new ExternalContext(roles: [(c.functionKey()): [CandidateRole.ACCESSORY]])

        when:
        def candidates = extractor.extractAll([c], state, ctx)

        then:
        candidates*.role == [CandidateRole.ACCESSORY]
        candidates[0].provenance.any { it.contains('COMPLEX_SUBUNIT') && it.contains('structured context') }
    }

    @Unroll
    def "malformed supplied role #bad is rejected instead of silently changing candidate identity"() {
        when:
        extractor.extractAll([observation(bad)], state)

        then:
        def error = thrown(IllegalArgumentException)
        error.message.contains('candidate_role')
        error.message.contains('p|CUSTOM|reaction-subunit:path:rxn')

        where:
        bad << ['CATALYTIC_TYPO', '', ['COMPLEX_SUBUNIT'], 42]
    }

    def "missing role keeps conservative defaults without reading a name prefix"() {
        given:
        def c = observation(null)

        expect:
        extractor.extractAll([c], state)*.role == [CandidateRole.UNSPECIFIED]
        extractor.extractAll([claim('p', 'EC:1.1.1.1', 0.01d)], state)*.role == [CandidateRole.CATALYTIC]
    }

    def "raw export and direct integration preserve the supplied role and exact feature vector"() {
        given:
        def claims = [observation('COMPLEX_SUBUNIT')]
        def file = new File(scratch, 'claims.jsonl')
        ClaimExtractor.writeClaimsJsonl(claims, file)

        when:
        def reread = new ClaimExtractor().readClaimsJsonl(file)
        def direct = extractor.extractAll(reread, state)
        def original = extractor.extractAll(claims, state)

        then:
        direct*.role == [CandidateRole.COMPLEX_SUBUNIT]
        direct*.candidateKey() == original*.candidateKey()
        direct[0].values == original[0].values
        direct[0].provenance == original[0].provenance
    }
}
