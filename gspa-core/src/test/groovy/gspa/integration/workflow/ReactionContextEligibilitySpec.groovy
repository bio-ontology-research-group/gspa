package gspa.integration.workflow

import gspa.integration.ClaimExtractor
import spock.lang.Specification
import spock.lang.TempDir

class ReactionContextEligibilitySpec extends Specification {
    @TempDir File dir

    def 'real workflow records RLC as skipped when no gap inputs exist'() {
        given:
        def fixture = new ContextFixtures(dir)
        def files = fixture.genome('g')
        def config = new IntegrationWorkflowConfig(lite: true, enabledPriors: 'none',
            genomeId: 'g', reactionGraphFile: fixture.reactions, reactionEcAliasesFile: fixture.aliases,
            genomeLayoutFile: files.layout, rlcSuggester: true)

        when:
        def result = new IntegrationWorkflow(config).run(new ClaimExtractor().readClaimsJsonl(files.claims))

        then:
        result.audit.stages.reaction_context.status == 'skipped'
        result.audit.stages.reaction_context.requested
        result.audit.stages.reaction_context.reason == 'no metabolic gaps to evaluate'
        result.audit.stages.evidence_integration.status == 'executed'
        result.integrated.suggestions.empty
        !result.integrated.provenance.empty
    }
}
