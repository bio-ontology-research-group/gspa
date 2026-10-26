package gspa.integration.cohort

import gspa.ontology.ReactionGraph
import groovy.transform.CompileStatic

/**
 * Direct EC -> reaction candidate generation. Does not use EC -> GO: an EC
 * call maps to every reaction carrying that EC in the {@link ReactionGraph}.
 * The reaction count is kept so downstream features can see how
 * (un)specific the EC support is; the confidence is not divided.
 */
@CompileStatic
class EcCandidateExpander {

    static List<CohortInput.CandidateClaim> expand(Collection<CohortInput.EcCall> calls, ReactionGraph graph) {
        List<CohortInput.CandidateClaim> out = new ArrayList<>()
        for (CohortInput.EcCall call : calls) {
            Set<String> rxns = graph.reactionsForEc(call.ecNumber)
            // Sorted for deterministic output.
            for (String r : new TreeSet<String>(rxns)) {
                out.add(new CohortInput.CandidateClaim(
                        genomeId: call.genomeId, proteinId: call.proteinId, reactionId: r,
                        baselineConfidence: call.confidence, ecNumber: call.ecNumber, ecDirect: true,
                        ecReactionCount: rxns.size(),
                        provenance: call.provenance ?: 'ec-direct'))
            }
        }
        out
    }
}
