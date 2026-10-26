package gspa.integration.contextbridge

import com.fasterxml.jackson.databind.DeserializationFeature
import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.learned.ExternalContext

/**
 * Reading and combining typed {@link ExternalContext} inputs. The explicit file is a strict JSON
 * rendering of the class (unknown keys are errors, values are validated); it is how
 * neighbourhood, FBA, predictor-status, role and assembly-completeness context that no producer in
 * the repository computes yet reaches the learned stage, without a code change.
 */
class ExternalContextSupport {

    private static final ObjectMapper MAPPER = new ObjectMapper()
        .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, true)

    static ExternalContext read(File file) {
        if (file == null || !file.isFile()) throw new IllegalArgumentException("external context file does not exist: ${file}")
        ExternalContext ctx
        try {
            ctx = MAPPER.readValue(file, ExternalContext)
        } catch (com.fasterxml.jackson.core.JsonProcessingException e) {
            throw new IllegalArgumentException("external context ${file} is not valid: ${e.originalMessage}", e)
        }
        try { ctx.validate() } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("external context ${file}: ${e.message}", e)
        }
        ctx
    }

    /**
     * Union of two contexts. A function key supplied by both for the same kind of context is an
     * error: neither source silently overrides the other. Either argument may be null.
     */
    static ExternalContext merge(ExternalContext explicit, ExternalContext produced, String producedLabel = 'cohort inference') {
        if (explicit == null) return produced
        if (produced == null) return explicit
        ExternalContext out = new ExternalContext()
        out.predictorStatus = explicit.predictorStatus
        out.assemblyCompleteness = explicit.assemblyCompleteness
        out.neighborhood.putAll(explicit.neighborhood)
        out.fba.putAll(explicit.fba)
        out.roles.putAll(explicit.roles)
        out.candidates.addAll(explicit.candidates)
        out.cohort.putAll(explicit.cohort)
        produced.cohort.each { String k, ExternalContext.Cohort c ->
            if (out.cohort.containsKey(k)) {
                throw new IllegalArgumentException("cohort features for ${k} are supplied both by the explicit context file and by ${producedLabel}")
            }
            out.cohort.put(k, c)
        }
        produced.neighborhood.each { String k, ExternalContext.Neighborhood n ->
            if (out.neighborhood.containsKey(k)) throw new IllegalArgumentException("neighborhood features for ${k} are supplied twice")
            out.neighborhood.put(k, n)
        }
        produced.fba.each { String k, ExternalContext.FbaConditional f ->
            if (out.fba.containsKey(k)) throw new IllegalArgumentException("conditional FBA features for ${k} are supplied twice")
            out.fba.put(k, f)
        }
        produced.roles.each { String k, List r -> if (out.roles.put(k, r) != null) throw new IllegalArgumentException("roles for ${k} are supplied twice") }
        out.candidates.addAll(produced.candidates)
        if (produced.predictorStatus != null) {
            if (out.predictorStatus != null) throw new IllegalArgumentException('predictor status is supplied twice')
            out.predictorStatus = produced.predictorStatus
        }
        if (produced.assemblyCompleteness != null) {
            if (out.assemblyCompleteness != null) throw new IllegalArgumentException('assembly completeness is supplied twice')
            out.assemblyCompleteness = produced.assemblyCompleteness
        }
        out.validate()
        out
    }
}
