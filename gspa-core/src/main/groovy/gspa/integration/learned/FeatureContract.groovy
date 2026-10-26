package gspa.integration.learned

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

import java.security.MessageDigest

/**
 * Typed, versioned feature contract for the learned context model: the ordered list of
 * feature columns, their provenance group and documentation. A trained artifact is bound
 * to a contract by id and fingerprint; extractor, trainer and inference must agree on it.
 *
 * <h3>Conventions (apply to every feature)</h3>
 * <ul>
 *   <li><b>Missing is not absent.</b> Whenever information may be unavailable the value
 *       column is paired with an availability/reported indicator; an unavailable value is
 *       stored as 0 and the indicator is 0. The model can therefore tell "no neighbours
 *       measured" from "neighbours measured and none support".</li>
 *   <li><b>No double counting.</b> Predictor evidence enters only as the strongest raw score
 *       per source and per evidence (correlation) group, so repeated or duplicated hits from
 *       one predictor add nothing.</li>
 *   <li><b>No posteriors as evidence.</b> Nothing here is read from the integrator's
 *       posterior log-odds, likelihoods, pinned floors, orthogroup consensus posteriors or
 *       previously calibrated probabilities. Metabolic-gap information is a binary listing
 *       indicator, never a Bayes factor.</li>
 *   <li><b>Conditional means conditional.</b> FBA features describe behaviour under the
 *       tested objective/media only, and are suppressed when the objective itself requires
 *       the candidate.</li>
 * </ul>
 */
class FeatureContract {

    static final String CANDIDATE_ID = 'gspa.context-features/2'
    static final String LOCUS_GROUP_ID = 'gspa.locus-group-features/1'

    /** Evidence (correlation) groups that carry functional sequence evidence. */
    static final List<String> SEQUENCE_GROUPS = ['homology', 'ml_protein_seq', 'structure', 'motif', 'domain_specific']

    /** Sources with their own columns; every other source maps to {@code other}. */
    static final List<String> SOURCES = [
        'diamond', 'mmseqs2', 'interproscan', 'hmmer', 'eggnog-mapper', 'foldseek',
        'deepgo-plusplus', 'deepgo-plusplus-light', 'proteinfer', 'clean', 'deepfri', 'esm2-go', 'other',
    ]

    final String id
    final List<FeatureSpec> features

    FeatureContract(String id, List<FeatureSpec> features) {
        this.id = id
        this.features = Collections.unmodifiableList(new ArrayList<>(features))
        if (this.features*.name.toSet().size() != this.features.size()) {
            throw new IllegalArgumentException("duplicate feature names in contract ${id}")
        }
    }

    List<String> names() { features*.name }

    int size() { features.size() }

    int indexOf(String name) { features.findIndexOf { it.name == name } }

    /** SHA-256 over id + ordered (name, group): descriptions may be edited without invalidating artifacts. */
    String fingerprint() {
        def md = MessageDigest.getInstance('SHA-256')
        md.update(id.getBytes('UTF-8'))
        features.each { md.update(("\n" + it.name + '\t' + it.group.name()).getBytes('UTF-8')) }
        md.digest().collect { String.format('%02x', it & 0xff) }.join()
    }

    /** Canonical JSON (stable key order); the Python trainer reads this file. */
    String toJson() {
        JsonOutput.prettyPrint(JsonOutput.toJson([
            id         : id,
            fingerprint: fingerprint(),
            features   : features.collect { [name: it.name, group: it.group.name(), description: it.description] },
        ])) + '\n'
    }

    static FeatureContract fromJson(String json) {
        def m = new JsonSlurper().parseText(json) as Map
        def specs = (m.features as List<Map>).collect {
            new FeatureSpec(name: it.name as String, group: FeatureGroup.valueOf(it.group as String),
                description: it.description as String)
        }
        def c = new FeatureContract(m.id as String, specs)
        if (m.fingerprint != null && m.fingerprint != c.fingerprint()) {
            throw new IllegalArgumentException("contract ${c.id}: stated fingerprint does not match its feature list")
        }
        c
    }

    // ------------------------------------------------------------------ the production contract

    static final FeatureContract CURRENT = buildCandidateContract()

    /** Group-level features of a reaction-locus group, input of the null (no-candidate) outcome. */
    static final FeatureContract LOCUS_GROUP = buildLocusGroupContract()

    private static FeatureContract buildCandidateContract() {
        List<FeatureSpec> f = []
        def add = { String n, FeatureGroup g, String d -> f << new FeatureSpec(name: n, group: g, description: d) }

        add('fn.type_ec', FeatureGroup.FUNCTION_IDENTITY, '1 if the function is an EC number')
        add('fn.type_go', FeatureGroup.FUNCTION_IDENTITY, '1 if the function is a GO term')
        add('fn.type_other', FeatureGroup.FUNCTION_IDENTITY, '1 for any other function namespace (KEGG, COG, AMR, ...)')
        add('fn.go_aspect_mf', FeatureGroup.FUNCTION_IDENTITY, '1 if a GO claim is explicitly labelled MF')
        add('fn.go_aspect_bp', FeatureGroup.FUNCTION_IDENTITY, '1 if a GO claim is explicitly labelled BP')
        add('fn.go_aspect_cc', FeatureGroup.FUNCTION_IDENTITY, '1 if a GO claim is explicitly labelled CC')
        add('fn.go_aspect_unknown', FeatureGroup.FUNCTION_IDENTITY, '1 for GO claims with no aspect label (never assumed BP)')
        add('fn.in_complex_term', FeatureGroup.FUNCTION_IDENTITY, '1 if the GO term is a protein-complex term in state.complexClassification')
        add('fn.complex_heteromeric', FeatureGroup.FUNCTION_IDENTITY, '1 if that complex term is classified heteromeric (n)')
        add('fn.complex_ambiguous', FeatureGroup.FUNCTION_IDENTITY, '1 if that complex term is classified ambiguous (a): neither homomeric nor heteromeric is asserted')

        CandidateRole.values().each {
            add(it.featureName(), FeatureGroup.ROLE, "1 if the candidate role is ${it.name()}; roles are distinct candidates")
        }

        SEQUENCE_GROUPS.each { g ->
            add("seq.group.${g}.max_raw", FeatureGroup.BASELINE_SEQUENCE,
                "strongest raw predictor score in evidence group '${g}' (0 if none reported)")
            add("seq.group.${g}.reported", FeatureGroup.BASELINE_SEQUENCE,
                "1 if any predictor in group '${g}' reported this function for this protein")
        }
        SOURCES.each { s ->
            add("seq.src.${s}.max_raw", FeatureGroup.BASELINE_SEQUENCE,
                "strongest raw score reported by source '${s}' for this function (0 if it reported none)")
            add("seq.src.${s}.reported", FeatureGroup.BASELINE_SEQUENCE,
                "1 if source '${s}' reported this function; 0 does NOT mean the source ran (see seq.cov.*)")
        }
        add('seq.n_sources_reporting', FeatureGroup.BASELINE_SEQUENCE, 'number of distinct sources reporting this function (not number of hits)')
        add('seq.n_groups_reporting', FeatureGroup.BASELINE_SEQUENCE, 'number of distinct evidence groups reporting this function')
        add('seq.max_raw_any', FeatureGroup.BASELINE_SEQUENCE, 'strongest raw score over all baseline sequence claims')
        add('seq.best_pident_available', FeatureGroup.BASELINE_SEQUENCE, '1 if any claim carries a percent-identity metadata value')
        add('seq.best_pident', FeatureGroup.BASELINE_SEQUENCE, 'best percent identity / 100 over claims carrying it (0 if unavailable)')

        add('seq.cov.status_available', FeatureGroup.RUN_COVERAGE, '1 if predictor completion statuses were supplied')
        add('seq.cov.n_sources_not_fully_run', FeatureGroup.RUN_COVERAGE,
            'number of predictors with outcome FAILED, UNAVAILABLE or PARTIAL (absence from them is uninformative); 0 if statuses not supplied')
        add('seq.cov.n_sources_ran_clean', FeatureGroup.RUN_COVERAGE,
            'number of predictors with outcome SUCCESS or SUCCESS_ZERO; 0 if statuses not supplied')

        add('ctx.operon.available', FeatureGroup.CONTEXT_NEIGHBORHOOD, '1 if operon calls were supplied (state.operons non-empty)')
        add('ctx.operon.member', FeatureGroup.CONTEXT_NEIGHBORHOOD, '1 if the protein is in a called operon')
        add('ctx.operon.log_size', FeatureGroup.CONTEXT_NEIGHBORHOOD, 'log(1 + number of genes in the protein\'s operon); 0 if not a member')
        add('ctx.operon.n_members_with_seq_support', FeatureGroup.CONTEXT_NEIGHBORHOOD,
            'other operon members with any baseline sequence claim of any function')
        add('ctx.operon.n_members_dark', FeatureGroup.CONTEXT_NEIGHBORHOOD, 'other operon members with no claim at all')
        add('ctx.layout.available', FeatureGroup.CONTEXT_NEIGHBORHOOD, '1 if gene coordinates (state.genomeLayout) are available')
        add('ctx.layout.n_neighbors_10kb', FeatureGroup.CONTEXT_NEIGHBORHOOD, 'other genes within +-10 kb (midpoint) on the same contig')
        add('ctx.layout.same_strand_fraction', FeatureGroup.CONTEXT_NEIGHBORHOOD, 'fraction of those neighbours on the candidate\'s strand (0 if none)')
        add('ctx.competitor.max_raw', FeatureGroup.CONTEXT_NEIGHBORHOOD,
            'strongest raw baseline score for the SAME function on another local protein (same operon or +-10 kb); paralog/decoy signal')
        add('ctx.competitor.n_with_claim', FeatureGroup.CONTEXT_NEIGHBORHOOD,
            'number of other local proteins with any claim for the same function')
        add('ctx.gap.available', FeatureGroup.CONTEXT_NEIGHBORHOOD, '1 if metabolic gaps were supplied (state.metabolicGaps non-empty)')
        add('ctx.gap.function_listed', FeatureGroup.CONTEXT_NEIGHBORHOOD,
            '1 if an open (not already closed) metabolic gap lists this EC/GO function; a listing indicator, NOT a Bayes factor')
        add('ctx.gap.gapseq_guessed_listed', FeatureGroup.CONTEXT_NEIGHBORHOOD,
            '1 if such a gap was one for which gapseq itself proposed a gene (a different, weaker kind of listing)')
        add('ctx.nbr.available', FeatureGroup.CONTEXT_NEIGHBORHOOD, '1 if structured neighbourhood features were supplied for this candidate')
        add('ctx.nbr.support_max', FeatureGroup.CONTEXT_NEIGHBORHOOD,
            'supplied: strongest raw baseline support among neighbours for functions in the candidate\'s reaction/pathway neighbourhood')
        add('ctx.nbr.n_supporting', FeatureGroup.CONTEXT_NEIGHBORHOOD, 'supplied: number of neighbours with such support')
        add('ctx.nbr.conflict_max', FeatureGroup.CONTEXT_NEIGHBORHOOD,
            'supplied: strongest raw support among neighbours for a function incompatible with the candidate (unfavourable context)')
        add('ctx.nbr.reaction_distance_available', FeatureGroup.CONTEXT_NEIGHBORHOOD, '1 if a reaction-graph distance was supplied')
        add('ctx.nbr.reaction_distance', FeatureGroup.CONTEXT_NEIGHBORHOOD, 'supplied: reaction-graph distance to the nearest supported neighbour reaction (0 if unavailable)')

        add('cohort.available', FeatureGroup.CONTEXT_COHORT, '1 if cohort recurrence was supplied for this candidate')
        add('cohort.log_n_genomes', FeatureGroup.CONTEXT_COHORT, 'log(1 + genomes in the cohort); a count, never normalised to certainty')
        add('cohort.frac_with_orthogroup', FeatureGroup.CONTEXT_COHORT, 'fraction of cohort genomes with the candidate\'s orthogroup (0 if unavailable)')
        add('cohort.frac_with_function_seq', FeatureGroup.CONTEXT_COHORT, 'fraction of those with baseline sequence support for the function (not posteriors)')
        add('cohort.frac_with_context', FeatureGroup.CONTEXT_COHORT, 'fraction of those where the function is also supported by neighbourhood context')
        add('cohort.n_genomes_incomplete', FeatureGroup.CONTEXT_COHORT, 'supplied count of cohort genomes with incomplete assembly (absence uninformative)')

        CohortFeatureColumns.declare(add, FeatureGroup.CONTEXT_COHORT)

        add('asm.available', FeatureGroup.ASSEMBLY, '1 if the genome, the protein\'s contig sequence and valid CDS coordinates are available; the four asm.* values below are 0 otherwise')
        add('asm.log_contig_length', FeatureGroup.ASSEMBLY, 'log(1 + length) of the protein\'s contig')
        add('asm.log_distance_to_contig_end', FeatureGroup.ASSEMBLY, 'log(1 + bp from the CDS to the nearest contig end); small = truncation risk')
        add('asm.n_contigs', FeatureGroup.ASSEMBLY, 'log(1 + number of contigs): fragmentation')
        add('asm.log_n50', FeatureGroup.ASSEMBLY, 'log(1 + N50) of the assembly')
        add('asm.completeness_available', FeatureGroup.ASSEMBLY, '1 if an external completeness estimate (e.g. CheckM2) was supplied')
        add('asm.completeness', FeatureGroup.ASSEMBLY, 'supplied completeness / 100 (0 if unavailable)')

        add('fba.available', FeatureGroup.FBA_CONDITIONAL, '1 if conditional FBA features were supplied and usable for this candidate')
        add('fba.suppressed_objective_forced', FeatureGroup.FBA_CONDITIONAL,
            '1 if the tested objective itself required the candidate, so its FBA values were NOT used as evidence')
        add('fba.n_conditions', FeatureGroup.FBA_CONDITIONAL, 'log(1 + number of media/objective conditions tested)')
        add('fba.required_fraction', FeatureGroup.FBA_CONDITIONAL, 'fraction of tested conditions in which the candidate\'s reaction carried flux (conditional on those conditions)')
        add('fba.optional_fraction_available', FeatureGroup.FBA_CONDITIONAL, '1 if the share of tested conditions with optional (non-required) flux was reported')
        add('fba.optional_fraction', FeatureGroup.FBA_CONDITIONAL, 'fraction of tested conditions in which the reaction could carry flux but was not required (conditional on those conditions; 0 if unavailable)')
        add('fba.knockout_growth_effect', FeatureGroup.FBA_CONDITIONAL, 'supplied relative growth change on knockout (conditional; sign preserved)')

        new FeatureContract(CANDIDATE_ID, f)
    }

    private static FeatureContract buildLocusGroupContract() {
        List<FeatureSpec> f = []
        def add = { String n, String d -> f << new FeatureSpec(name: n, group: FeatureGroup.CONTEXT_NEIGHBORHOOD, description: d) }
        add('grp.max_margin', 'largest per-candidate margin (logit of the frozen candidate model) in the supplied group')
        add('grp.second_margin_available', '1 if the group has at least two candidates')
        add('grp.second_margin', 'second-largest margin (0 if unavailable)')
        add('grp.log1p_n_candidates', 'log(1 + number of supplied candidates); a learned weight, never a fixed cap')
        add('grp.any_seq_support', '1 if any candidate has baseline sequence evidence')
        add('grp.asm_fragmented', '1 if the genome has more than one contig (a missing locus may sit on an unassembled contig)')
        new FeatureContract(LOCUS_GROUP_ID, f)
    }
}
