package gspa.integration.learned

/**
 * The contract columns that carry {@link ExternalContext.Weighted} cohort features, declared once
 * so the contract ({@link FeatureContract}) and the extractor ({@link ContextFeatureExtractor})
 * cannot disagree. Every nullable field belongs to an availability group; the group's indicator
 * column is 1 only when the group's lead field was supplied, and the value columns of an
 * unavailable group are 0 (missing is not absent, as everywhere in the contract).
 *
 * <p>Nothing here recomputes, rounds or thresholds a value: it is a typed copy of the cohort
 * engine's output. Counts enter as {@code log(1 + n)} so that they stay bounded; fractions and
 * signed log-ratios enter unchanged.</p>
 */
class CohortFeatureColumns {

    static class Entry {
        String column            // contract column
        String property          // ExternalContext.Weighted property
        String group             // availability group; null = available whenever Weighted is present
        boolean log1p
        String description
    }

    /** group name -> lead property (non-null means the group is available). */
    static final Map<String, String> GROUP_LEADS = [
        states        : 'fracSupported',
        support_raw   : 'supportFractionRaw',
        interval      : 'supportFractionLo90',
        unknown_sens  : 'supportUnknownAsNegative',
        bg            : 'bgRate',
        complementary : 'complementarity',
        perm          : 'permP',
        multicopy     : 'multicopyFraction',
        organism_rxn  : 'peerReactionInOrganismFraction',
        local         : 'localAnchorCount',
        local_enrich  : 'localEnrichmentLog',
    ].asImmutable()

    static final String PREFIX = 'cohort.w.'

    static final List<Entry> ENTRIES = build([
        // peer set and uncertainty (always reported with a weighted block)
        [null, 'effectiveN', true, 'log(1 + effective N): summed weight of informative independent clade/strain units'],
        [null, 'phylogenyKnownFraction', false, 'clade/strain-weighted share of peers with resolved clade metadata (1 = all known; unknown never counts as independent)'],
        [null, 'coverage', false, 'share of peer clades that were informative'],
        [null, 'supportFraction', false, 'weighted Jeffreys-smoothed share of informative peers whose family member sits in a context-positive window'],
        [null, 'peerOrganismSetKnownFraction', false, 'clade/strain-weighted share of peers with a known organism-specific reaction set'],
        [null, 'membershipWeight', false, 'soft family membership weight of the candidate protein'],
        [null, 'nFamilyAlternatives', true, 'log(1 + further family assignments of the protein); paralog/subfamily ambiguity'],
        [null, 'nReactionAlternatives', true, 'log(1 + further reaction hypotheses of the protein); alternatives are reported, not summed'],
        [null, 'ecReactionCount', true, 'log(1 + reactions the supporting EC call maps to); 0 when the candidate did not come from a direct EC call'],
        // observation states
        ['states', 'fracSupported', false, 'weighted share of peers: family member in a context-positive window'],
        ['states', 'fracNegative', false, 'weighted share of peers: family present, fully observed windows lack the signature'],
        ['states', 'fracUnknownContext', false, 'weighted share of peers: family present but the window is cut or of unknown extent'],
        ['states', 'fracTrueAbsent', false, 'weighted share of peers where absence is certified (complete assembly AND all-protein family search)'],
        ['states', 'fracUnobservedAbsent', false, 'weighted share of peers with no member but absence NOT certified (unknown, not absence)'],
        ['states', 'fracPathwayWithoutFamily', false, 'weighted share of certified-absent peers that carry the pathway elsewhere'],
        // raw (unsmoothed) support and sensitivity to unobservable windows
        ['support_raw', 'supportFractionRaw', false, 'unsmoothed weighted support fraction among informative peers'],
        ['unknown_sens', 'supportUnknownAsNegative', false, 'support fraction if unknown-context windows were negative'],
        ['unknown_sens', 'supportUnknownAsSupport', false, 'support fraction if unknown-context windows were support'],
        // clade-adjusted uncertainty
        ['interval', 'supportFractionLo90', false, 'approximate 90% lower bound on the support fraction (only with clade metadata)'],
        ['interval', 'supportFractionHi90', false, 'approximate 90% upper bound on the support fraction (only with clade metadata)'],
        // matched background
        ['bg', 'bgRate', false, 'rate of the signature in matched same-genome background windows (never total genome presence)'],
        ['bg', 'enrichmentLog', false, 'log ratio of support to matched background'],
        ['bg', 'supportMinusBg', false, 'support fraction minus matched background rate'],
        // complementary partial contexts
        ['complementary', 'complementarity', false, 'joint minus mean single coverage of the reaction signature across independent lineages'],
        ['complementary', 'complementaryUnannotated', false, 'weighted share of signature reactions with no confident baseline candidate in a peer; UNANNOTATED label only, never evidence of loss'],
        ['complementary', 'jointCoverage', false, 'share of signature reactions jointly covered by independent units next to the family'],
        // permutation null
        ['perm', 'permP', false, 'add-one permutation p of the support statistic against matched windows (only with clade metadata)'],
        ['perm', 'permExcess', false, 'observed minus permutation-null mean of the support statistic'],
        // copy number
        ['multicopy', 'multicopyFraction', false, 'weighted share of informative peers with several family members (only with certified family search)'],
        // organism-specific plausibility
        ['organism_rxn', 'peerReactionInOrganismFraction', false, 'weighted share of peers whose own plausible reaction set contains the reaction'],
        // the candidate's own window
        ['local', 'localAnchorWeight', false, 'sum of hop-decayed baseline confidences of anchor reactions in the candidate window'],
        ['local', 'localAnchorCount', true, 'log(1 + distinct anchor reactions in the candidate window)'],
        ['local', 'localAnchorCurated', true, 'log(1 + anchors linked through curated interaction edges)'],
        ['local', 'localAnchorUpstream', true, 'log(1 + upstream anchors, respecting reaction direction)'],
        ['local', 'localAnchorDownstream', true, 'log(1 + downstream anchors, respecting reaction direction)'],
        ['local', 'localContextPositive', false, 'candidate window is context-positive (enough anchors)'],
        ['local', 'localWindowTruncated', false, 'candidate window is cut by a known contig end'],
        ['local', 'localWindowUnknown', false, 'candidate window extent is unknown'],
        ['local_enrich', 'localEnrichmentLog', false, 'log ratio of the window signature to matched same-genome windows'],
    ])

    /** Boolean flags of the block, each its own column. */
    static final List<List> FLAGS = [
        ['insufficient_peers', 'insufficientPeers', 'peers exist but none was informative (effective N is 0)'],
        ['adjusted_available', 'adjustedAvailable', 'clade/strain-adjusted outputs exist (metadata known); 0 = intervals and permutation p are withheld'],
        ['layout_available', 'layoutAvailable', 'gene layout of the candidate genome was available to the cohort engine'],
        ['ec_direct', 'ecDirect', 'the candidate arose from a direct EC -> reaction expansion (no EC -> GO step)'],
        ['family_known', 'familyKnown', 'the protein has a family assignment'],
        ['role_evidence_missing', 'roleEvidenceMissing', 'no protein-intrinsic role evidence at all (missing, NOT non-catalytic)'],
    ]

    /** Raw role evidence columns: value and a reported flag per role. */
    static final List<CandidateRole> EVIDENCE_ROLES = [CandidateRole.CATALYTIC, CandidateRole.ACCESSORY,
        CandidateRole.TRANSPORTER, CandidateRole.REGULATOR, CandidateRole.COMPLEX_SUBUNIT].asImmutable()

    static final String AVAILABLE = PREFIX + 'available'

    /** The same engine computes local evidence even when no independent peers exist. */
    static final Set<String> LOCAL_GROUPS = ['local', 'local_enrich'] as Set
    static final Set<String> LOCAL_PROPERTIES = ['membershipWeight', 'nFamilyAlternatives',
        'nReactionAlternatives', 'ecReactionCount'] as Set
    static final Set<String> LOCAL_FLAGS = ['layoutAvailable', 'ecDirect', 'familyKnown', 'roleEvidenceMissing'] as Set

    static boolean localEntry(Entry e) { e.group in LOCAL_GROUPS || e.property in LOCAL_PROPERTIES }

    static Set<String> peerProperties() {
        (ENTRIES.findAll { !localEntry(it) }*.property + ['nPeerGenomes', 'bgWindows']) as Set<String>
    }

    static String roleColumn(CandidateRole r) { PREFIX + 'role_evidence.' + r.name().toLowerCase(Locale.ROOT) }
    static String roleReportedColumn(CandidateRole r) { roleColumn(r) + '_reported' }
    static String groupColumn(String group) { PREFIX + group + '_available' }

    private static List<Entry> build(List<List> rows) {
        rows.collect { List r ->
            new Entry(group: r[0] as String, property: r[1] as String, log1p: r[2] as boolean, description: r[3] as String,
                column: PREFIX + snake(r[1] as String))
        }.asImmutable()
    }

    static String snake(String camel) {
        camel.replaceAll(/([a-z0-9])([A-Z])/, '$1_$2').toLowerCase(Locale.ROOT)
    }

    /** Declare every column of the block on a contract under {@code group}. */
    static void declare(Closure add, FeatureGroup featureGroup) {
        add(AVAILABLE, featureGroup, '1 if a joint (phylogeny-weighted, background-matched) cohort analysis supplied features for this candidate')
        GROUP_LEADS.each { String g, String lead ->
            add(groupColumn(g), featureGroup, "1 if the '${g}' feature group (${lead}) was computed for this candidate; its values are 0 otherwise")
        }
        FLAGS.each { List f -> add(PREFIX + f[0], featureGroup, f[2] as String) }
        ENTRIES.each { Entry e -> add(e.column, featureGroup, e.description + (e.group ? " [group ${e.group}]" : '')) }
        EVIDENCE_ROLES.each { CandidateRole r ->
            add(roleColumn(r), featureGroup, "raw protein-intrinsic ${r.name()} evidence in [0,1] (0 if not reported)")
            add(roleReportedColumn(r), featureGroup, "1 if raw ${r.name()} evidence was reported for this protein")
        }
    }

    /** Write the block's values for one candidate through {@code set(column, value)}. */
    static void extract(ExternalContext.Weighted w, Closure set, boolean localOnly = false) {
        if (!localOnly) set(AVAILABLE, 1)
        GROUP_LEADS.each { String g, String lead ->
            if (localOnly && !(g in LOCAL_GROUPS)) return
            if (w.getProperty(lead) != null) set(groupColumn(g), 1)
        }
        FLAGS.each { List f ->
            if ((!localOnly || f[1] in LOCAL_FLAGS) && w.getProperty(f[1] as String)) set(PREFIX + f[0], 1)
        }
        ENTRIES.each { Entry e ->
            if (localOnly && !localEntry(e)) return
            Double v = w.getProperty(e.property) as Double
            if (v == null) return
            if (e.group != null && w.getProperty(GROUP_LEADS[e.group]) == null) return
            set(e.column, e.log1p ? Math.log(1.0d + v) : v)
        }
        EVIDENCE_ROLES.each { CandidateRole r ->
            Double v = w.roleEvidence?.get(r)
            if (v != null) {
                set(roleColumn(r), v)
                set(roleReportedColumn(r), 1)
            }
        }
    }
}
