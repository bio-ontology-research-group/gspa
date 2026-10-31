package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.CalibrationTable
import gspa.integration.ClaimExtractor
import gspa.integration.contextbridge.CohortContextBuilder
import gspa.integration.contextbridge.CohortStage
import gspa.integration.learned.ExternalContext
import spock.lang.Specification
import spock.lang.TempDir

/**
 * The cohort stage of the production workflow, on SYNTHETIC genomes: two-pass batch (stage-0 baseline
 * of every genome, one joint inference, held-out features per genome), reaction ambiguity reported
 * rather than summed, conservative handling of unknown metadata, honest absence calls, and
 * order/duplicate invariance. These check mechanisms, not accuracy.
 */
class CohortWorkflowSpec extends Specification {
    @TempDir File dir
    ContextFixtures fx
    final ObjectMapper mapper = new ObjectMapper()
    static final String KEY = 'p3|EC|EC:2.2.2.2'

    def setup() { fx = new ContextFixtures(dir) }

    private Map<String, Map<String, File>> genomes(List<String> ids, Map<String, Map> options = [:]) {
        ids.collectEntries { [it, fx.genome(it, options[it] ?: [:])] }
    }

    private static Map<String, Map> cladesFor(List<String> ids) {
        ids.collectEntries { id -> [id, [clade: 'clade_' + id, strainCluster: 'strain_' + id]] }
    }

    /** Pass 1 exactly as the CLI does it: one joint inference over all genomes. */
    private CohortContextBuilder.Result cohortOf(List<String> ids, Map<String, Map<String, File>> files, IntegrationWorkflowConfig base,
                                                 List<String> order = null) {
        def shared = IntegrationWorkflow.loadSharedReferences(base)
        def workers = (order ?: ids).collect { fx.worker(base, it, files[it], new File(dir, "${it}.tsv")) }
        new CohortStage().run(base, CohortStage.targetsOf(workers), shared.reactionGraph)
    }

    private IntegrationWorkflow.Result runGenome(String id, Map<String, File> files, IntegrationWorkflowConfig base,
                                                 CohortContextBuilder.Result cohort) {
        def wc = fx.worker(base, id, files, new File(dir, "${id}.out.tsv"))
        def wf = new IntegrationWorkflow(wc)
        wf.sharedReferences = IntegrationWorkflow.loadSharedReferences(base)
        wf.batchCohort = cohort
        wf.run(new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(files.claims))
    }

    private ExternalContext.Weighted weighted(CohortContextBuilder.Result r, String genome, String key = KEY) {
        r.contexts.get(genome)?.cohort?.get(key)?.weighted
    }

    def 'a batch has one snapshot, one joint inference and a held-out context per genome, from raw EC claims with no GO step'() {
        given:
        def ids = ['g1', 'g2', 'g3', 'g4']
        def files = genomes(ids)
        def base = fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(cladesFor(ids)))

        when:
        def r = cohortOf(ids, files, base)

        then: 'the baseline is the stage-0 raw snapshot and every genome has a report'
        r.snapshotId.startsWith('stage0_raw_max_score:')
        r.reports.keySet() == ids as Set
        r.reports.g1.baseline.label == 'stage0_raw_max_score'
        r.reports.g1.baseline.note.contains('not the legacy prior baseline')

        and: 'protein ids collide across genomes, yet every genome has its own features (namespaced in memory)'
        ids.every { weighted(r, it) != null }
        r.reports.g1.held_out_genomes == ['g1']
        r.reports.g1.n_peer_genomes == 3

        and: 'three independent clades with real metadata: adjusted outputs exist, support is high, weak candidate kept'
        def w = weighted(r, 'g1')
        w.effectiveN == 3.0d
        w.adjustedAvailable
        w.supportFractionLo90 != null && w.supportFractionHi90 != null
        w.supportFraction > 0.8d
        w.ecDirect
        w.layoutAvailable
        w.familyKnown
        w.roleEvidenceMissing && w.roleEvidence.isEmpty()
        r.contexts.g1.cohort[KEY].reactionId in ['r2', 'r2b']
    }

    def 'an EC that maps to several reactions is one candidate: alternatives are reported with their own support and never summed'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def r = cohortOf(ids, genomes(ids), fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(cladesFor(ids))))
        def co = r.contexts.g1.cohort[KEY]

        expect: 'one entry for the key; the other reaction is an alternative with support and effective N of its own'
        co.reactionAlternatives.size() == 1
        [co.reactionId, co.reactionAlternatives[0].reactionId].sort() == ['r2', 'r2b']
        co.reactionAlternatives[0].effectiveN != null
        and: 'the candidate knows it is ambiguous and how many reactions its EC maps to'
        co.weighted.ecReactionCount == 2.0d
        and: 'the report lists the ambiguity with the representative rule'
        def amb = r.reports.g1.reaction_ambiguity.find { it.function_key == KEY }
        amb.alternatives.size() == 1
        amb.representative_rule.contains('highest weighted support')
        and: 'the effective N is the independent-unit weight, not doubled by the second reaction'
        co.weighted.effectiveN == 2.0d
    }

    def 'independent curated edges and supplied role evidence are read strictly and reach the features'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def files = genomes(ids)
        File edges = new File(dir, 'edges.tsv')
        edges.text = 'reaction_a\treaction_b\tsource\tdirected\nr2\tr4\tcurated-test-edge\tfalse\n'
        File roles = new File(dir, 'roles.tsv')
        roles.text = 'genome_id\tprotein_id\tsource\trole\tscore\ng1\tp3\tdomain-test\tTRANSPORTER\t0.8\n'
        def base = fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(cladesFor(ids)))

        when:
        def plain = weighted(cohortOf(ids, files, base), 'g1')
        def rich = weighted(cohortOf(ids, files, fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(cladesFor(ids)),
            cohortInteractionEdgesFile: edges, cohortRoleEvidenceFile: roles)), 'g1')

        then: 'the curated edge makes the unrelated reaction part of the signature and p6 an anchor'
        plain.localAnchorCurated == 0.0d
        rich.localAnchorCurated == 1.0d
        and: 'supplied protein-intrinsic role evidence is raw, with its own availability; the direct EC call adds none'
        plain.roleEvidenceMissing && plain.roleEvidence.isEmpty()
        !rich.roleEvidenceMissing
        rich.roleEvidence == [(gspa.integration.learned.CandidateRole.TRANSPORTER): 0.8d]
    }

    def 'malformed side inputs fail with the line, never skipped'() {
        given:
        def ids = ['g1', 'g2']
        File bad = new File(dir, 'bad.tsv')
        bad.text = text

        when:
        cohortOf(ids, genomes(ids), fx.config(cohortInference: true, "${option}": bad))

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains(fragment)

        where:
        option                         | text                                                              | fragment
        'cohortRoleEvidenceFile'       | 'genome_id\tprotein_id\tsource\trole\tscore\ng1\tp3\ts\tCATALYTIC\t1.5\n'   | 'score must be in [0,1]'
        'cohortRoleEvidenceFile'       | 'genome_id\tprotein_id\tsource\trole\tscore\ng1\tp3\ts\tWIZARD\t0.5\n'      | 'unknown role'
        'cohortInteractionEdgesFile'   | 'reaction_a\treaction_b\tsource\tdirected\nr1\tr2\ts\tmaybe\n'              | 'directed must be true or false'
        'cohortReactionDirectionsFile' | 'reaction_id\tdirection\nr1\tSIDEWAYS\n'                                    | 'unknown direction'
        'cohortOrthogroupsFile'        | 'genome_id\tprotein_id\torthogroup_id\ng1\tp3\tOG\tsub\tbig\n'              | 'is not a number'
    }

    def 'genome order and reference order do not change any feature'() {
        given:
        def ids = ['g1', 'g2', 'g3', 'g4']
        def files = genomes(ids)
        def base = fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(cladesFor(ids)))

        when:
        def a = cohortOf(ids, files, base)
        def b = cohortOf(ids, files, base, ids.reverse())

        then:
        a.snapshotId == b.snapshotId
        ids.every { mapper.writeValueAsString(a.contexts[it]) == mapper.writeValueAsString(b.contexts[it]) }
    }

    def 'a duplicated strain adds no independent evidence'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def meta = cladesFor(ids)
        def withDup = ids + ['g2dup']
        def files = genomes(withDup)
        def metaDup = meta + [g2dup: [clade: 'clade_g2', strainCluster: 'strain_g2']]

        when:
        def plain = cohortOf(ids, files, fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(meta, 'a.json')))
        def dup = cohortOf(withDup, files, fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(metaDup, 'b.json')))
        def a = weighted(plain, 'g1'), b = weighted(dup, 'g1')

        then:
        ['effectiveN', 'supportFraction', 'supportFractionLo90', 'supportFractionHi90', 'fracSupported',
         'enrichmentLog', 'complementarity', 'supportMinusBg'].every { a[it] == b[it] }
        a.effectiveN == 2.0d
    }

    def 'unknown metadata is read conservatively: one unresolved block, adjusted outputs withheld, never invented clades'() {
        given:
        def ids = ['g1', 'g2', 'g3', 'g4']
        def files = genomes(ids)

        when: 'only the target has metadata; its three peers are unknown'
        def r = cohortOf(ids, files, fx.config(cohortInference: true,
            cohortMetadataFile: fx.metadata([g1: [clade: 'A', strainCluster: 'S1']])))
        def w = weighted(r, 'g1')

        then: 'three unknown genomes are ONE block of one clade weight, and nothing adjusted is reported'
        w.effectiveN <= 1.0d
        !w.adjustedAvailable
        w.supportFractionLo90 == null && w.supportFractionHi90 == null && w.permP == null
        w.phylogenyKnownFraction == 0.0d
        r.reports.g1.metadata.genomes_without_metadata.size() == 3
        and: 'the descriptive unadjusted share stays available'
        w.supportFractionRaw != null

        when: 'no metadata at all: every genome shares the unresolved block, so nobody has an independent peer'
        def none = cohortOf(ids, files, fx.config(cohortInference: true))

        then: 'local context survives; no independent recurrence is invented'
        ids.every { none.contexts[it].cohort.values().every { it.localOnly && it.weighted.effectiveN == null } }
        ids.every { none.reports[it].no_peers }
    }

    def 'a single run with reference genomes equals the same genome inside a batch (metadata present or absent)'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def files = genomes(ids)
        File manifest = new File(dir, 'refs.tsv')
        manifest.text = 'genome_id\tclaims_path\tlayout_path\torthogroups_path\n' +
            ['g2', 'g3'].collect { "${it}\t${files[it].claims}\t${files[it].layout}\t${files[it].orthogroups}" }.join('\n') + '\n'

        when:
        def batch = cohortOf(ids, files, fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(cladesFor(ids))))
        def singleBase = fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(cladesFor(ids)), cohortManifestFile: manifest)
        def wc = fx.worker(singleBase, 'g1', files.g1, new File(dir, 'single.tsv'))
        def single = new IntegrationWorkflow(wc).run(new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(files.g1.claims))

        then:
        single.audit.report().stages.find { it.name == 'cohort_context' }.route == 'single_run'
        mapper.writeValueAsString(single.context.cohort) == mapper.writeValueAsString(batch.contexts.g1.cohort)

        when: 'the same without metadata: both routes agree that there is no independent peer'
        def batchNone = cohortOf(ids, files, fx.config(cohortInference: true))
        def singleNone = new IntegrationWorkflow(fx.worker(fx.config(cohortInference: true, cohortManifestFile: manifest), 'g1', files.g1,
            new File(dir, 'single2.tsv'))).run(new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(files.g1.claims))

        then:
        batchNone.contexts.g1.cohort.values().every { it.localOnly }
        mapper.writeValueAsString(singleNone.context.cohort) == mapper.writeValueAsString(batchNone.contexts.g1.cohort)
        singleNone.cohortReport.no_peers
    }

    def 'absence is certified only for a complete assembly with an all-protein family search'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def files = genomes(ids, [g3: [candidate: false, withNeighbours: false]])
        def meta = { boolean certified ->
            fx.metadata(cladesFor(ids) + [g3: [clade: 'clade_g3', strainCluster: 'strain_g3', assemblyStatus: 'COMPLETE',
                familySearchComplete: certified, contigs: [c1: [length: 9000, closed: true]]]], "m-${certified}.json")
        }

        when:
        def sure = weighted(cohortOf(ids, files, fx.config(cohortInference: true, cohortMetadataFile: meta(true))), 'g1')
        def unsure = weighted(cohortOf(ids, files, fx.config(cohortInference: true, cohortMetadataFile: meta(false))), 'g1')

        then:
        sure.fracTrueAbsent > 0d
        sure.fracUnobservedAbsent == 0d
        unsure.fracTrueAbsent == 0d
        unsure.fracUnobservedAbsent > 0d
    }

    def 'a genome without peers retains local features while recurrence remains unavailable'() {
        given:
        def r = cohortOf(['g1'], genomes(['g1']), fx.config(cohortInference: true))

        expect:
        r.contexts.g1.cohort.values().every { it.localOnly && it.weighted.effectiveN == null }
        r.reports.g1.no_peers
        r.reports.g1.candidates.keys_without_peers == 4
        r.reports.g1.candidates.keys_with_cohort_features == 0
        r.reports.g1.candidates.keys_with_local_features == 3 // isolated r4 has no testable signature
    }

    def 'metadata is validated strictly'() {
        given:
        def ids = ['g1', 'g2']
        def files = genomes(ids)
        File bad = new File(dir, 'bad.json')
        bad.text = text

        when:
        cohortOf(ids, files, fx.config(cohortInference: true, cohortMetadataFile: bad))

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains(fragment)

        where:
        text                                                                                         | fragment
        '{"genomes":{"g1":{"cladd":"x"}}}'                                                           | 'not valid'
        '{"genomes":{"g1":{"assemblyStatus":"FINISHED"}}}'                                           | 'assemblyStatus'
        '{"genomes":{"nope":{"clade":"x"}}}'                                                         | "not in the cohort"
        '{"genomes":{"g1":{"familySearchComplete":true,"familySearchCoverage":0.5}}}'                | 'contradicts'
        '{"genomes":{"g1":{"assemblyStatus":"COMPLETE","contigs":{"c1":{"closed":false}}}}}'        | 'COMPLETE contradicts'
        '{"genomes":{"g1":{"clade":" "}}}'                                                           | 'blank'
    }

    def 'claims, layout and orthogroups must share the FAA ids: a mismatch fails loudly instead of joining zero rows'() {
        given:
        def files = genomes(['g1', 'g2'])
        files.g2.layout.text = 'protein_id\tcontig\tstart\tend\tstrand\nGENE_1\tc1\t1000\t1900\t+\n'

        when:
        cohortOf(['g1', 'g2'], files, fx.config(cohortInference: true))

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('not in its layout')
        e.message.contains('FAA sequence ids')
    }

    /** g2 loses the layout row of its candidate p3 (a real protein with no usable coordinate). */
    private Map<String, Map<String, File>> withUnplacedCandidate(List<String> ids, String genome = 'g2') {
        def files = genomes(ids)
        files[genome].layout.text = files[genome].layout.readLines().findAll { !it.startsWith('p3\t') }.join('\n') + '\n'
        files
    }

    def 'a declared unplaced protein is accepted: family present, context unknown, no coordinate invented, other genomes unaffected'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def files = withUnplacedCandidate(ids)
        def meta = cladesFor(ids) + [g2: [clade: 'clade_g2', strainCluster: 'strain_g2', unplacedProteins: ['p3']]]
        def base = fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(meta))

        when:
        def r = cohortOf(ids, files, base)

        then: 'the unplaced genome still gets a context; its own local context is unavailable, not zero'
        def own = weighted(r, 'g2')
        own != null
        !own.layoutAvailable
        own.localAnchorCount == null
        own.familyKnown
        and: 'the raw candidate and its family are retained'
        r.contexts.g2.cohort[KEY] != null
        r.reports.g2.orthogroup_input.proteins_with_family == 3
        and: 'the unplaced accounting is reported, not silent'
        r.reports.g2.unplaced_proteins.declared == 1
        r.reports.g2.unplaced_proteins.with_ec_calls == 1
        r.reports.g2.unplaced_proteins.in_families == 1
        r.reports.g1.unplaced_proteins.declared == 0
        and: 'as a peer of g1 the unplaced member is present with unknown context: one of two peers'
        def peer = weighted(r, 'g1')
        peer.layoutAvailable
        peer.fracUnknownContext > 0d
        peer.fracTrueAbsent == 0d
        peer.fracNegative == 0d
    }

    def 'an unplaced family member in a complete, certified assembly is present with unknown context, never true absence'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def meta = { boolean unplaced ->
            def g3 = [clade: 'clade_g3', strainCluster: 'strain_g3', assemblyStatus: 'COMPLETE', familySearchComplete: true,
                      contigs: [c1: [length: 9000, closed: true]]]
            if (unplaced) g3.unplacedProteins = ['p3']
            fx.metadata(cladesFor(ids) + [g3: g3], "m-${unplaced}.json")
        }

        when: 'g3 truly lacks the candidate'
        def absent = weighted(cohortOf(ids, genomes(ids, [g3: [candidate: false, withNeighbours: false]]),
            fx.config(cohortInference: true, cohortMetadataFile: meta(false))), 'g1')
        and: 'g3 has the candidate but no coordinate for it'
        def unobserved = weighted(cohortOf(ids, withUnplacedCandidate(ids, 'g3'),
            fx.config(cohortInference: true, cohortMetadataFile: meta(true))), 'g1')

        then: 'complete assembly + certified search: genuine absence only in the first case'
        absent.fracTrueAbsent > 0d
        absent.fracUnknownContext == 0d
        unobserved.fracTrueAbsent == 0d
        unobserved.fracUnobservedAbsent == 0d
        unobserved.fracUnknownContext > 0d
    }

    def 'without the declaration the same genome is still rejected, and typos in the declaration are rejected'() {
        given:
        def ids = ['g1', 'g2']
        def files = withUnplacedCandidate(ids)
        File bad = new File(dir, 'unplaced.json')
        bad.text = text

        when:
        cohortOf(ids, files, fx.config(cohortInference: true, cohortMetadataFile: bad))

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains(fragment)

        where:
        text                                                                    | fragment
        '{"genomes":{"g2":{"clade":"x"}}}'                                      | 'not in its layout'
        '{"genomes":{"g2":{"unplacedProteins":["p9"]}}}'                        | 'not in its layout'
        '{"genomes":{"g2":{"unplacedProteins":["p3","p2"]}}}'                   | 'have a coordinate in the layout'
        '{"genomes":{"g2":{"unplacedProteins":["p3","p3"]}}}'                   | 'twice'
        '{"genomes":{"g2":{"unplacedProteins":["p3",""]}}}'                     | 'blank'
        '{"genomes":{"g2":{"unplacedProteins":["p3"],"unplaced":["p3"]}}}'      | 'not valid'
    }

    def 'the workflow records the cohort stage, its baseline label and the legacy scoring mode'() {
        given:
        def ids = ['g1', 'g2', 'g3']
        def files = genomes(ids)
        def base = fx.config(cohortInference: true, cohortMetadataFile: fx.metadata(cladesFor(ids)))
        def cohort = cohortOf(ids, files, base)

        when:
        def result = runGenome('g1', files.g1, base, cohort)
        def report = result.audit.report()
        def stage = report.stages.find { it.name == 'cohort_context' }

        then:
        stage.status == 'executed'
        stage.route == 'batch'
        stage.baseline == 'stage0_raw_max_score'
        stage.baseline_snapshot_id == cohort.snapshotId
        result.cohortReport.genome_id == 'g1'
        result.context.cohort.containsKey(KEY)
        and: 'no learned model was requested: the ledger says the scores are the legacy heuristic ones'
        report.scoring_mode == 'legacy_heuristic'
        report.stages.find { it.name == 'learned_integration' }.status == 'skipped'
        report.stages.find { it.name == 'learned_integration' }.reason.contains('legacy heuristic')
        report.baseline_definitions.keySet() == ['stage0_raw_max_score', 'legacy_prior_baseline'] as Set
        report.changes_measured_from == 'legacy_prior_baseline'
        and: 'the ledger file carries the same stages'
        mapper.readValue(new File(dir, 'g1.out.tsv.stages.json'), Map).stages.find { it.name == 'cohort_context' }.status == 'executed'
    }
}
