package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper

/**
 * SYNTHETIC files for the production context workflow tests: a five-reaction network, genomes with
 * the same eight protein ids each (so ids collide across genomes, as real FAA seqids do), EC claims
 * with raw scores, layouts, orthogroups and cohort metadata. They exercise mechanisms; they say
 * nothing about accuracy.
 *
 * <pre>
 * r1: a -> b   EC 1.1.1.1        p2 carries it
 * r2: b -> c   EC 2.2.2.2        the candidate p3 (weak raw score) - EC maps to TWO reactions (r2, r2b)
 * r2b: b -> f  EC 2.2.2.2
 * r3: c -> d   EC 3.3.3.3        p4 carries it
 * r4: x -> y   EC 4.4.4.4        unrelated, p6
 * </pre>
 */
class ContextFixtures {
    final File dir
    final File reactions
    final File aliases
    final ObjectMapper mapper = new ObjectMapper()

    ContextFixtures(File dir) {
        this.dir = dir
        reactions = new File(dir, 'reactions.tsv')
        reactions.text = 'id\tstoichiometry\tec_numbers\n' +
            'r1\t-1:a:0:0;1:b:0:0\t1.1.1.1\n' +
            'r2\t-1:b:0:0;1:c:0:0\t2.2.2.2\n' +
            'r2b\t-1:b:0:0;1:f:0:0\t2.2.2.2\n' +
            'r3\t-1:c:0:0;1:d:0:0\t3.3.3.3\n' +
            'r4\t-1:x:0:0;1:y:0:0\t4.4.4.4\n'
        aliases = new File(dir, 'aliases.tsv')
        aliases.text = 'MS ID\tOld MS ID\tExternal ID\tSource\n' +
            'r1\t\t1.1.1.1\tt\nr2|r2b\t\t2.2.2.2\tt\nr3\t\t3.3.3.3\tt\nr4\t\t4.4.4.4\tt\n'
    }

    /**
     * Writes claims/layout/orthogroups for one genome and returns their paths. Options: {@code candidate}
     * (default true: p3 carries a weak EC 2.2.2.2 claim and an orthogroup row), {@code family} (orthogroup of
     * p3), {@code withNeighbours} (default true), {@code score} (raw score of the candidate).
     */
    Map<String, File> genome(String id, Map o = [:]) {
        boolean candidate = o.containsKey('candidate') ? o.candidate : true
        boolean neighbours = o.containsKey('withNeighbours') ? o.withNeighbours : true
        double score = (o.score ?: 0.4d) as double
        List<Map> claims = []
        def claim = { String p, String type, String fn, double raw, String src = 'diamond' ->
            claims << [protein_id: p, function_type: type, function_id: fn, source: src, evidence_type: 'SEQUENCE_SIMILARITY',
                       raw_score: raw, calibrated_prob: raw] + (type == 'GO' ? [go_aspect: 'MF'] : [:])
        }
        if (neighbours) {
            claim('p2', 'EC', 'EC:1.1.1.1', 0.9d)
            claim('p4', 'EC', 'EC:3.3.3.3', 0.9d)
        }
        claim('p6', 'EC', 'EC:4.4.4.4', 0.9d)
        if (candidate) {
            claim('p3', 'EC', 'EC:2.2.2.2', score)
            claim('p3', 'EC', 'EC:2.2.2.2', score / 2, 'interproscan')   // a second raw claim: only the strongest counts
            claim('p3', 'GO', 'GO:0003824', 0.7d)                          // a non-EC claim: ignored by the EC baseline
        }
        File claimsFile = new File(dir, "${id}.claims.jsonl")
        claimsFile.text = claims.collect { mapper.writeValueAsString(it) }.join('\n') + '\n'

        File layout = new File(dir, "${id}.layout.tsv")
        layout.text = 'protein_id\tcontig\tstart\tend\tstrand\n' +
            (1..8).collect { "p${it}\tc1\t${it * 1000}\t${it * 1000 + 900}\t+" }.join('\n') + '\n'

        File og = new File(dir, "${id}.orthogroups.tsv")
        List<String> rows = ['p2\tOG_up', 'p4\tOG_down']
        if (candidate || o.family) rows << "p3\t${o.family ?: 'OG_cand'}"
        og.text = rows.join('\n') + '\n'
        [claims: claimsFile, layout: layout, orthogroups: og]
    }

    /** Cohort metadata JSON. {@code genomes} maps id -> Map of metadata keys. */
    File metadata(Map<String, Map> genomes, String name = 'cohort-metadata.json', Map extra = [:]) {
        File f = new File(dir, name)
        f.text = mapper.writeValueAsString(extra + [genomes: genomes])
        f
    }

    IntegrationWorkflowConfig config(Map overrides = [:]) {
        def c = new IntegrationWorkflowConfig(lite: true, enabledPriors: 'none', reactionGraphFile: reactions,
            reactionEcAliasesFile: aliases, rlcCurrencyPercentile: 100.0d)
        overrides.each { k, v -> c[k as String] = v }
        c
    }

    /** One genome's own configuration inside a batch. */
    IntegrationWorkflowConfig worker(IntegrationWorkflowConfig base, String id, Map files, File out) {
        def w = base.copy()
        w.genomeId = id
        w.claimsFile = files.claims
        w.genomeLayoutFile = files.layout
        w.orthogroupsFile = files.orthogroups
        w.outFile = out
        w
    }
}
