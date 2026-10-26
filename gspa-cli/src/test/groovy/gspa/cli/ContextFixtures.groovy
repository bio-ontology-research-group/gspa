package gspa.cli

import com.fasterxml.jackson.databind.ObjectMapper

/**
 * SYNTHETIC files for the CLI context tests (same construction as the core workflow tests): a small
 * reaction network where EC 2.2.2.2 (the weak candidate p3) maps to two reactions, genomes that share
 * their eight protein ids, EC claims with raw scores, layouts, orthogroups and cohort metadata.
 * Mechanism checks only; no accuracy claim.
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

    Map<String, File> genome(String id) {
        List<Map> claims = []
        def claim = { String p, String fn, double raw ->
            claims << [protein_id: p, function_type: 'EC', function_id: fn, source: 'diamond',
                       evidence_type: 'SEQUENCE_SIMILARITY', raw_score: raw, calibrated_prob: raw]
        }
        claim('p2', 'EC:1.1.1.1', 0.9d)
        claim('p3', 'EC:2.2.2.2', 0.4d)
        claim('p4', 'EC:3.3.3.3', 0.9d)
        claim('p6', 'EC:4.4.4.4', 0.9d)
        File claimsFile = new File(dir, "${id}.claims.jsonl")
        claimsFile.text = claims.collect { mapper.writeValueAsString(it) }.join('\n') + '\n'
        File layout = new File(dir, "${id}.layout.tsv")
        layout.text = 'protein_id\tcontig\tstart\tend\tstrand\n' +
            (1..8).collect { "p${it}\tc1\t${it * 1000}\t${it * 1000 + 900}\t+" }.join('\n') + '\n'
        File og = new File(dir, "${id}.orthogroups.tsv")
        og.text = 'p2\tOG_up\np4\tOG_down\np3\tOG_cand\n'
        [claims: claimsFile, layout: layout, orthogroups: og]
    }

    File metadata(List<String> ids, String name = 'cohort-metadata.json') {
        File f = new File(dir, name)
        f.text = mapper.writeValueAsString([genomes: ids.collectEntries { [it, [clade: 'clade_' + it, strainCluster: 'strain_' + it]] }])
        f
    }

    /** A batch manifest row with the layout and orthogroup columns filled in. */
    String row(String id, Map<String, File> files, File out) {
        [id, files.claims, out, '-', '-', '-', '-', '-', '-', '-', '-', files.layout, files.orthogroups].join('\t')
    }

    File manifest(List<String> ids, Map<String, Map<String, File>> files, String name = 'batch.tsv') {
        File f = new File(dir, name)
        f.text = ids.collect { row(it, files[it], new File(dir, "${it}.integrated.tsv")) }.join('\n') + '\n'
        f
    }

    /** Options that turn the cohort on for the shared reaction network. */
    List<String> cohortArgs(File metadata) {
        ['--lite', '--enable-priors', 'none', '--reaction-graph', reactions.path, '--reaction-ec-aliases', aliases.path,
         '--rlc-currency-pct', '100', '--cohort', '--cohort-metadata', metadata.path]
    }
}
