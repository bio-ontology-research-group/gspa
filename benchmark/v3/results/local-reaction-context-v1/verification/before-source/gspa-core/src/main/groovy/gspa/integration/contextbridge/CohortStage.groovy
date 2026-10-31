package gspa.integration.contextbridge

import gspa.integration.CalibrationTable
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.integration.cohort.CohortInput
import gspa.integration.workflow.IntegrationWorkflowConfig
import gspa.model.GenomeLayout
import gspa.model.GenomeLayoutLoader
import gspa.ontology.ReactionGraph

/**
 * The production cohort stage: turns workflow configuration plus the genomes of a run into
 * {@link CohortContextBuilder} members and builds the per-genome learned-stage context in ONE joint
 * inference. Used identically by the batch pass (all manifest genomes, then each genome's learned
 * stage consumes its own held-out features) and by a single run (the genome plus the reference
 * genomes of {@code cohortManifestFile}). With no peers the result carries no features for the genome
 * and says so; a one-genome run is routed to the cohort-free feature set, never read as negative.
 *
 * <p>Reference manifest TSV (header required): {@code genome_id  claims_path  [layout_path]  [orthogroups_path]};
 * relative paths resolve against the manifest. Reference genomes are peers only.</p>
 */
class CohortStage {

    /** One genome that is built a context for. */
    static class Target {
        String genomeId
        /** Raw claims, read lazily so a batch never holds every genome's claims at once. */
        Closure<List<EvidenceClaim>> claims
        GenomeLayout layout
        File layoutFile
        /** The genome's own {@code --orthogroups} file (protein, cluster); ignored when a cohort orthogroups file is given. */
        File orthogroupsFile
    }

    Closure<Void> info = { String m -> }

    /** Targets of a batch from the per-genome configurations: claims, layout and orthogroup files of each genome. */
    static List<Target> targetsOf(Collection<IntegrationWorkflowConfig> perGenome) {
        perGenome.collect { IntegrationWorkflowConfig wc ->
            if (wc.claimsFile == null) throw new IllegalArgumentException("genome '${wc.genomeId}' has no claims file")
            new Target(genomeId: wc.genomeId, layoutFile: wc.genomeLayoutFile, orthogroupsFile: wc.orthogroupsFile,
                claims: { -> new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(wc.claimsFile) })
        }
    }

    CohortContextBuilder.Result run(IntegrationWorkflowConfig c, List<Target> targets, ReactionGraph graph) {
        if (graph == null) throw new IllegalArgumentException('cohort inference needs the reaction graph (--reaction-graph); it was not loaded')
        List<Map> references = c.cohortManifestFile != null ? readReferenceManifest(c.cohortManifestFile) : []
        Set<String> targetIds = targets*.genomeId as Set
        references.each { Map r ->
            if (targetIds.contains(r.genomeId)) {
                throw new IllegalArgumentException("genome '${r.genomeId}' is both a target of this run and a reference genome in ${c.cohortManifestFile}")
            }
        }
        Set<String> cohortIds = (targetIds + references*.genomeId) as Set
        CohortOrthogroupReader.Rows cohortRows = c.cohortOrthogroupsFile != null ?
            CohortOrthogroupReader.readCohortFile(c.cohortOrthogroupsFile, cohortIds) : null

        def settings = new CohortContextBuilder.Settings(config: c.cohortConfig, graph: graph)
        if (c.cohortMetadataFile != null) settings.metadata = CohortMetadata.read(c.cohortMetadataFile)
        if (c.cohortRoleEvidenceFile != null) settings.roleEvidence = CohortSideInputs.readRoleEvidence(c.cohortRoleEvidenceFile)
        if (c.cohortInteractionEdgesFile != null) settings.edges = CohortSideInputs.readEdges(c.cohortInteractionEdgesFile)
        if (c.cohortReactionDirectionsFile != null) settings.directions = CohortSideInputs.readDirections(c.cohortReactionDirectionsFile)

        List<CohortContextBuilder.Member> members = []
        targets.each { Target t ->
            members << member(t.genomeId, t.claims.call(), t.layout ?: loadLayout(t.layoutFile, t.genomeId), t.orthogroupsFile,
                cohortRows, true, settings.metadata)
        }
        references.each { Map r ->
            List<EvidenceClaim> claims = new ClaimExtractor(new CalibrationTable()).readClaimsJsonl(r.claims as File)
            members << member(r.genomeId as String, claims, loadLayout(r.layout as File, r.genomeId as String),
                r.orthogroups as File, cohortRows, false, settings.metadata)
        }
        info.call("  Cohort: ${members.size()} genomes (${targets.size()} targets, ${references.size()} reference); baseline = stage-0 raw scores".toString())
        CohortContextBuilder.Result result = new CohortContextBuilder().build(members, settings)
        info.call("  Cohort baseline snapshot ${result.snapshotId}; contexts for ${result.contexts.size()} of ${targets.size()} target genomes".toString())
        result
    }

    private static GenomeLayout loadLayout(File file, String genomeId) {
        file == null ? null : GenomeLayoutLoader.load(file, genomeId)
    }

    private static CohortContextBuilder.Member member(String genomeId, List<EvidenceClaim> claims, GenomeLayout layout,
                                                       File perGenomeOrthogroups, CohortOrthogroupReader.Rows cohortRows, boolean target,
                                                       CohortMetadata metadata) {
        def baseline = RawBaselineSnapshot.forGenome(genomeId, claims)
        def m = new CohortContextBuilder.Member(genomeId: genomeId, baseline: baseline, layout: layout, target: target)
        if (cohortRows != null) {
            m.families = cohortRows.byGenome.get(genomeId) ?: []
            m.orthogroupNotes = [source: 'cohort_orthogroups_file', rows_in_file: cohortRows.rows,
                                 rows_of_genomes_outside_cohort: cohortRows.rowsOfGenomesOutsideCohort,
                                 genomes_outside_cohort: cohortRows.genomesOutsideCohort.toList(),
                                 proteins_with_family: m.families*.proteinId.toSet().size()]
        } else if (perGenomeOrthogroups != null) {
            Set<String> known = layout != null ? (layout.byProtein.keySet() as Set<String>) : (baseline.calls*.proteinId as Set<String>)
            // Declared unplaced proteins belong to the genome too; dropping their rows would turn a present family member into an absent one.
            if (layout != null) known = known + ((metadata?.genomes?.get(genomeId)?.unplacedProteins ?: []) as Set<String>)
            def rows = CohortOrthogroupReader.readPerGenomeFile(perGenomeOrthogroups, genomeId, known)
            m.families = rows.byGenome.get(genomeId) ?: []
            m.orthogroupNotes = [source: 'per_genome_orthogroups_file', rows_in_file: rows.rows,
                                 rows_of_other_genomes: rows.rowsOfProteinsOutsideGenome,
                                 proteins_with_family: m.families*.proteinId.toSet().size()]
        } else {
            m.orthogroupNotes = [source: 'none', note: 'no family assignment: cohort recurrence is unavailable for this genome']
        }
        m
    }

    static List<Map> readReferenceManifest(File file) {
        if (!file.isFile()) throw new IllegalArgumentException("cohort manifest does not exist: ${file}")
        List<Map> out = []
        Set<String> seen = new LinkedHashSet<>()
        boolean header = false
        int lineNo = 0
        file.eachLine('UTF-8') { String raw ->
            lineNo++
            if (!raw.trim() || raw.trim().startsWith('#')) return
            String[] f = raw.split('\t', -1)
            if (!header) {
                header = true
                if (f.length < 2 || f[0].trim().toLowerCase() != 'genome_id' || !(f[1].trim().toLowerCase() in ['claims_path', 'claims'])) {
                    throw new IllegalArgumentException("${file}: the first line must be the header 'genome_id<TAB>claims_path[<TAB>layout_path<TAB>orthogroups_path]'")
                }
                return
            }
            if (f.length < 2 || f.length > 4 || !f[0].trim() || !f[1].trim()) {
                throw new IllegalArgumentException("${file}:${lineNo}: expected genome_id, claims_path and at most layout_path, orthogroups_path")
            }
            String id = f[0].trim()
            if (!seen.add(id)) throw new IllegalArgumentException("${file}:${lineNo}: duplicate reference genome '${id}'")
            Closure<File> path = { String p ->
                File x = new File(p.trim())
                x.isAbsolute() ? x : new File(file.absoluteFile.parentFile, p.trim())
            }
            File claims = path(f[1])
            if (!claims.isFile()) throw new IllegalArgumentException("${file}:${lineNo}: claims file does not exist: ${claims}")
            File layout = f.length > 2 && f[2].trim() && f[2].trim() != '-' ? path(f[2]) : null
            File og = f.length > 3 && f[3].trim() && f[3].trim() != '-' ? path(f[3]) : null
            [layout, og].each { File x -> if (x != null && !x.isFile()) throw new IllegalArgumentException("${file}:${lineNo}: file does not exist: ${x}") }
            out << [genomeId: id, claims: claims, layout: layout, orthogroups: og]
        }
        if (!header) throw new IllegalArgumentException("${file} is empty")
        out
    }
}
