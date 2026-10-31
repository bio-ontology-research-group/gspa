package gspa.integration.workflow

import com.fasterxml.jackson.databind.ObjectMapper
import gspa.integration.ClaimExtractor
import gspa.integration.EvidenceClaim
import gspa.integration.IntegrationState
import gspa.metrics.TaxonInference
import gspa.model.GenomeLayoutLoader
import gspa.ontology.AnatomicalConstraints
import gspa.ontology.ReactionGraphLoader

/**
 * Loads the reference data of {@link IntegrationWorkflowConfig} into an
 * {@link IntegrationState}. Heavyweight immutable references (GO ontology and
 * reasoner, anatomical constraints, complex classes, pathway/module databases,
 * reaction graph) are loaded once by {@link #loadShared} and copied by reference
 * into each genome's state by {@link #wire}; everything genome-specific (taxon
 * checker, operons, gaps, orthogroups, layout, essential profile) is built per
 * genome so parallel genomes cannot affect one another.
 *
 * <p>Configuration names are resolved on the {@link IntegrationWorkflowConfig}
 * ({@code c}) through {@code propertyMissing}, so the loading code reads like the
 * option it implements. A reference that was requested but cannot be loaded
 * raises an error; it is never reported as an absent prior.</p>
 */
class IntegrationReferenceLoader {
    final IntegrationWorkflowConfig c
    Closure<Void> info = { String m -> }

    IntegrationReferenceLoader(IntegrationWorkflowConfig config) {
        this.c = config
    }

    def propertyMissing(String name) { c.getProperty(name) }

    /** Progress lines (the CLI prints them; library callers usually log or ignore them). */
    private void println(Object message) { info.call(String.valueOf(message)) }

    /** Heavyweight references only, as a state that per-genome states are wired from. */
    IntegrationState loadShared() {
        def sharedState = new IntegrationState(new gspa.model.Genome(id: 'shared_reference'))
        wireReferenceData(sharedState, [])
        sharedState
    }

    /** Per-genome state: shared references by reference, genome-specific inputs freshly loaded. */
    void wire(IntegrationState state, List<EvidenceClaim> claims, IntegrationState shared) {
        wireReferenceData(state, claims, shared)
    }

    /**
     * Wire optional reference data into the IntegrationState. Each piece is
     * independent: if a file isn't provided, the corresponding prior
     * no-ops.
     */
    private void wireReferenceData(IntegrationState state, List<EvidenceClaim> claims = [], IntegrationState shared = null) {
        if (shared != null) {
            state.goOntology = shared.goOntology
            state.goReasoner = shared.goReasoner
            state.anatomicalConstraints = shared.anatomicalConstraints
            state.complexClassification = new LinkedHashMap(shared.complexClassification)
            state.pathwayDatabase = shared.pathwayDatabase
            state.moduleDefinitions = shared.moduleDefinitions
            state.reactionGraph = shared.reactionGraph
        }
        if (shared == null) {
            // GO ontology + ELK reasoner
            if (goOwlFile != null) {
                try {
                    def goOntology = new gspa.ontology.GoOntology()
                    goOntology.loadOwl(goOwlFile)
                    state.goOntology = goOntology
                    if (!lite) {
                        def reasoner = new gspa.ontology.GoReasoner(goOntology)
                        if (reasonerCacheDir != null) {
                            reasoner.cacheDir = reasonerCacheDir
                        }
                        reasoner.initialize()
                        state.goReasoner = reasoner
                    }
                    println "  Loaded GO ontology from ${goOwlFile}"
                } catch (Exception e) {
                    throw new IllegalStateException("Failed to load GO ontology ${goOwlFile}: ${e.message}", e)
                }
            }

            // Bundled anatomical constraints & complex terms
            try {
                state.anatomicalConstraints = AnatomicalConstraints.loadBundled()
                println "  Anatomical constraints: ${state.anatomicalConstraints.constrainedGoTermCount()} GO terms, ${state.anatomicalConstraints.constrainedAnatomyCount()} anatomical structures"
            } catch (Exception e) {
                throw new IllegalStateException("Failed to load bundled anatomical constraints: ${e.message}", e)
            }

            try {
                if (complexTermsFile != null && complexTermsFile.exists()) {
                    complexTermsFile.eachLine { line ->
                        if (line.startsWith('GO_term') || line.startsWith('#') || line.trim().isEmpty()) return
                        String[] f = line.split('\t')
                        if (f.length >= 2) state.complexClassification[f[0].trim()] = f[1].trim()
                    }
                } else {
                    state.loadBundledComplexTerms()
                }
                println "  Complex classification: ${state.complexClassification.size()} GO complex terms"
            } catch (Exception e) {
                throw new IllegalStateException("Failed to load complex classification: ${e.message}", e)
            }

        }

        // Essential functions
        if (essentialFunctionsFile != null) {
            try {
                state.essentialFunctions = gspa.config.EssentialFunctions
                    .loadFromTsv(essentialFunctionsFile, essentialProfile)
                println "  Essential functions: ${state.essentialFunctions.goTerms.size()} GO terms"
            } catch (Exception e) {
                throw new IllegalStateException("Failed to load essential functions ${essentialFunctionsFile}: ${e.message}", e)
            }
        } else if (essentialProfile != null) {
            try {
                state.essentialFunctions = gspa.config.EssentialFunctions.loadPreset(essentialProfile)
                println "  Essential functions (preset=${essentialProfile}): ${state.essentialFunctions.goTerms.size()} GO terms"
            } catch (Exception e) {
                throw new IllegalStateException("Failed to load essential preset '${essentialProfile}': ${e.message}", e)
            }
        }

        if (shared == null) {
            // Pathway database — supports stacking multiple sources (KEGG main +
            // KEGG Modules + MetaCyc + BioCyc) so a bacterial genome's operons
            // get reasonable enrichment hits at each granularity.
            if ((pathwaysFile != null && ec2goFile != null) || (moduleFiles != null && !moduleFiles.isEmpty())) {
                try {
                    if (pathwaysFile != null && ec2goFile != null) {
                        state.pathwayDatabase = gspa.ontology.PathwayLoader.loadPathways(
                            pathwaysFile,
                            gspa.ontology.PathwayLoader.loadEc2Go(ec2goFile),
                        )
                        println "  Pathways: ${state.pathwayDatabase.pathways.size()} loaded from ${pathwaysFile.name}"
                    } else {
                        state.pathwayDatabase = new gspa.ontology.PathwayDatabase()
                    }
                    for (File mod : (moduleFiles ?: [])) {
                        if (mod != null && mod.exists()) {
                            gspa.ontology.PathwayLoader.loadPathwaysInto(state.pathwayDatabase, mod)
                            println "  Pathways: now ${state.pathwayDatabase.pathways.size()} (added ${mod.name})"
                        }
                    }
                } catch (Exception e) {
                    throw new IllegalStateException("Failed to load pathway database: ${e.message}", e)
                }
            }
            if (moduleDefsFile != null && moduleDefsFile.exists()) {
                try {
                    Map<String, String> ec2goMap = state.pathwayDatabase?.ec2go ?: (ec2goFile != null ? gspa.ontology.PathwayLoader.loadEc2Go(ec2goFile) : [:])
                    state.moduleDefinitions = gspa.ontology.ModuleDefinitions.load(moduleDefsFile, ec2goMap)
                    state.moduleDefinitions.excludeOrganisms.addAll(moduleExcludeOrganisms ?: [])
                    println "  Modules: ${state.moduleDefinitions.modules.size()} KEGG module definitions from ${moduleDefsFile.name}"
                } catch (Exception e) {
                    throw new IllegalStateException("Failed to load module definitions ${moduleDefsFile}: ${e.message}", e)
                }
            }

        }

        // SAT consistency checker: wire when --taxonomy, --taxon, --infer-taxon, or --enforce-stage2 is given.
        boolean wantSat = (taxonomyFile != null) || (organismTaxon != null && !organismTaxon.trim().isEmpty()) || inferTaxon || enforceStage2
        if (wantSat) {
            try {
                def taxonConstraints = new gspa.ontology.TaxonConstraints()
                if (taxonConstraintsFile != null && taxonConstraintsFile.exists()) {
                    if (taxonConstraintsFile.name.toLowerCase().endsWith('.obo')) {
                        taxonConstraints.loadFromObo(taxonConstraintsFile)
                    } else {
                        taxonConstraints.loadFromTsv(taxonConstraintsFile)
                    }
                } else {
                    File bundledTc = classpathToTempFile('/taxon-constraints/go-taxon-constraints.tsv', 'go-taxon-constraints', '.tsv')
                    if (bundledTc != null && bundledTc.exists()) {
                        taxonConstraints.loadFromTsv(bundledTc)
                    } else if (state.goOntology != null) {
                        taxonConstraints.loadFromGoOntology(state.goOntology)
                    }
                }
                boolean wantAnat = enabledPriors == null || enabledPriors.trim().toLowerCase() == 'all' ||
                    enabledPriors.split(',').collect { it.trim().toLowerCase() }.contains('anatomical')
                if (state.anatomicalConstraints != null && wantAnat) {
                    state.anatomicalConstraints.augmentTaxonConstraints(taxonConstraints)
                }
                if (taxonConstraints.constrainedTermCount() > 0) {
                    def checker = new gspa.ontology.SatConsistencyChecker(taxonConstraints)
                    if (taxonomyFile != null && taxonomyFile.exists()) {
                        // Detect 3-column TSV (Term, Relationship, Other) vs 2-column legacy
                        String firstData = taxonomyFile.readLines().find { !it.trim().isEmpty() && !it.startsWith('#') } ?: ''
                        if (firstData.split('\t').length >= 3) {
                            checker.loadTaxonomyTsv(taxonomyFile)
                        } else {
                            checker.loadTaxonomyHierarchy(taxonomyFile)
                        }
                    } else {
                        File bundledHier = classpathToTempFile('/taxon-constraints/ncbi-taxon-hierarchy.tsv', 'ncbi-taxon-hierarchy', '.tsv')
                        if (bundledHier != null && bundledHier.exists()) {
                            checker.loadTaxonomyTsv(bundledHier)
                        }
                    }
                    String resolvedTax = resolveOrganismTaxon(organismTaxon)
                    if (resolvedTax) {
                        checker.organismTaxon = resolvedTax
                        println "  Asserted organism taxon: ${resolvedTax}"
                    } else if (inferTaxon && claims) {
                        Map<String, Double> termScore = [:].withDefault { 0.0d }
                        for (EvidenceClaim c : claims) {
                            if (c.functionType == 'GO' && c.calibratedProb > termScore[c.functionId]) {
                                termScore[c.functionId] = c.calibratedProb
                            }
                        }
                        def inf = new TaxonInference(checker: checker, goOntology: state.goOntology, minScore: 0.5d).infer(termScore)
                        if (inf?.taxon) {
                            checker.organismTaxon = inf.taxon
                            println "  Inferred organism taxon: ${inf.taxon} (${inf.label})"
                        }
                    }
                    state.satConsistencyChecker = checker
                    println "  SAT checker: ${taxonConstraints.onlyInTaxon.size()} only_in + ${taxonConstraints.neverInTaxon.size()} never_in (organismTaxon=${checker.organismTaxon ?: 'none'})"
                }
            } catch (Exception e) {
                throw new IllegalStateException("Failed to wire SAT consistency checker: ${e.message}", e)
            }
        } else {
            println "  SAT checker: disabled (pass --taxon, --infer-taxon, or --taxonomy to enable)"
        }

        // Operons TSV: one line per operon, tab-separated protein IDs.
        if (operonsFile != null && operonsFile.exists()) {
            List<List<String>> operons = []
            operonsFile.eachLine { line ->
                line = line.trim()
                if (line.isEmpty() || line.startsWith('#')) return
                operons << (line.split('\t') as List)
            }
            state.operons = operons
            println "  Operons: ${operons.size()}"
        }

        // Metabolic gaps JSONL
        if (gapsFile != null && gapsFile.exists()) {
            def mapper = new ObjectMapper()
            List gaps = []
            gapsFile.eachLine { line ->
                line = line.trim()
                if (line.isEmpty() || line.startsWith('#')) return
                Map rec = mapper.readValue(line, Map)
                gaps << new gspa.integration.MetabolicGap(
                    pathwayId: rec.pathway_id as String,
                    reactionId: rec.reaction_id as String,
                    ecNumber: rec.ec_number as String,
                    goTerm: rec.go_term as String,
                    gapseqGuessed: (rec.gapseq_guessed ?: false) as boolean,
                )
            }
            state.metabolicGaps = gaps
            state.metabolicGapsAvailable = true
            println "  Metabolic gaps: ${gaps.size()}"
        }

        // Phase 11 cross-genome: orthogroup membership + consensus.
        if (orthogroupsFile != null && orthogroupsFile.exists()) {
            Map<String, String> orthogroups = new LinkedHashMap<>()
            orthogroupsFile.eachLine { line ->
                line = line.trim()
                if (line.isEmpty() || line.startsWith('#')) return
                String[] parts = line.split('\t')
                if (parts.length >= 2) {
                    // Strip optional "namespace:" prefix from the rep so the
                    // orthogroup ID matches the bare-accession format used by
                    // build_catalog.py and ReactionLocusCatalog.
                    String rep = parts[1]
                    int colon = rep.indexOf(':')
                    if (colon >= 0) rep = rep.substring(colon + 1)
                    orthogroups[parts[0]] = rep
                }
            }
            state.orthogroupMap = orthogroups
            println "  Orthogroups: ${orthogroups.size()} protein assignments"
        }
        if (clusterConsensusFile != null && clusterConsensusFile.exists()) {
            Map<String, Double> consensus = new LinkedHashMap<>()
            clusterConsensusFile.eachLine { line ->
                line = line.trim()
                if (line.isEmpty() || line.startsWith('#')) return
                String[] parts = line.split('\t')
                if (parts.length >= 4) {
                    String key = "${parts[0]}|${parts[1]}|${parts[2]}".toString()
                    try {
                        consensus[key] = Double.parseDouble(parts[3])
                    } catch (NumberFormatException ignored) {}
                }
            }
            state.orthogroupConsensus = consensus
            println "  Cluster consensus entries: ${consensus.size()}"
        }

        // Phase 12 RLGC: reaction graph + genome layout.
        if (shared == null && reactionGraphFile != null && reactionGraphFile.exists()) {
            try {
                state.reactionGraph = ReactionGraphLoader.load(
                    reactionGraphFile, diffusionMetsFile, rlcCurrencyPercentile,
                    reactionEcAliasesFile)
                println "  Reaction graph: ${state.reactionGraph.reactions.size()} reactions, ${state.reactionGraph.currencyMetabolites.size()} currency metabolites, ${state.reactionGraph.ecToReactions.size()} ECs"
            } catch (Exception e) {
                throw new IllegalStateException("Failed to load reaction graph ${reactionGraphFile}: ${e.message}", e)
            }
        }
        if (genomeLayoutFile != null && genomeLayoutFile.exists()) {
            try {
                state.genomeLayout = GenomeLayoutLoader.load(genomeLayoutFile)
                println "  Genome layout: ${state.genomeLayout.size()} loci across ${state.genomeLayout.byContig.size()} contig(s)"
            } catch (Exception e) {
                throw new IllegalStateException("Failed to load genome layout ${genomeLayoutFile}: ${e.message}", e)
            }
        }
    }
    private static String resolveOrganismTaxon(String raw) {
        if (!raw) return null
        raw = raw.trim()
        switch (raw.toLowerCase()) {
            case ['bacteria', 'bacterium', 'bacterial']: return 'NCBITaxon_2'
            case ['archaea', 'archaeal']:                return 'NCBITaxon_2157'
            case ['eukaryote', 'eukaryota', 'eukaryotic']: return 'NCBITaxon_2759'
            case ['virus', 'viruses', 'viral']:          return 'NCBITaxon_10239'
        }
        if (raw.startsWith('NCBITaxon')) return raw.replace(':', '_')
        if (raw ==~ /\d+/) return "NCBITaxon_${raw}"
        return raw
    }

    private File classpathToTempFile(String resourcePath, String prefix, String suffix) {
        def stream = getClass().getResourceAsStream(resourcePath)
        if (stream == null) return null
        try {
            File tmp = File.createTempFile(prefix, suffix)
            tmp.deleteOnExit()
            tmp.withOutputStream { out -> stream.transferTo(out) }
            return tmp
        } finally {
            stream.close()
        }
    }
}
