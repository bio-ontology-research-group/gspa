/* One joint cohort pass through the production CLI. No inference is implemented here. */

def cohortColumns() {
    ['genome_id', 'claims_path', 'out_path', 'taxon', 'essential_profile', 'operons_path',
     'gaps_path', 'max_iter', 'enforce_stage2', 'enforce_mode', 'enable_priors',
     'genome_layout_path', 'orthogroups_path', 'cluster_consensus_path', 'fba_medium', 'fba_gram']
}

def cohortBatch(rows) {
    if (!rows) throw new IllegalArgumentException('cohort input is empty')
    def ids = rows.collect { it.genome_id?.toString() }
    if (ids.any { it == null || !(it ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,63}/) } || ids.unique(false).size() != ids.size())
        throw new IllegalArgumentException('cohort needs distinct, safe genome IDs')
    def paths = []
    def records = rows.sort(false) { it.genome_id }.collect { row ->
        def unknown = row.keySet() - cohortColumns()
        if (unknown) throw new IllegalArgumentException("unknown cohort columns: ${unknown}")
        def record = [:]
        cohortColumns().each { key ->
            if (key == 'out_path') return // Nextflow owns the output directory.
            def value = row[key]
            if (key in ['claims_path', 'genome_layout_path'] && (value == null || value.toString() in ['', '-']))
                throw new IllegalArgumentException("${row.genome_id}: ${key} is required")
            if (value != null && !(value.toString() in ['', '-'])) {
                if (value.toString() =~ /[\t\r\n\x00]/)
                    throw new IllegalArgumentException("${row.genome_id}: invalid ${key}")
                if (key.endsWith('_path')) {
                    // File collections are hashed as unordered bags. Bind identity as a value too:
                    // swapping the same two files between genomes must invalidate the task.
                    record[key] = [index: paths.size(), source: value.toString()]
                    paths << value
                } else record[key] = value.toString()
            }
        }
        record
    }
    tuple(records, paths)
}

def cohortResources(resources) {
    def allowed = ['goOwlFile', 'ec2goFile', 'pathwaysFile', 'thetaFile', 'learnedModelFile',
        'externalContextFile', 'reactionGraphFile', 'reactionEcAliasesFile', 'diffusionMetsFile',
        'cohortMetadataFile', 'cohortOrthogroupsFile', 'cohortRoleEvidenceFile',
        'cohortInteractionEdgesFile', 'cohortReactionDirectionsFile', 'moduleDefsFile', 'moduleFiles',
        'cohortPathwaysFile',
        'taxonomyFile', 'taxonConstraintsFile', 'complexTermsFile', 'essentialFunctionsFile', 'overlayFile', 'reasonerCacheDir']
    if (resources.keySet() - allowed)
        throw new IllegalArgumentException("unknown cohort resources: ${resources.keySet() - allowed}")
    ['reactionGraphFile', 'reactionEcAliasesFile', 'cohortMetadataFile', 'cohortOrthogroupsFile'].each {
        if (!resources[it]) throw new IllegalArgumentException("cohort resource ${it} is required")
    }
    def keys = []
    def paths = []
    resources.keySet().sort().each { key ->
        def values = resources[key] instanceof List ? resources[key] : [resources[key]]
        if (key != 'moduleFiles' && values.size() != 1)
            throw new IllegalArgumentException("${key} requires one file")
        values.each { keys << [key: key, source: it.toString()]; paths << it }
    }
    tuple(keys, paths)
}

process INTEGRATE_COHORT {
    tag 'joint-cohort'
    container 'eclipse-temurin:21-jre'
    cache 'deep'
    errorStrategy 'terminate'
    publishDir params.outdir, mode: 'copy', saveAs: { 'cohort' }

    input:
    tuple val(records), path(sample_files, stageAs: 'samples??/*')
    tuple val(resource_bindings), path(resource_files, stageAs: 'references??/*')
    path gspa_jar, stageAs: 'engine/gspa.jar'
    val settings

    output:
    path 'results', emit: results

    script:
    def files = sample_files instanceof List ? sample_files : [sample_files]
    def refs = resource_files instanceof List ? resource_files : [resource_files]
    def config = new LinkedHashMap(settings)
    resource_bindings.eachWithIndex { binding, index ->
        def key = binding.key
        if (key == 'moduleFiles') {
            if (!config.moduleFiles) config.moduleFiles = []
            config.moduleFiles << refs[index].toString()
        } else if (key == 'reasonerCacheDir') config[key] = 'reasoner-cache'
        else if (key != 'overlayFile') config[key] = refs[index].toString()
    }
    config.cohortInference = true
    config.provenanceFile = 'provenance.json'
    config.suggestionsOut = 'suggestions.tsv'
    if (config.remove('featureExport') == true) config.featureExportDir = 'results/features'
    def header = cohortColumns()
    def lines = records.collect { row ->
        header.collect { key ->
            if (key == 'out_path') return "results/${row.genome_id}/${row.genome_id}_integrated.tsv"
            if (key.endsWith('_path') && row.containsKey(key)) return files[row[key].index as int].toString()
            row.get(key, '-')
        }.join('\t')
    }
    def manifest = ([header.join('\t')] + lines).join('\n')
    def configJson = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(config))
    def overlayIndex = resource_bindings.findIndexOf { it.key == 'overlayFile' }
    def overlay = overlayIndex >= 0 ?
        "--workflow-overlay '" + refs[overlayIndex].toString().replace("'", "'\"'\"'") + "'" : ''
    def cacheIndex = resource_bindings.findIndexOf { it.key == 'reasonerCacheDir' }
    // Never let the engine update a staged symlink to the caller's reference cache.
    def prepareCache = cacheIndex >= 0 ?
        "mkdir -p reasoner-cache; cp -a -- '" + refs[cacheIndex].toString().replace("'", "'\"'\"'") + "/.' reasoner-cache/" : ''
    def dirs = records.collect { "'results/${it.genome_id}'" }.join(' ')
    """
    mkdir -p ${dirs} .tmp
    export TMPDIR="\$PWD/.tmp"
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.io.tmpdir=\$TMPDIR"
    ${prepareCache}
    cat > batch.tsv <<'GSPA_BATCH'
${manifest}
GSPA_BATCH
    cat > workflow.json <<'GSPA_WORKFLOW'
${configJson}
GSPA_WORKFLOW
    java -jar engine/gspa.jar integrate --workflow-config workflow.json --batch-manifest batch.tsv --batch-threads ${task.cpus} ${overlay}
    cp batch.tsv workflow.json results/
    """
}
