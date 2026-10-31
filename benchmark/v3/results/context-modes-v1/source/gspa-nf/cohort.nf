#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { INTEGRATE_COHORT; cohortBatch; cohortResources } from './modules/cohort'

/* Integration-only entry point; the main annotation pipeline uses the same process. */
workflow {
    if (!params.cohort_manifest || !params.cohort_resources || !params.gspa_jar)
        error 'requires --cohort_manifest TSV, --cohort_resources JSON and --gspa_jar'
    def manifest = file(params.cohort_manifest, checkIfExists: true)
    def resolvePath = { value, parent ->
        def p = java.nio.file.Paths.get(value.toString())
        file(p.isAbsolute() ? p : parent.resolve(p), checkIfExists: true)
    }
    def resourcesFile = file(params.cohort_resources, checkIfExists: true)
    def resourceSpec = new groovy.json.JsonSlurper().parse(resourcesFile.toFile())
    def resources = resourceSpec.collectEntries { key, value ->
        [key, value instanceof List ? value.collect { resolvePath.call(it, resourcesFile.parent) } : resolvePath.call(value, resourcesFile.parent)]
    }
    def settings = [lite: true, essentialProfile: params.essential_profile ?: 'bacteria',
                    enabledPriors: params.enable_priors ?: 'essentiality,coherence,gap_filling,genomic_context',
                    cohortConfig: [holdout: params.cohort_holdout ?: 'CLADE', contextMode: params.cohort_context_mode ?: 'UNIVERSAL'],
                    featureExport: params.feature_export ?: false]
    def batch = Channel.fromPath(manifest).splitCsv(header: true, sep: '\t').map { row ->
        row.collectEntries { key, value ->
            [key, key.endsWith('_path') && key != 'out_path' && value && value != '-' ? resolvePath.call(value, manifest.parent) : value]
        }
    }.toList().map { cohortBatch(it) }
    INTEGRATE_COHORT(batch, cohortResources(resources), file(params.gspa_jar, checkIfExists: true), settings)
}
