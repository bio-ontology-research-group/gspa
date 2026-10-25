package gspa.config

import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory
import com.fasterxml.jackson.databind.DeserializationFeature
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Loads and merges GSPA configuration from multiple YAML sources.
 * Merge order: built-in defaults -> kingdom preset -> user config file -> CLI overrides.
 */
class ConfigLoader {

    private static final Logger log = LoggerFactory.getLogger(ConfigLoader)

    private static final ObjectMapper mapper = new ObjectMapper(new YAMLFactory())
        .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, false)

    /**
     * Load configuration from a YAML file.
     */
    static GspaConfig loadFromFile(File configFile) {
        log.info("Loading configuration from: ${configFile}")
        mapper.readValue(configFile, GspaConfig)
    }

    /**
     * Load configuration from a YAML string.
     */
    static GspaConfig loadFromString(String yaml) {
        mapper.readValue(yaml, GspaConfig)
    }

    /**
     * Load a kingdom preset from built-in resources.
     */
    static GspaConfig loadPreset(String presetName) {
        String resourcePath = "config/${presetName}.yaml"
        def stream = ConfigLoader.classLoader.getResourceAsStream(resourcePath)
        if (stream == null) {
            log.warn("No built-in preset found for: ${presetName}, using defaults")
            return new GspaConfig()
        }
        log.info("Loading preset: ${presetName}")
        stream.withCloseable { mapper.readValue(it, GspaConfig) }
    }

    /**
     * Build the final configuration by merging sources.
     *
     * @param userConfigFile Optional user config file path
     * @param kingdom Optional kingdom override
     * @param cliOverrides Map of CLI flag overrides
     */
    static GspaConfig buildConfig(File userConfigFile = null, String kingdom = null,
                                   Map<String, Object> cliOverrides = [:], String workflowPreset = null) {
        // Merge sparse YAML maps BEFORE constructing typed configs. Deserializing
        // a partial YAML into GspaConfig first fills absent keys with class defaults
        // and would silently erase settings inherited from the preset.
        Map values = mapper.convertValue(new GspaConfig(), Map)
        if (kingdom) {
            values.input.kingdom = kingdom
            def stream = ConfigLoader.classLoader.getResourceAsStream("config/${kingdom}.yaml")
            if (stream != null) {
                Map preset = stream.withCloseable { mapper.readValue(it, Map) }
                values = deepMerge(values, preset, true)
            }
        }
        if (workflowPreset) {
            values = deepMerge(values, WorkflowPresets.annotationDefaults(workflowPreset, kingdom), true)
        }
        if (userConfigFile != null) {
            if (!userConfigFile.isFile()) {
                throw new IllegalArgumentException("Configuration file does not exist: ${userConfigFile}")
            }
            values = deepMerge(values, mapper.readValue(userConfigFile, Map) ?: [:], true)
        }
        GspaConfig config = mapper.convertValue(values, GspaConfig)
        applyOverrides(config, cliOverrides)

        config
    }

    /**
     * Merge two configs, with 'override' taking precedence for non-null/non-default values.
     */
    static GspaConfig merge(GspaConfig base, GspaConfig override) {
        // Simple merge: serialize both to maps, deep-merge, deserialize back
        def baseMap = mapper.convertValue(base, Map)
        def overrideMap = mapper.convertValue(override, Map)
        def merged = deepMerge(baseMap, overrideMap)
        mapper.convertValue(merged, GspaConfig)
    }

    private static Map deepMerge(Map base, Map override, boolean replaceNulls = false) {
        Map result = new LinkedHashMap(base)
        override.each { key, value ->
            if (value instanceof Map && result[key] instanceof Map) {
                result[key] = deepMerge(result[key] as Map, value as Map, replaceNulls)
            } else if (value != null || replaceNulls) {
                result[key] = value
            }
        }
        result
    }

    private static void applyOverrides(GspaConfig config, Map<String, Object> overrides) {
        overrides.each { key, value ->
            switch (key) {
                case 'kingdom':
                    config.input.kingdom = value as String
                    break
                case 'input-type':
                    config.input.type = value as String
                    break
                case 'output-dir':
                    config.output.outputDir = value as String
                    break
                case 'go-owl':
                    config.quality.goOwlFile = value as String
                    break
                case 'database':
                    config.predictors.similarity.database = value as String
                    config.predictors.similarity.enabled = true
                    break
                case 'evalue':
                    config.predictors.similarity.evalue = value as double
                    break
                case 'no-coherence':
                    config.quality.coherence.process = false
                    config.quality.coherence.pathway = false
                    config.quality.coherence.complex = false
                    break
                case 'no-consistency':
                    config.quality.consistency.taxonConstraints = false
                    break
            }
        }
    }

    /**
     * Write a config to YAML file.
     */
    static void writeConfig(GspaConfig config, File outputFile) {
        mapper.writeValue(outputFile, config)
    }

    static String toYaml(GspaConfig config) {
        mapper.writeValueAsString(config)
    }
}
