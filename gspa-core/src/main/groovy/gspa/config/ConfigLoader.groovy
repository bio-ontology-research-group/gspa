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
        String resourcePath = "/config/${presetName}.yaml"
        def stream = ConfigLoader.classLoader.getResourceAsStream(resourcePath)
        if (stream == null) {
            log.warn("No built-in preset found for: ${presetName}, using defaults")
            return new GspaConfig()
        }
        log.info("Loading preset: ${presetName}")
        mapper.readValue(stream, GspaConfig)
    }

    /**
     * Build the final configuration by merging sources.
     *
     * @param userConfigFile Optional user config file path
     * @param kingdom Optional kingdom override
     * @param cliOverrides Map of CLI flag overrides
     */
    static GspaConfig buildConfig(File userConfigFile = null, String kingdom = null,
                                   Map<String, Object> cliOverrides = [:]) {
        // Start with defaults
        GspaConfig config = new GspaConfig()

        // Apply kingdom preset if specified
        if (kingdom) {
            config.input.kingdom = kingdom
            def preset = loadPreset(kingdom)
            config = merge(config, preset)
        }

        // Apply user config file
        if (userConfigFile?.exists()) {
            def userConfig = loadFromFile(userConfigFile)
            config = merge(config, userConfig)
        }

        // Apply CLI overrides
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

    private static Map deepMerge(Map base, Map override) {
        Map result = new LinkedHashMap(base)
        override.each { key, value ->
            if (value instanceof Map && result[key] instanceof Map) {
                result[key] = deepMerge(result[key] as Map, value as Map)
            } else if (value != null) {
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
}
