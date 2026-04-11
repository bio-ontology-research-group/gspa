package gspa.predictor

import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Registry for discovering and managing predictors.
 * Uses Java ServiceLoader for plugin discovery.
 */
class PredictorRegistry {

    private static final Logger log = LoggerFactory.getLogger(PredictorRegistry)

    private Map<String, Predictor> predictors = [:]

    /**
     * Register a predictor manually.
     */
    void register(Predictor predictor) {
        predictors[predictor.name] = predictor
        log.info("Registered predictor: ${predictor.name}")
    }

    /**
     * Discover predictors via ServiceLoader.
     */
    void discoverPredictors() {
        ServiceLoader.load(Predictor).each { predictor ->
            register(predictor)
        }
    }

    /**
     * Get a predictor by name.
     */
    Predictor get(String name) {
        predictors[name]
    }

    /**
     * Get all registered predictors.
     */
    Collection<Predictor> getAll() {
        predictors.values()
    }

    /**
     * Get all available predictors (tools installed and runnable).
     */
    List<Predictor> getAvailable() {
        predictors.values().findAll { it.isAvailable() }
    }

    /**
     * Get predictors filtered by enabled/disabled lists from config.
     */
    List<Predictor> getEnabled(Set<String> enable, Set<String> disable) {
        def available = getAvailable()
        if (enable && !enable.isEmpty()) {
            available = available.findAll { enable.contains(it.name) }
        }
        if (disable && !disable.isEmpty()) {
            available = available.findAll { !disable.contains(it.name) }
        }
        available
    }

    /**
     * List all registered predictors with their availability status.
     */
    Map<String, Boolean> status() {
        predictors.collectEntries { name, predictor ->
            [(name): predictor.isAvailable()]
        }
    }
}
