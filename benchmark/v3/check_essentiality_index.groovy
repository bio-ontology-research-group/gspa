/* Run with java -cp <production.jar> groovy.ui.GroovyMain this-file output.json.
 * Synthetic semantic/performance replay; not an annotation-accuracy benchmark.
 */
import groovy.json.JsonOutput
import gspa.config.EssentialFunctions
import gspa.integration.IntegrationState
import gspa.integration.prior.EssentialityPrior
import gspa.model.Genome
import gspa.ontology.GoOntology

class ReplayOntology extends GoOntology {
    Map<String, Set<String>> descendants = [
        'GO:E1': ['GO:D2', 'GO:D1'] as Set,
        'GO:E2': ['GO:D1', 'GO:D3'] as Set,
        'GO:E3': ['GO:D4'] as Set]
    @Override Set<String> getDescendants(String id) { descendants.get(id) ?: ([] as Set) }
    @Override Set<String> getAncestors(String id) {
        descendants.findAll { k, v -> v.contains(id) }.keySet() as Set
    }
    @Override Set<String> propagateAnnotations(Set<String> ids) {
        Set<String> result = new LinkedHashSet<>(ids)
        ids.each { result.addAll(getAncestors(it)) }
        result
    }
    @Override String getAspect(String id) { id == 'GO:D4' ? 'MF' : 'BP' }
}

def clean
clean = { Object value ->
    if (value instanceof Double && !Double.isFinite(value)) return value.toString()
    if (value instanceof Map) return value.collectEntries { k, v -> [k, clean(v)] }
    if (value instanceof Collection) return value.collect { clean(it) }
    value
}
def logit = { double p -> Math.log(p / (1.0d - p)) }
def ontology = new ReplayOntology()
List<Map> cases = []
for (String mode : ['legacy', 'targeted']) {
    for (int seed = 0; seed < 40; seed++) {
        def random = new Random(seed)
        def state = new IntegrationState(new Genome(id: 'synthetic'))
        state.essentialFunctions = new EssentialFunctions(functions: [
            'GO:E1': 'test', 'GO:E2': 'test', 'GO:E3': 'test'])
        state.goOntology = seed % 3 == 0 ? null : ontology
        List<String> keys = []
        List<String> terms = ['GO:E1', 'GO:E2', 'GO:E3', 'GO:D1', 'GO:D2', 'GO:D3', 'GO:D4', 'GO:other']
        for (int i = 0; i < 40; i++) keys << "sp|protein${i % 10}|GO|${terms[i % terms.size()]}".toString()
        for (int i = 0; i < 25; i++) keys << "ec${i}|EC|1.1.1.${i}".toString()
        keys << 'malformed'
        Collections.shuffle(keys, random)
        Map<String, Double> values = new LinkedHashMap<>()
        keys.each { values[it] = logit([0.10d, 0.16d, 0.3d, 0.4d, 0.45d, 0.5d, 0.8d][random.nextInt(7)]) }
        state.posteriorLogOdds = new LinkedHashMap<>(values)
        state.likelihoodLogOdds = seed % 4 == 0 ? [:] : values
        state.claimCounts = seed % 7 == 0 ? [:] : values.collectEntries { k, v -> [k, 1 + random.nextInt(3)] }
        state.operons = [['sp|protein1', 'sp|protein2', 'sp|protein3'], ['sp|protein4', 'sp|protein5']]
        state.pinnedFloors['sp|protein1|GO|GO:D1'] = logit(0.9d)
        def prior = new EssentialityPrior(rescueMode: mode, maxRescuePerCategory: 1 + seed % 2)
        for (int iteration = 0; iteration < 3; iteration++) {
            if (iteration == 1) {
                // Same-size key replacement and changed likelihoods must not reuse an old candidate index.
                state.posteriorLogOdds.remove(keys[0])
                state.posteriorLogOdds['new|GO|GO:E3'] = logit(0.45d)
                state.likelihoodLogOdds = new LinkedHashMap<>(state.posteriorLogOdds)
            }
            if (iteration == 2) {
                state.posteriorLogOdds['new|GO|GO:E3'] = logit(0.8d)
                state.likelihoodLogOdds = new LinkedHashMap<>(state.posteriorLogOdds)
            }
            prior.beginIteration(state)
            cases << clean([mode: mode, seed: seed, iteration: iteration,
                boosts: state.posteriorLogOdds.keySet().collectEntries { k -> [k, prior.logOddsBoost('', k, state)] },
                floors: new LinkedHashMap<>(state.pinnedFloors), rescues: state.rescueLog.collect { new LinkedHashMap<>(it) }])
        }
    }
}
List<Map> timings = []
for (int round = 0; round < 4; round++) {
    def state = new IntegrationState(new Genome(id: 'synthetic-ec'))
    state.essentialFunctions = new EssentialFunctions(functions: (0..<32).collectEntries { ["GO:E${it}".toString(), 'test'] })
    (0..<100000).each { int i -> state.posteriorLogOdds["p${i}|EC|1.1.1.1".toString()] = -0.5d }
    state.claimCounts = ['p0|EC|1.1.1.1': 1]
    def prior = new EssentialityPrior()
    long start = System.nanoTime()
    prior.beginIteration(state)
    double seconds = (System.nanoTime() - start) / 1e9d
    assert state.rescueLog.isEmpty() && state.pinnedFloors.isEmpty()
    assert prior.logOddsBoost('p0', 'p0|EC|1.1.1.1', state) == 0.0d
    timings << [round: round, candidates: 100000, essential_categories: 32, elapsed_seconds: seconds]
}
new File(args[0]).text = JsonOutput.prettyPrint(JsonOutput.toJson([
    schema: 'gspa.essentiality-index-replay/1', scope: 'Synthetic semantic and cost check only.',
    cases: cases, timings: timings])) + '\n'
