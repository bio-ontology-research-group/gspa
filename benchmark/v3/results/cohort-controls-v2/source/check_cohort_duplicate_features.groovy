/* Label-free mechanism replay: production cohort features, bridge and fitted
 * model coefficients. Only feature deltas are attributed; no synthetic vector
 * is claimed to be a supported biological model input or a calibrated outcome.
 * java -cp engine.jar groovy.ui.GroovyMain script fixture.groovy model.json out.json
 */
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import gspa.integration.cohort.*
import gspa.integration.contextbridge.CohortContextBuilder
import gspa.integration.learned.CohortFeatureColumns

Class fixtures = new GroovyClassLoader(getClass().classLoader).parseClass(new File(args[0]))
Map model = new JsonSlurper().parse(new File(args[1])) as Map
List cases = []
for (int copies : [0, 1, 3, 5]) {
    def input = fixtures.newInput()
    fixtures.addGenome(input, [id: 'G0', clade: 'ct', cluster: 'st', status: 'COMPLETE', closed: true,
                              family: 'OG1', candidate: 'T'])
    fixtures.addGenome(input, [id: 'P0', clade: 'c0', cluster: 's0', status: 'COMPLETE', closed: true,
                              family: 'OG1', anchors: ['N1', 'N2']])
    fixtures.addGenome(input, [id: 'P1', clade: 'c1', cluster: 's1', status: 'COMPLETE', closed: true,
                              family: 'OG1', anchors: []])
    copies.times { int k ->
        fixtures.addGenome(input, [id: "D${k}", clade: 'c0', cluster: 's0', status: 'COMPLETE', closed: true,
                                  family: 'OG1', anchors: ['N1', 'N2']])
    }
    def config = new CohortConfig(permutations: 200, seed: 20261026, holdout: HoldoutLevel.CLADE)
    def record = new CohortInference(input, fixtures.graph(), config).run().record('G0', 'G0_20', 'T')
    Map<String, Double> columns = [:]
    CohortFeatureColumns.extract(CohortContextBuilder.weighted(record), { String n, Object v -> columns[n] = v as double })
    cases << [copies: copies, features: record.features, model_columns: columns]
}
Map base = cases.first().model_columns
cases.drop(1).each { Map c ->
    Map changes = [:]
    (base.keySet() + c.model_columns.keySet()).sort().each { String name ->
        double a = base.get(name) ?: 0d, b = c.model_columns.get(name) ?: 0d
        if (Math.abs(a - b) > 1e-12d) {
            int i = model.contract.features.indexOf(name)
            assert i >= 0
            double delta = (b - a) * (model.logistic.weights[i] as double) / (model.normalization.scale[i] as double)
            changes[name] = [before: a, after: b, model_margin_delta: delta]
        }
    }
    c.changed_columns = changes
    c.total_margin_delta_from_changed_columns = changes.values().sum { it.model_margin_delta } ?: 0d
}
new File(args[2]).text = JsonOutput.prettyPrint(JsonOutput.toJson([
    interpretation: 'synthetic duplicated-evidence mechanism; not calibrated biological predictions',
    cases: cases])) + '\n'
