// P10 driver: runs the unmodified `gspa integrate` command of the jar on the
// classpath, with exactly one change: the AnomalyExtractor's global cap of
// 150 anomalies per genome is lifted to $P10_MAX_ANOMALIES (default 100000).
//
// Why: AnomalyExtractor emits anomalies in a fixed order (PATHWAY_GAP, then
// INCOMPLETE_COMPLEX, BORDERLINE_ESSENTIAL, BORDERLINE_CONFLICT) and stops at
// 150. In large bacterial proteomes the pathway gaps alone exhaust the cap, so
// the later categories are never emitted. Lifting the cap changes neither the
// order nor the content of any record: the capped (production) output is the
// first 150 lines of the uncapped output, which trigger_rates_v2.py relies on.
//
// Usage:
//   P10_MAX_ANOMALIES=100000 java -Xmx12g -cp gspa-unified.jar groovy.ui.GroovyMain \
//       p10_integrate.groovy --claims ... --out ... --anomalies-out ... [integrate flags]
import gspa.agent.AnomalyExtractor

int cap = Integer.parseInt(System.getenv('P10_MAX_ANOMALIES') ?: '100000')
def ctor = AnomalyExtractor.getConstructor()
AnomalyExtractor.metaClass.constructor = { ->
    def e = ctor.newInstance()
    e.maxAnomalies = cap
    e
}
System.err.println("[p10_integrate] AnomalyExtractor.maxAnomalies=${cap}")
int rc = new picocli.CommandLine(new gspa.cli.IntegrateCommand()).execute(args as String[])
System.exit(rc)
