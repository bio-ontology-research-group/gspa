import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput
import gspa.integration.IntegrationRunAudit
import java.security.DigestOutputStream
import java.security.MessageDigest

def method = IntegrationRunAudit.getDeclaredMethod('changes', Map, Map)
method.accessible = true
method.invoke(null, [:], [:]) // load the actual implementation before the measurement
int count = args[0] as int
Map<String,Double> after = new LinkedHashMap<>()
for (int i=0; i<count; i++) after.put('protein.'+i+'|EC|1.1.1.1', (i % 1001)/1000d)
def runtime=Runtime.runtime
3.times { System.gc(); Thread.sleep(100) }
long baseline=runtime.totalMemory()-runtime.freeMemory()
long started=System.nanoTime()
def changes=method.invoke(null, [:], after)
long captureNanos=System.nanoTime()-started
3.times { System.gc(); Thread.sleep(100) }
long retained=runtime.totalMemory()-runtime.freeMemory()
def digest=MessageDigest.getInstance('SHA-256')
new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(
    new DigestOutputStream(OutputStream.nullOutputStream(),digest), changes)
println JsonOutput.toJson([count:changes.size(),json_sha256:digest.digest().encodeHex().toString(),
 baseline_heap_bytes:baseline,retained_heap_bytes:retained,delta_heap_bytes:retained-baseline,
 capture_seconds:captureNanos/1e9d,
 interpretation:'Controlled actual changes() snapshot; heap used after explicit GC. Not whole-application peak RSS, runtime or a biological result.'])
