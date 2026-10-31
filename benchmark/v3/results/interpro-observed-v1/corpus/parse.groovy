import gspa.predictor.domain.InterProScanPredictor
import groovy.json.JsonOutput
def rows = new InterProScanPredictor().parseOutput(new File(args[0]))
new File(args[1]).withWriter('UTF-8') { out ->
    rows.each { protein, annotations -> annotations.each { a ->
        out.writeLine(JsonOutput.toJson([protein:protein, type:a.type.name(), value:a.value, score:a.score, metadata:a.metadata]))
    } }
}
