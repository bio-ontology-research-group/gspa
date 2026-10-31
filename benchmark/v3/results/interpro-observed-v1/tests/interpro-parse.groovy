import gspa.predictor.domain.InterProScanPredictor
import groovy.json.JsonOutput
File fixture = new File(args[0])
def result = new InterProScanPredictor().parseOutput(fixture)
def rows = result.collectMany { protein, anns -> anns.collect { [protein:protein, type:it.type.name(), value:it.value, score:it.score, metadata:it.metadata] } }
new File(args[1]).text = JsonOutput.prettyPrint(JsonOutput.toJson(rows)) + '\n'
println(JsonOutput.toJson([rows:rows.size(), negativeScores:rows.count { it.score < 0 }]))
