import groovy.json.JsonOutput
import gspa.predictor.pathway.GapseqPredictor

def predictor = new GapseqPredictor(target: 'proteome')
def annotations = predictor.parseOutput(new File(args[0]))
def rows = annotations.collectMany { protein, values -> values.collect { a ->
    [protein_id:protein,type:a.type.name(),value:a.value,source:a.source,score:a.score,
     evidence_type:a.evidenceType?.name(),metadata:a.metadata]
} }.sort { a,b -> (a.protein_id+'|'+a.type+'|'+a.value) <=> (b.protein_id+'|'+b.type+'|'+b.value) }
new File(args[1]).text = JsonOutput.prettyPrint(JsonOutput.toJson([proteins:annotations.size(),annotations:rows]))+'\n'
println JsonOutput.toJson([proteins:annotations.size(),annotations:rows.size()])
