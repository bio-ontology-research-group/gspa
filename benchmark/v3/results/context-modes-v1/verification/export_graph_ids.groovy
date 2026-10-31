import gspa.ontology.ReactionGraphLoader
import groovy.json.JsonOutput
def graph = ReactionGraphLoader.load(new File(args[0]), null, 100d)
new File(args[1]).text = JsonOutput.toJson(graph.reactions.keySet().sort()) + '\n'
println "production graph reactions: ${graph.reactions.size()}"
