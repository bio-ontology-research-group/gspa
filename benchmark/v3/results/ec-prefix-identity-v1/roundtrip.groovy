import groovy.json.JsonOutput
import gspa.integration.ClaimExtractor
import gspa.integration.IntegrationState
import gspa.integration.learned.ContextFeatureExtractor
import gspa.model.Genome

def x=new ClaimExtractor()
def before=x.readClaimsJsonl(new File(args[0]))
def out=new File(args[1])
ClaimExtractor.writeClaimsJsonl(before,out)
def after=x.readClaimsJsonl(out)
def f=new ContextFeatureExtractor()
def state=new IntegrationState(new Genome(id:'g'))
def a=f.extractAll(before,state)[0]
def b=f.extractAll(after,state)[0]
println JsonOutput.prettyPrint(JsonOutput.toJson([before_contracts:a.sequenceContracts,after_contracts:b.sequenceContracts,
    before_metadata_types:before.collect{it.metadata.predictor_contract.getClass().name},
    after_metadata_types:after.collect{it.metadata.predictor_contract.getClass().name},
    same_sources:a.sequenceSources==b.sequenceSources,same_feature_bits:a.digest()==b.digest()]))
