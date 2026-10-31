package gspa.predictor.pathway

import groovy.json.JsonSlurper
import gspa.config.GspaConfig
import gspa.integration.EvidenceType
import gspa.model.*
import gspa.predictor.*
import spock.lang.Specification
import spock.lang.TempDir

class PathwayObservationSpec extends Specification {
    @TempDir File temp
    static class Replay extends GapseqPredictor {
        String output
        String pathways
        boolean available = true
        int code = 0
        @Override boolean isAvailable() { available }
        @Override protected ProcessResult execute(List<String> command, File dir) {
            if (code == 0) new File(dir,'proteins-all-Reactions.tbl').text=output
            if (code == 0 && pathways != null) new File(dir,'proteins-all-Pathways.tbl').text=pathways
            new ProcessResult(exitCode:code,stdout:'',stderr:code ? 'injected tool failure' : '')
        }
    }
    String table() { getClass().getResource('/observed-output/gapsmith/training-subset-Reactions.tbl').text }

    def 'completed pathway observations automatically reach shared integration as candidate gaps'() {
        given:
        File fasta = new File(temp, 'proteins.faa')
        fasta.text = ['ABB58505.1','ABB58544.1','ABB56978.1','ABB58152.1','ABB57067.1'].collect {
            ">${it}\n" + 'MA'*350 + '\n'
        }.join('')
        def cfg = new GspaConfig()
        cfg.input.proteinsOnly = true; cfg.input.proteinFasta = fasta.path
        cfg.predictors.similarity.enabled = false; cfg.predictors.domains.enabled = false
        cfg.predictors.pathway.enabled = false; cfg.predictors.operons.enabled = false
        cfg.integration.enabled = true
        cfg.integration.workflow.lite = true
        cfg.integration.workflow.enabledPriors = 'none'
        def pipeline = new AnnotationPipeline(cfg)
        pipeline.genomePredictors << new Replay(output:table(),
            pathways:'ID\tPrediction\n2FLHYD\tTRUE\nARGSYN-SUCC-PWY\tTRUE\n')

        when:
        def genome = pipeline.annotate(fasta)

        then:
        pipeline.lastIntegration.state.metabolicGaps.size() == 2
        pipeline.lastIntegration.state.metabolicGaps*.reactionId.sort() == ['ALPHA-FUCOSIDASE', 'RXN-9778']
        pipeline.lastIntegration.state.metabolicGapsAvailable
        pipeline.pathwaySearchReports.gapsmith.pathway_rows.size() == 2

        when: 'standalone integration of the same genome does not reuse prior observations'
        def fresh = pipeline.integrate(genome, pipeline.rawClaims)

        then:
        fresh.state.metabolicGaps.isEmpty()
        !fresh.state.metabolicGapsAvailable
    }

    def 'pipeline retains raw pathway observations, sequence claims and clears them between genomes'() {
        given:
        File fasta=new File(temp,'proteins.faa')
        fasta.text=['ABB58505.1','ABB58544.1','ABB56978.1','ABB58152.1','ABB57067.1'].collect {
            ">${it}\n"+'MA'*350+'\n'
        }.join('')
        def cfg=new GspaConfig()
        cfg.input.proteinsOnly=true; cfg.input.proteinFasta=fasta.path
        cfg.predictors.similarity.enabled=false;cfg.predictors.domains.enabled=false
        cfg.predictors.pathway.enabled=false;cfg.predictors.operons.enabled=false
        def pipeline=new AnnotationPipeline(cfg)
        def predictor=new Replay(output:table())
        pipeline.genomePredictors << predictor
        when:
        def genome=pipeline.annotate(fasta)
        pipeline.writeOutput(genome,new File(temp,'first'))
        def report=new JsonSlurper().parse(new File(temp,'first/proteins_gapsmith_observations.json'))
        then:
        report.reaction_rows.size()==8
        report.outcome=='SUCCESS'
        pipeline.rawClaims.size()==5
        pipeline.rawClaims.every { it.evidenceType==EvidenceType.SEQUENCE_SIMILARITY }
        when:
        predictor.available=false
        def later=pipeline.annotate(fasta)
        pipeline.writeOutput(later,new File(temp,'first'))
        then:
        pipeline.pathwaySearchReports.isEmpty()
        pipeline.rawClaims.isEmpty()
        !new File(temp,'first/proteins_gapsmith_observations.json').exists()
    }

    def 'unknown output protein fails the actual run contract rather than being silently dropped'() {
        given:
        def predictor=new Replay(output:table())
        def genome=new Genome(id:'g',contigs:[new Contig(id:'c',proteins:[new Protein(id:'different',sequence:'MA')])])
        when:
        predictor.predictGenome(genome)
        then:
        thrown(ToolOutputContractException)
        predictor.lastRunStatus.outcome==PredictorOutcome.FAILED
        predictor.reactionObservations.size()==8
    }

    def 'nucleotide hits cannot be asserted as protein annotations without coordinate attribution'() {
        given:
        def predictor=new Replay(tool:'gapseq',target:'genome',output:table())
        def genome=new Genome(id:'g',contigs:[new Contig(id:'c',sequence:'ATGAAA')])
        when:
        predictor.predictGenome(genome)
        then:
        thrown(ToolOutputContractException)
        predictor.lastRunStatus.outcome==PredictorOutcome.FAILED
    }

    def 'protein batch rejects nucleotide mode before process execution'() {
        given:
        def predictor=new Replay(tool:'gapseq',target:'genome',output:table())
        when:
        predictor.predictBatch([new Protein(id:'p',sequence:'MA')])
        then:
        thrown(IllegalArgumentException)
        predictor.lastRunStatus.outcome==PredictorOutcome.FAILED
        predictor.reactionObservations.isEmpty()
    }

    def 'failed repeated protein batch clears the previous observations'() {
        given:
        def predictor=new Replay(output:table(), pathways:'ID\tPrediction\n2FLHYD\tTRUE\n')
        def proteins=['ABB58505.1','ABB58544.1','ABB56978.1','ABB58152.1','ABB57067.1'].collect { new Protein(id:it,sequence:'MA') }
        predictor.predictBatch(proteins)
        when:
        predictor.code=2
        predictor.predictBatch(proteins)
        then:
        predictor.lastRunStatus.outcome==PredictorOutcome.FAILED
        predictor.reactionObservations.isEmpty()
        predictor.reactionTableSha256==null
        predictor.pathwayObservations==null
        predictor.pathwayTableSha256==null
    }

    def 'malformed pathway tables fail the adapter instead of certifying gap observations'() {
        given:
        def predictor = new Replay(output:table(), pathways:pathways)
        def proteins = ['ABB58505.1','ABB58544.1','ABB56978.1','ABB58152.1','ABB57067.1'].collect {
            new Protein(id:it, sequence:'MA')
        }

        when:
        predictor.predictBatch(proteins)

        then:
        thrown(ToolOutputContractException)
        predictor.lastRunStatus.outcome == PredictorOutcome.FAILED
        predictor.reactionTableSha256 == null
        predictor.pathwayTableSha256 == null

        where:
        pathways << ['', 'ID\tUnknown\nP\tTRUE\n', 'ID\tPrediction\nP\tmaybe\n',
                      'ID\tPrediction\nP\tTRUE\nP\tFALSE\n', 'ID\tPrediction\nP\tTRUE\textra\n']
    }
}
