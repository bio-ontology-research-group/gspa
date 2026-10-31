package gspa.predictor.pathway

import gspa.config.GspaConfig
import gspa.model.*
import gspa.predictor.*
import spock.lang.Specification
import spock.lang.TempDir

class PathwayInputContractSpec extends Specification {
    @TempDir File temp

    static class Capture extends GapseqPredictor {
        String input
        List<String> command
        @Override protected ProcessResult execute(List<String> cmd, File dir) {
            command = cmd
            input = new File(cmd.last()).text
            new File(dir, 'proteins-all-Reactions.tbl').text =
                'pathway\trxn\tname\tec\tkeyrea\tfile\tdbhit\tspont\ttype\tsrc\tis_complex\tsubunit_count\tsubunits\tqseqid\tpident\tevalue\tbitscore\tqcovs\tstitle\tsstart\tsend\tcomplex\texception\tstatus\tsubunits_found\tsubunit_undefined_found\tcomplex.status\tpathway.status\n'
            new ProcessResult(exitCode: 0, stdout: '', stderr: '')
        }
    }

    private Genome genome() {
        new Genome(id:'g', contigs:[new Contig(id:'c', sequence:'ATGAAGTAA', proteins:[
            new Protein(id:'p1',sequence:'MK'),new Protein(id:'p2',sequence:'MA')])])
    }

    def 'default gapsmith search receives proteins even when nucleotide contigs exist'() {
        given:
        def p = new Capture()
        when:
        p.predictGenome(genome())
        then:
        p.input == '>p1\nMK\n>p2\nMA\n'
        p.lastRunStatus.outcome == PredictorOutcome.SUCCESS_ZERO
    }

    def 'gapsmith rejects nucleotide mode before launching a protein aligner'() {
        given:
        def p = new Capture(target:'genome')
        when:
        p.predictGenome(genome())
        then:
        thrown(IllegalArgumentException)
        p.command == null
        p.lastRunStatus.outcome == PredictorOutcome.FAILED
    }

    def 'missing protein sequences are partial coverage'() {
        given:
        def g = genome()
        g.proteins[1].sequence = null
        def p = new Capture(target:'proteome')
        when:
        p.predictGenome(g)
        then:
        p.input == '>p1\nMK\n'
        p.lastRunStatus.outcome == PredictorOutcome.PARTIAL
        p.lastRunStatus.skippedInputs == ['p2']
    }

    def 'gapseq taxonomy uses -t and protein mode is explicit'() {
        given:
        def p = new GapseqPredictor(tool:'gapseq',target:'proteome',taxonomy:'Archaea')
        when:
        def command = p.buildCommand(new File(temp,'input.faa'),temp)
        then:
        command[command.indexOf('-t')+1] == 'Archaea'
        command[command.indexOf('-M')+1] == 'prot'
        !command.contains('-b')
    }

    def 'both pathway tools receive the configured thread budget'() {
        given:
        def p = new GapseqPredictor(tool:tool,target:'proteome',threads:3)
        when:
        def command = p.buildCommand(new File(temp,'input.faa'),temp)
        then:
        command[command.indexOf('-K')+1] == '3'
        where:
        tool << ['gapsmith','gapseq']
    }

    def 'pipeline honors explicit pathway tool configuration'() {
        given:
        def config = new GspaConfig()
        config.predictors.pathway.tool = 'gapseq'
        def pipeline = new AnnotationPipeline(config)
        when:
        pipeline.configure()
        then:
        pipeline.genomePredictors.find { it instanceof GapseqPredictor }.tool == 'gapseq'
    }

    def 'pipeline forwards pathway references target executable taxonomy and threads'() {
        given:
        def config = new GspaConfig()
        config.predictors.pathway.with {
            tool='gapsmith'; target='proteome'; dataDir='/reference/gapsmith'
            executablePath='/tools/gapsmith'; taxonomy='Archaea'; threads=3
        }
        def pipeline = new AnnotationPipeline(config)
        when:
        pipeline.configure()
        def p = pipeline.genomePredictors.find { it instanceof GapseqPredictor }
        then:
        p.tool == 'gapsmith'
        p.target == 'proteome'
        p.dataDir == '/reference/gapsmith'
        p.executablePath == '/tools/gapsmith'
        p.taxonomy == 'Archaea'
        p.threads == 3
    }

    def 'unknown pathway tool is rejected rather than silently invoking gapsmith'() {
        when:
        new GapseqPredictor(tool:'typo').buildCommand(new File(temp,'input.faa'),temp)
        then:
        thrown(IllegalArgumentException)
    }
}
