import com.fasterxml.jackson.databind.ObjectMapper
import gspa.config.GspaConfig
import gspa.predictor.AnnotationPipeline
import gspa.predictor.genecalling.PyrodigalCaller
import gspa.cli.AntiSmashClaimsCommand

File input = new File(args[0])
String executable = args[1]
boolean meta = args[2] == 'meta'
File out = new File(args[3]); out.mkdirs()
def config = new GspaConfig()
config.integration.enabled = false
config.input.type = meta ? 'metagenome' : 'genome'
def pipeline = new AnnotationPipeline(config)
pipeline.parallel = false
pipeline.geneCaller = new PyrodigalCaller(metaMode:meta, executablePath:executable)
def genome = pipeline.annotate(input)
new File(out, 'genes.tsv').withWriter('UTF-8') { w ->
    w.writeLine('protein_id\tcontig_id\tstart\tend\tstrand\tsequence')
    genome.proteins.each { p -> w.writeLine([p.id,p.contig.id,p.start,p.end,p.sourceFeature.strand.symbol,p.sequence].join('\t')) }
}
new File(out, 'result.json').text = new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString([
    protein_count: genome.proteinCount,
    feature_count: genome.contigs.sum { it.features.size() },
    statuses: pipeline.runStatuses.collect { it.toMap() }
]) + '\n'
if (args.length > 4) {
    def joined = AntiSmashClaimsCommand.loadCoordinates(input, new File(args[4]), new File(args[5]))
    new File(out, 'strict-genes.tsv').withWriter('UTF-8') { w ->
        w.writeLine('protein_id\tcontig_id\tstart\tend\tstrand\tsequence')
        joined.proteins.each { p -> w.writeLine([p.id,p.contig.id,p.start,p.end,p.sourceFeature.strand.symbol,p.sequence].join('\t')) }
    }
}
