package gspa.predictor.genecalling

import gspa.io.FastaReader
import gspa.io.GeneModelReader
import gspa.model.*
import gspa.predictor.*
import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Timeout

class GeneCallerContractSpec extends Specification {
    @TempDir File dir
    String prefix = '##gff-version 3\n# Sequence Data: seqnum=1;seqlen=120;seqhdr="chr"\n'
    String cds = 'chr\tfixture\tCDS\t1\t30\t.\t+\t0\tID=chr_1\n' +
        'chr\tfixture\tCDS\t61\t90\t.\t-\t0\tID=chr_2\n'
    String faa = '>chr_1 # 1 # 30 # 1 # ID=1_1;partial=00\nMMMMMMMMMM\n' +
        '>chr_2 # 61 # 90 # -1 # ID=1_2;partial=00\nAAAAAAAAAA\n'

    File file(String name, String text) {
        File f = new File(dir, name); f.text = text; f
    }
    File fixture(String name) {
        File f = new File(dir, name)
        f.bytes = getClass().getResourceAsStream('/genecalling/pyrodigal-3.7.1/' + name).bytes
        f
    }
    Genome genome() {
        def g = new Genome(id: 'sample')
        g.addContig(new Contig(id: 'chr', sequence: 'ATG' * 40))
        g
    }
    File executable(String body) {
        File f = file('tool-' + UUID.randomUUID(), '#!' + (System.getenv('GSPA_TEST_PYTHON') ?: '/usr/bin/python3') + '\n' + body)
        f.setExecutable(true)
        f
    }
    File replay(String gff, String proteins, int exit = 0) {
        executable("""import sys
from pathlib import Path
a=sys.argv
if '-V' in a or '-v' in a:
 print('pyrodigal fixture');sys.exit(0)
Path(a[a.index('-o')+1]).write_bytes(bytes.fromhex('${gff.bytes.encodeHex()}'))
Path(a[a.index('-a')+1]).write_bytes(bytes.fromhex('${proteins.bytes.encodeHex()}'))
print('native diagnostic',file=sys.stderr)
sys.exit(${exit})
""")
    }
    PredictorRunStatus status(GeneCaller caller) {
        caller instanceof StatusReportingPredictor ? caller.lastRunStatus : null
    }

    def 'both callers join reversed proteins by identity and retain public coordinate-order IDs'() {
        given:
        def g = genome()
        def reversed = faa.split('(?=>chr_)').reverse().join('')
        def caller = type.newInstance(executablePath: replay(prefix + cds, reversed).path)
        when:
        caller.callGenes(g)
        then:
        g.proteins*.id == ['sample_00001', 'sample_00002']
        g.proteins*.start == [1, 61]
        g.proteins*.sequence == ['MMMMMMMMMM', 'AAAAAAAAAA']
        g.proteins.every { it.sourceFeature.id == it.id && it.sourceFeature.locusTag == it.id && it.contig.is(g.contigs[0]) }
        status(caller)?.outcome == PredictorOutcome.SUCCESS
        status(caller)?.exitCode == 0
        where:
        type << [PyrodigalCaller, ProdigalCaller]
    }

    def 'malformed or incomplete output never mutates the original genome'() {
        given:
        String gff = prefix + cds
        String proteins = faa
        int exit = 0
        switch (fault) {
            case 'missing protein': proteins = faa.split('>chr_2')[0]; break
            case 'extra protein': proteins += '>unknown\nMMMM\n'; break
            case 'duplicate protein': proteins += faa; break
            case 'duplicate CDS': gff += cds; break
            case 'wrong coordinate': proteins = faa.replace('# 61 #', '# 62 #'); break
            case 'wrong strand': proteins = faa.replace('# -1 #', '# 1 #'); break
            case 'bad GFF': gff += 'broken\n'; break
            case 'missing coverage': gff = '##gff-version 3\n' + cds; break
            case 'wrong coverage': gff = gff.replace('seqlen=120', 'seqlen=119'); break
            case 'unknown contig': gff = gff.replace('chr\t', 'other\t'); break
            case 'native error': exit = 7; break
            case 'garbage FAA': proteins = 'not fasta\n'; gff = prefix; break
        }
        def g = genome()
        def sentinel = new Feature(type: FeatureType.RRNA, start: 100, end: 110)
        g.contigs[0].addFeature(sentinel)
        def caller = new PyrodigalCaller(executablePath: replay(gff, proteins, exit).path)
        when:
        caller.callGenes(g)
        then:
        g.proteinCount == 0
        g.contigs[0].features.size() == 1 && g.contigs[0].features[0].is(sentinel)
        status(caller)?.outcome == PredictorOutcome.FAILED
        status(caller)?.exitCode == exit
        where:
        fault << ['missing protein','extra protein','duplicate protein','duplicate CDS','wrong coordinate',
            'wrong strand','bad GFF','missing coverage','wrong coverage','unknown contig','native error','garbage FAA']
    }

    def 'observed native zero output remains successful with no false gene-absence inference'() {
        given:
        def g = new Genome(id: 'short')
        g.addContig(new Contig(id: 'short', sequence: 'ATGAAATAA'))
        def caller = new PyrodigalCaller(metaMode: true,
            executablePath: replay(fixture('short-meta.gff').text, fixture('short-meta.faa').text).path)
        when:
        caller.callGenes(g)
        then:
        g.proteinCount == 0
        g.contigs[0].features.empty
        status(caller)?.outcome == PredictorOutcome.SUCCESS_ZERO
        status(caller)?.exitCode == 0
        status(caller)?.detail.contains('does not establish')
    }

    def 'observed numeric-contig native metadata does not create false ambiguity'() {
        given:
        File fasta = fixture('numeric.fna'), gff = fixture('numeric-meta.gff'), proteins = fixture('numeric-meta.faa')
        def g = FastaReader.readGenome(fasta, 'sample')
        def caller = new PyrodigalCaller(metaMode: true, executablePath: replay(gff.text, proteins.text).path)
        when:
        def joined = GeneModelReader.read(fasta, gff, proteins)
        caller.callGenes(g)
        then:
        joined.proteinCount == 291
        g.proteinCount == 291
        status(caller)?.outcome == PredictorOutcome.SUCCESS
        g.proteins.collect { [it.contig.id,it.start,it.end,it.strand,it.sequence] } ==
            joined.contigs.collectMany { c -> c.proteins.sort { it.start }.collect { [c.id,it.start,it.end,it.strand,it.sequence] } }
    }

    def 'a skipped sequence is partial and stale success cannot survive the next empty run'() {
        given:
        def g = genome()
        g.addContig(new Contig(id: 'unloaded'))
        def caller = new PyrodigalCaller(executablePath: replay(prefix + cds, faa).path)
        when:
        caller.callGenes(g)
        def first = status(caller)
        def empty = new Genome(id: 'empty')
        empty.addContig(new Contig(id: 'no_sequence'))
        caller.callGenes(empty)
        then:
        first?.outcome == PredictorOutcome.PARTIAL
        first?.skippedInputs == ['unloaded']
        first?.annotationCount == 2
        status(caller)?.outcome == PredictorOutcome.FAILED
        !status(caller).is(first)
        status(caller)?.exitCode == null
    }

    def 'successful process with missing outputs fails instead of reporting zero'() {
        given:
        def caller = new PyrodigalCaller(executablePath: executable('pass\n').path)
        when:
        caller.callGenes(genome())
        then:
        status(caller)?.outcome == PredictorOutcome.FAILED
        status(caller)?.exitCode == 0
    }

    def 'existing CDS is preserved and a missing executable is unavailable'() {
        given:
        def caller = new PyrodigalCaller(executablePath: new File(dir, 'not-installed').path)
        def g = genome()
        when:
        caller.callGenes(g)
        then:
        status(caller)?.outcome == PredictorOutcome.UNAVAILABLE
        !caller.isAvailable()
        caller.unavailableReason()
        when:
        g.contigs[0].addFeature(new Feature(type: FeatureType.CDS, start:1, end:30))
        caller.callGenes(g)
        then:
        status(caller)?.outcome == PredictorOutcome.FAILED
        status(caller)?.detail.contains('existing CDS')
        g.contigs[0].features.size() == 1
    }

    def 'availability requires successful recognizable native version output'() {
        given:
        def caller = new ProdigalCaller(executablePath: executable("print('prodigal fixture')\nraise SystemExit(7)\n").path)
        expect:
        !caller.isAvailable()
        caller.unavailableReason()
        when:
        caller.executablePath = replay(prefix, '').path
        then:
        caller.isAvailable()
        caller.unavailableReason() == null
    }

    @Timeout(10)
    def 'timeout bounds output consumption and terminates descendants holding the log open'() {
        given:
        File pid = new File(dir, 'child.pid')
        def tool = executable("""import subprocess,time
from pathlib import Path
p=subprocess.Popen(['sleep','60'])
Path('${pid.path}').write_text(str(p.pid))
print('before timeout',flush=True)
p.wait()
""")
        def caller = new PyrodigalCaller(executablePath: tool.path, timeoutMillis: 500L)
        long start = System.currentTimeMillis()
        when:
        caller.callGenes(genome())
        then:
        System.currentTimeMillis() - start < 4000
        status(caller)?.outcome == PredictorOutcome.FAILED
        status(caller)?.exitCode == null
        status(caller)?.detail.contains('timed out')
        status(caller)?.detail.contains('before timeout')
        pid.isFile()
        !ProcessHandle.of(pid.text.toLong()).map { it.isAlive() }.orElse(false)
        cleanup:
        if (pid.isFile()) ProcessHandle.of(pid.text.toLong()).ifPresent { it.destroyForcibly() }
    }
}
