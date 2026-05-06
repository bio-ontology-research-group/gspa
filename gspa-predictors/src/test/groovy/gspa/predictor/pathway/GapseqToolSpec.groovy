package gspa.predictor.pathway

import spock.lang.Specification

/**
 * CLI-shape tests for the {@code tool} field (gapsmith vs gapseq).
 * Verifies the two tools produce the correct, incompatible flag syntax
 * so downstream runs don't silently misparse.
 */
class GapseqToolSpec extends Specification {

    def "default tool is gapsmith"() {
        given:
        def gp = new GapseqPredictor()

        expect:
        gp.tool == 'gapsmith'
        gp.name == 'gapsmith'
        gp.executable == 'gapsmith'
    }

    def "gapsmith buildCommand uses -t/-o and places fasta last"() {
        given:
        def gp = new GapseqPredictor()
        def inFa = new File('/tmp/x.faa')
        def outDir = new File('/tmp/out')

        when:
        def cmd = gp.buildCommand(inFa, outDir)

        then:
        cmd[0] == 'gapsmith'
        cmd[1] == 'find'
        cmd.contains('-p') && cmd[cmd.indexOf('-p') + 1] == 'all'
        cmd.contains('-t') && cmd[cmd.indexOf('-t') + 1] == 'Bacteria'
        cmd.contains('-o') && cmd[cmd.indexOf('-o') + 1] == outDir.absolutePath
        cmd.last() == inFa.absolutePath
        !cmd.contains('-b')
    }

    def "gapsmith --data-dir is threaded in at the top level"() {
        given:
        def gp = new GapseqPredictor(dataDir: '/ref/gapsmith-data')

        when:
        def cmd = gp.buildCommand(new File('/tmp/x.faa'), new File('/tmp/out'))

        then:
        cmd[0] == 'gapsmith'
        cmd[1] == '--data-dir'
        cmd[2] == '/ref/gapsmith-data'
        cmd[3] == 'find'
    }

    def "tool=gapseq falls back to legacy -b flag and no -o"() {
        given:
        def gp = new GapseqPredictor(tool: 'gapseq')

        expect:
        gp.name == 'gapseq'
        gp.executable == 'gapseq'

        when:
        def cmd = gp.buildCommand(new File('/tmp/x.faa'), new File('/tmp/out'))

        then:
        cmd[0] == 'gapseq'
        cmd.contains('-b') && cmd[cmd.indexOf('-b') + 1] == 'Bacteria'
        !cmd.contains('-t')
        !cmd.contains('-o')
        cmd.last() == '/tmp/x.faa'
    }

    def "explicit executablePath overrides the binary name"() {
        given:
        def gp = new GapseqPredictor(executablePath: '/opt/gapsmith-0.2/gapsmith')

        when:
        def cmd = gp.buildCommand(new File('/tmp/x.faa'), new File('/tmp/out'))

        then:
        cmd[0] == '/opt/gapsmith-0.2/gapsmith'
    }
}
