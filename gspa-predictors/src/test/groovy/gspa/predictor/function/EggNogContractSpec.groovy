package gspa.predictor.function

import groovy.json.JsonSlurper
import gspa.config.ConfigLoader
import gspa.model.AnnotationType
import gspa.predictor.AnnotationPipeline
import gspa.predictor.ToolOutputContractException
import spock.lang.Specification
import spock.lang.TempDir
import java.nio.file.Path
import java.security.MessageDigest

class EggNogContractSpec extends Specification {
    @TempDir Path tmp

    static class Versioned extends EggNogMapperPredictor {
        String reported
        @Override String getVersion() { reported }
    }

    def 'version-specific commands preserve v2 donor filtering and allow v3 curated donors'() {
        given:
        def p = new Versioned(reported: version, threads: 2, dataDir: '/data/egg')
        when:
        def command = p.buildCommand(tmp.resolve('input.faa').toFile(), tmp.toFile())
        then:
        command.contains('--go_evidence') == filter
        command[command.indexOf('--cpu') + 1] == '2'
        command[command.indexOf('--data_dir') + 1] == '/data/egg'
        where:
        version | filter
        'emapper-2.1.13' | true
        'emapper-3.0.0-beta6' | false
    }

    def 'unknown versions and unsupported v3 modes fail before execution'() {
        when:
        new Versioned(reported: version, searchMode: mode).buildCommand(tmp.resolve('input').toFile(), tmp.toFile())
        then:
        thrown(IllegalArgumentException)
        where:
        version | mode
        'unknown' | 'diamond'
        'emapper-4.0.0' | 'diamond'
        'emapper-3.0.0' | 'hmmer'
    }

    def 'version probe reads beyond database warnings and does not guess on failed execution'() {
        given:
        File tool = tmp.resolve('emapper').toFile()
        tool.text = "#!/bin/sh\nprintf 'database warning\\nemapper-3.0.0-beta6 / Expected DB: 7\\n'\nexit ${exit}\n"
        tool.setExecutable(true)
        def p = new EggNogMapperPredictor(executablePath: tool.absolutePath)
        expect:
        p.version == expected
        where:
        exit | expected
        0 | 'emapper-3.0.0-beta6'
        1 | 'unknown'
    }

    def 'typed configuration reaches eggNOG registration'() {
        given:
        def cfg = ConfigLoader.loadFromString('''predictors:
  similarity: {enabled: false}
  domains: {enabled: false}
  operons: {enabled: false}
  pathway: {enabled: false}
  enable: [eggnog-mapper]
  eggnog: {threads: 2, dataDir: /data/egg, executablePath: /tools/emapper.py, searchMode: mmseqs}
''')
        def pipeline = new AnnotationPipeline(cfg)
        when:
        pipeline.configure()
        def p = pipeline.registry.get('eggnog-mapper')
        then:
        p.threads == 2
        p.dataDir == '/data/egg'
        p.executablePath == '/tools/emapper.py'
        p.searchMode == 'mmseqs'
    }

    def 'observed v2 and v3 tables preserve identities and uncalibrated score provenance'() {
        given:
        File dir = new File(getClass().getResource("/observed-output/eggnog/v${major}/provenance.json").toURI()).parentFile
        def provenance = new JsonSlurper().parse(new File(dir, 'provenance.json'))
        provenance.artifacts.each { name, hash ->
            assert MessageDigest.getInstance('SHA-256').digest(new File(dir, name).bytes).encodeHex().toString() == hash
        }
        when:
        def result = new EggNogMapperPredictor().parseOutput(dir)
        then:
        result.size() == proteins
        result.values().flatten().every { it.metadata.toolVersion.startsWith("emapper-${major}.") &&
            it.metadata.statisticType == 'evalue' && it.metadata.scoreCalibrated == false &&
            it.metadata.seedOrtholog && it.metadata.rawStatistic != null && it.metadata.rawBitScore != null }
        result.values().flatten().any { it.type == AnnotationType.EC }
        where:
        major | proteins
        2 | 2
        3 | 6
    }

    def 'v3 cannot silently fall back to the v2 headerless column positions'() {
        given:
        def source = getClass().getResource('/observed-output/eggnog/v3/eggnog_results.emapper.annotations').text
        tmp.resolve('eggnog_results.emapper.annotations').toFile().text = source.readLines().findAll {
            !it.startsWith('#query') && (keepVersion || !it.startsWith('#')) }.join('\n')
        when:
        new EggNogMapperPredictor().parseOutput(tmp.toFile())
        then:
        thrown(ToolOutputContractException)
        where:
        keepVersion << [true, false]
    }

    def 'optional provenance follows named columns and never treats E-value as seed identity'() {
        given:
        String header = '#query\tevalue\tscore\tGOs\tEC\tKEGG_ko\teggNOG_OGs\tCOG_category'
        String row = 'p\t1e-15\t50\tGO:0005524\t1.2.3.4\tko:K00001\tCOG001\tJ'
        tmp.resolve('eggnog_results.emapper.annotations').toFile().text = header + extraHeader + '\n' + row + extraRow + '\n'
        when:
        def annotations = new EggNogMapperPredictor().parseOutput(tmp.toFile()).p
        then:
        annotations.every { it.metadata.seedOrtholog == seed && it.score == 0.5d }
        annotations.every { it.metadata.annotationConfidence == confidence }
        where:
        extraHeader | extraRow | seed | confidence
        '' | '' | '' | null
        '\tseed_ortholog\tannotation_confidence' | '\t12345\t1,2,0' | '12345' | '1,2,0'
    }

    def 'duplicate score column names fail instead of choosing an arbitrary score'() {
        given:
        tmp.resolve('eggnog_results.emapper.annotations').toFile().text = '#query\tevalue\tevalue\tscore\tGOs\tEC\tKEGG_ko\teggNOG_OGs\tCOG_category\n'
        when:
        new EggNogMapperPredictor().parseOutput(tmp.toFile())
        then:
        thrown(ToolOutputContractException)
    }
}
