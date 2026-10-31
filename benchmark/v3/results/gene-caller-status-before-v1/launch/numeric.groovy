import com.fasterxml.jackson.databind.ObjectMapper
import gspa.cli.AntiSmashClaimsCommand
File input = new File(args[0]); File out = new File(args[1]); out.mkdirs()
def result = [:]
['single','meta'].each { mode ->
    try {
        def genome = AntiSmashClaimsCommand.loadCoordinates(new File(input,'inputs/numeric.fna'),new File(input,"numeric-${mode}.gff"),new File(input,"numeric-${mode}.faa"))
        result[mode] = [accepted:true,proteins:genome.proteinCount]
    } catch (IllegalArgumentException e) {
        result[mode] = [accepted:false,error:e.message]
    }
}
new File(out,'numeric-reader.json').text = new ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString(result)+'\n'
assert result.values().any { !it.accepted }
