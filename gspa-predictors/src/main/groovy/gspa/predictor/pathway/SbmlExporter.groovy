package gspa.predictor.pathway

import gspa.model.AnnotationType
import gspa.model.Genome
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Exports metabolic models in SBML format.
 * Creates a basic SBML Level 3 model from EC annotations.
 *
 * For full FBA, the SBML can be loaded into COBRApy.
 */
class SbmlExporter {

    private static final Logger log = LoggerFactory.getLogger(SbmlExporter)

    /**
     * Export a genome's metabolic annotations as an SBML model.
     */
    static void exportSbml(Genome genome, File outputFile) {
        def ecNumbers = genome.proteins.collectMany { protein ->
            protein.annotations.byType(AnnotationType.EC).collect { ann ->
                [proteinId: protein.id, ec: ann.value, score: ann.score]
            }
        }

        log.info("Exporting SBML for ${genome.id}: ${ecNumbers.size()} reactions")

        outputFile.withWriter { writer ->
            writer.writeLine('<?xml version="1.0" encoding="UTF-8"?>')
            writer.writeLine('<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">')
            writer.writeLine("  <model id=\"${xmlEsc(genome.id)}\" name=\"GSPA metabolic model for ${xmlEsc(genome.id)}\">")

            // Compartments
            writer.writeLine('    <listOfCompartments>')
            writer.writeLine('      <compartment id="c" name="cytoplasm" constant="true"/>')
            writer.writeLine('      <compartment id="e" name="extracellular" constant="true"/>')
            writer.writeLine('    </listOfCompartments>')

            // Reactions from EC annotations
            if (ecNumbers) {
                writer.writeLine('    <listOfReactions>')
                ecNumbers.eachWithIndex { entry, idx ->
                    String rxnId = "R_${entry.ec.replaceAll(/[^a-zA-Z0-9]/, '_')}_${idx}"
                    writer.writeLine("      <reaction id=\"${rxnId}\" name=\"${xmlEsc(entry.ec)}\" reversible=\"true\">")
                    writer.writeLine("        <notes><body xmlns=\"http://www.w3.org/1999/xhtml\">")
                    writer.writeLine("          <p>EC: ${xmlEsc(entry.ec)}</p>")
                    writer.writeLine("          <p>Gene: ${xmlEsc(entry.proteinId)}</p>")
                    writer.writeLine("          <p>Confidence: ${entry.score}</p>")
                    writer.writeLine("        </body></notes>")
                    writer.writeLine('      </reaction>')
                }
                writer.writeLine('    </listOfReactions>')
            }

            writer.writeLine('  </model>')
            writer.writeLine('</sbml>')
        }

        log.info("SBML model written to: ${outputFile}")
    }

    /**
     * Run FBA via COBRApy on exported SBML models.
     * Requires Python with cobra installed.
     */
    static Map<String, Object> runFba(File sbmlFile, String cobraScript = null) {
        log.info("Running FBA on: ${sbmlFile}")

        String script = cobraScript ?: """
import cobra
import json
import sys

model = cobra.io.read_sbml_model(sys.argv[1])
solution = model.optimize()
result = {
    'objective_value': solution.objective_value if solution.status == 'optimal' else None,
    'status': solution.status,
    'fluxes': {r.id: solution.fluxes[r.id] for r in model.reactions} if solution.status == 'optimal' else {}
}
print(json.dumps(result))
"""

        def tmpScript = File.createTempFile('gspa_fba_', '.py')
        try {
            tmpScript.text = script
            def proc = ['python3', tmpScript.absolutePath, sbmlFile.absolutePath].execute()
            proc.waitFor(30, java.util.concurrent.TimeUnit.MINUTES)

            if (proc.exitValue() == 0) {
                def result = new groovy.json.JsonSlurper().parseText(proc.text)
                log.info("FBA result: status=${result.status}, objective=${result.objective_value}")
                return result as Map
            } else {
                log.error("FBA failed: ${proc.errorStream.text}")
                return [status: 'error', error: proc.errorStream.text]
            }
        } finally {
            tmpScript.delete()
        }
    }

    private static String xmlEsc(String s) {
        s?.replaceAll('&', '&amp;')?.replaceAll('<', '&lt;')?.replaceAll('>', '&gt;')?.replaceAll('"', '&quot;') ?: ''
    }
}
