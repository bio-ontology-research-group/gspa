plugins {
    groovy
    `java-library`
}

dependencies {
    // Ontology
    api("net.sourceforge.owlapi:owlapi-distribution:${property("owlapiVersion")}")
    api("org.semanticweb.elk:elk-owlapi:${property("elkVersion")}")

    // Bioinformatics
    api("org.biojava:biojava-core:${property("biojavaVersion")}")

    // Graph algorithms (pathway coherence)
    api("org.jgrapht:jgrapht-core:${property("jgraphtVersion")}")

    // SAT solver (consistency checking + weighted MaxSAT promotion strategy)
    api("org.ow2.sat4j:org.ow2.sat4j.core:${property("sat4jVersion")}")
    api("org.ow2.sat4j:org.ow2.sat4j.maxsat:${property("sat4jVersion")}")

    // Configuration
    implementation("com.fasterxml.jackson.core:jackson-databind:${property("jacksonVersion")}")
    implementation("com.fasterxml.jackson.dataformat:jackson-dataformat-yaml:${property("jacksonVersion")}")

    // Parallel processing
    implementation("org.codehaus.gpars:gpars:1.2.1")

    // JSON support for reports
    implementation("org.apache.groovy:groovy-json:${property("groovyVersion")}")
}
