plugins {
    groovy
}

allprojects {
    group = "org.gspa"
    version = "1.3.0"

    repositories {
        mavenCentral()
    }
}

subprojects {
    apply(plugin = "groovy")

    java {
        sourceCompatibility = JavaVersion.VERSION_21
        targetCompatibility = JavaVersion.VERSION_21
    }

    dependencies {
        implementation("org.apache.groovy:groovy:${property("groovyVersion")}")
        implementation("org.slf4j:slf4j-api:${property("slf4jVersion")}")

        testImplementation("org.spockframework:spock-core:${property("spockVersion")}")
        testRuntimeOnly("ch.qos.logback:logback-classic:${property("logbackVersion")}")
    }

    tasks.withType<Test> {
        useJUnitPlatform()
    }
}

// ---------------------------------------------------------------------------
// Convenience tasks for the Nextflow sibling pipeline (gspa-nf/).
// These just shell out to `nextflow`; they skip cleanly if Nextflow isn't on
// PATH. The pipeline itself is NOT a Gradle subproject — it has no sources
// we compile and no tests we run; see gspa-nf/README.md.
// ---------------------------------------------------------------------------

fun hasNextflow(): Boolean = try {
    ProcessBuilder("sh", "-c", "command -v nextflow")
        .redirectErrorStream(true)
        .start().waitFor() == 0
} catch (_: Exception) { false }

tasks.register<Exec>("nfLint") {
    group = "nextflow"
    description = "Parse + lint gspa-nf/main.nf (requires nextflow on PATH)."
    workingDir = rootDir
    commandLine("nextflow", "inspect", "gspa-nf/main.nf")
    isIgnoreExitValue = !hasNextflow()
    onlyIf { hasNextflow() }
}

tasks.register("nfHelp") {
    group = "nextflow"
    description = "Print the gspa-nf quickstart (no Nextflow required)."
    doLast {
        println(file("gspa-nf/README.md").readText().lineSequence().take(30).joinToString("\n"))
    }
}
