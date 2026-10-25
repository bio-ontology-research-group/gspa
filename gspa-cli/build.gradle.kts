plugins {
    groovy
    application
    id("com.github.johnrengelman.shadow") version "8.1.1"
}

application {
    mainClass.set("gspa.cli.GspaMain")
}

dependencies {
    implementation(project(":gspa-core"))
    implementation(project(":gspa-predictors"))
    implementation("info.picocli:picocli-groovy:${property("picocliVersion")}")
    runtimeOnly("ch.qos.logback:logback-classic:${property("logbackVersion")}")
}

tasks.named<ProcessResources>("processResources") {
    val v = project.version.toString()
    inputs.property("version", v)
    filesMatching("version.properties") {
        filter(org.apache.tools.ant.filters.ReplaceTokens::class, "tokens" to mapOf("version" to v))
    }
}

// A stdlib-only Python zipapp keeps the standalone jar's downloader identical
// to scripts/gspa_resources.py and includes small models/reference tables offline.
val resourceManager by tasks.registering(Zip::class) {
    archiveFileName.set("gspa-resources.pyz")
    destinationDirectory.set(layout.buildDirectory.dir("resource-manager"))
    from(rootProject.file("scripts/gspa_resources.py")) { rename { "__main__.py" } }
    from(rootProject.file("gspa-core/src/main/resources/downloads/catalog.json"))
    from(rootProject.file("gspa-core/src/main/resources")) {
        include("taxon-constraints/**", "anatomical-constraints/**", "complex-terms/**", "essential-functions/**")
        into("bundled")
    }
    from(rootProject.file("deepgo-plusplus/models")) { include("*.json"); into("bundled/models") }
    from(rootProject.file("benchmark/neural")) {
        include("run_neural_predictors.py", "run_region_predictors.py", "run_term_predictors.py",
                "run_site_predictors.py", "run_genomic_predictors.py")
        into("bundled/sidecars")
    }
}

tasks.named<ProcessResources>("processResources") {
    from(resourceManager) { into("downloads") }
}

tasks.named<JavaExec>("run") {
    standardInput = System.`in`
}

tasks.named<com.github.jengelman.gradle.plugins.shadow.tasks.ShadowJar>("shadowJar") {
    archiveBaseName.set("gspa")
    archiveClassifier.set("")
    archiveVersion.set(project.version.toString())
    mergeServiceFiles()
}
