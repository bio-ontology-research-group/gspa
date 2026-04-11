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

tasks.named<JavaExec>("run") {
    standardInput = System.`in`
}

tasks.named<com.github.jengelman.gradle.plugins.shadow.tasks.ShadowJar>("shadowJar") {
    archiveBaseName.set("gspa")
    archiveClassifier.set("")
    archiveVersion.set(project.version.toString())
    mergeServiceFiles()
}
