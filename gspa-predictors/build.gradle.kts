plugins {
    groovy
    `java-library`
}

dependencies {
    api(project(":gspa-core"))
    implementation("org.apache.groovy:groovy-json:${property("groovyVersion")}")
}
