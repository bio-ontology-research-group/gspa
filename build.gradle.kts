plugins {
    groovy
}

allprojects {
    group = "org.gspa"
    version = "0.1.0-SNAPSHOT"

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
