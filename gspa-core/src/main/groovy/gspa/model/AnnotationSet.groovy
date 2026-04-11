package gspa.model

import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * A collection of annotations, typically for a single protein or an entire genome.
 * Provides convenience methods for filtering and querying annotations.
 */
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class AnnotationSet {

    List<Annotation> annotations = []

    AnnotationSet() {}

    AnnotationSet(List<Annotation> annotations) {
        this.annotations = new ArrayList<>(annotations)
    }

    void add(Annotation annotation) {
        annotations << annotation
    }

    void addAll(Collection<Annotation> anns) {
        annotations.addAll(anns)
    }

    List<Annotation> byType(AnnotationType type) {
        annotations.findAll { it.type == type }
    }

    List<Annotation> bySource(String source) {
        annotations.findAll { it.source == source }
    }

    List<Annotation> goAnnotations() {
        byType(AnnotationType.GO)
    }

    List<Annotation> ecAnnotations() {
        byType(AnnotationType.EC)
    }

    List<Annotation> pfamAnnotations() {
        byType(AnnotationType.PFAM)
    }

    /** Get GO annotations filtered by aspect (MF, BP, CC) */
    List<Annotation> goByAspect(String aspect) {
        goAnnotations().findAll { it.goAspect == aspect }
    }

    /** Get unique GO term IDs */
    Set<String> goTermIds() {
        goAnnotations().collect { it.value } as Set
    }

    /** Get unique annotation values of any type */
    Set<String> valuesOfType(AnnotationType type) {
        byType(type).collect { it.value } as Set
    }

    /** Get all distinct sources that contributed annotations */
    Set<String> sources() {
        annotations.collect { it.source } as Set
    }

    /** Filter annotations above a confidence threshold */
    AnnotationSet aboveThreshold(double threshold) {
        new AnnotationSet(annotations.findAll { it.score >= threshold })
    }

    int size() { annotations.size() }
    boolean isEmpty() { annotations.isEmpty() }

    /**
     * Merge with another annotation set, deduplicating by (type, value).
     * When duplicates exist, keep the one with higher score.
     */
    AnnotationSet merge(AnnotationSet other) {
        Map<String, Annotation> best = [:]
        (this.annotations + other.annotations).each { ann ->
            String key = "${ann.type}:${ann.value}"
            Annotation existing = best[key]
            if (!existing || ann.score > existing.score) {
                best[key] = ann
            }
        }
        new AnnotationSet(best.values().toList())
    }
}
