package gspa.model

import groovy.transform.Canonical

/**
 * A cluster of proteins grouped by similarity (typically MMseqs2 linclust
 * at a given identity / coverage threshold).
 *
 * The representative is always included in {@link #members}.
 */
@Canonical
class ProteinCluster {
    String clusterId
    ProteinRef representative
    List<ProteinRef> members = []

    boolean isSingleton() { members.size() <= 1 }

    int size() { members.size() }
}
