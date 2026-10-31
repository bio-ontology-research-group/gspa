# Supplied GFF3/protein route repair

Complete native inputs, four before-repair and ten after-repair CLI outputs,
845 final frozen sources, baseline sources, build/test XML, scripts and reviews.
Run `python3 verify_archive.py` here. Gzip files preserve exact original bytes;
`retained.json` binds their decompressed contents. Final `archive.sha256` binds
all retained files. JARs and Java runtime are omitted and explicitly hash-listed
by the verifier; the JARs remain in the remote job folders. Source tarballs and
base runtime packages remain remote, with package hashes retained. Generated
Gradle binary result stores are omitted; test XML and logs are retained.

This archive proves identity/coordinate/neighborhood preservation, not biological
function or EC-prediction accuracy. See `docs/benchmarks/v3-supplied-gene-models.md`.
