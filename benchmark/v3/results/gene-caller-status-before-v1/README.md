# Gene-caller pre-repair evidence

Four Slurm jobs exited zero: 10507 native Pyrodigal3.7.1 outcomes; 10509 unchanged production JAR/shared AnnotationPipeline (six calls); 10510 native calls on two real contigs renamed 2 and 1; 10511 unchanged strict reader rejecting those native outputs. No production repair or accuracy improvement is claimed.

Native and controlled cases are distinct: short/all-N/empty inputs and numeric header renaming test contracts, not biological gene truth. Reordering/removing a captured protein FASTA record is explicit output corruption, not native Pyrodigal behavior. Full original assembly calls remain correct. No native Prodigal binary was on node004 PATH.

Production JAR: /data/hohndor/gspa-v3/antismash-nextflow-v1/jobs/10481/gspa-3.0.0.jar, SHA256 069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303. Ten relevant sources match its 828-file frozen source manifest. JAR/runtime and wheels are not duplicated here; commands and wheel/source hashes remain in the retained scripts and observations. The inherited finish hook of 10509 also copied an unrelated base-runtime JAR; it was NOT used by probe.py and its jar.sha256 is deliberately omitted from this archive. Probe.py verifies the actual before.jar hash before every group of calls.

Files above 256 KiB are gzip-compressed. retained-bindings.json maps original relative names to retained paths and decompressed hashes. Run `python3 verify_archive.py` for independent hash, sequence/coordinate, count and status checks. No test suite was compiled/executed in these jobs. Groovy harnesses invoke the existing production classes.
