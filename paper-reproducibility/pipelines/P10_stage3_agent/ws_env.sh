# Source on ws (leechuck-office, 56 cores) to run P10 steps against the ~/p10_mirror copy of /data.
# Used on 2026-10-19 because the unimatrix01 debug nodes were fully allocated by other jobs
# (node007 additionally lacks the /data mount). Same jar, same JDK build, same inputs (byte copies).
export M=$HOME/p10_mirror
export JAVA=$M/data/hohndor/envs/jdk-21.0.10+7/bin/java
export JAR=$M/data/hohndor/gspa/bin/gspa-unified.jar
export REF=$M/data/hohndor/gspa/reference
export P10_PIPE=$M/data/hohndor/gspa/proteomes/pillar4_repro_2026-10-19/code
export P10_ROOT=$M/data/hohndor/gspa/proteomes/pillar4_repro_2026-10-19
export XMX=${XMX:-6g}
export P10_CLEAN_CLAIMS=1
sed "s#/data/#$M/data/#g" "$P10_PIPE/genomes_v2.tsv" > "$P10_PIPE/genomes_v2_node.tsv"
