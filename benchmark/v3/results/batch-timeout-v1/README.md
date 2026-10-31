# Verified-batch timeout process cleanup

The failing regression launches a Python wrapper, a Python engine and a child
that ignores SIGTERM. The old subprocess timeout kills only the wrapper and
leaves the child running. The repaired runner launches a separate process group,
terminates that group, allows five seconds for shutdown, then sends SIGKILL to
any remaining members. It preserves TimeoutExpired and writes no success receipt.
Processes that deliberately detach into a new session are outside this guarantee.

All 23 focused runner/capacity tests pass; the full benchmark run passes 525 tests
with three optional skips, including eight live Nextflow tests. The exact old
runner and failing log are retained. Historical receipts bind historical runner
hashes and must be checked using their archived sources.

The capacity driver now accepts a positive per-arm inference deadline. Its result
records distinguish that deadline from total elapsed time including verification.
The frozen launch in launch/ submitted Slurm job 10058: DeepProZyme full-matrix
stage1, 9,790,188 candidates, 64 GiB heap, 80 GiB allocation, four CPUs, four-hour
inference deadline, four-and-a-half-hour allocation. It retains the original
c469273224723bff080efd9440d1c1ae07fc0c943c2cacffa04a98003481f801 JAR to
match the completed baseline. Submission/running is not capacity acceptance.
At the observation recorded in job-observation.json it was still running; results are pending.
