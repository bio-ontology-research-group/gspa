# Sonnet review disposition

Both reviews used Claude Code `claude-sonnet-5-5` in read-only plan mode.
No reviewer approval flow was adopted; implementation and benchmark work were
already authorized. Reviewer-generated plan files outside the repository do not
control this task.

Accepted: complete-collection finalizer; source/JAR/parser/native-library hashing;
immutable attempts; strict small/stream metadata parity; learned-report size guard;
cell-identity checks; process-level RSS measurement; respecting reserved node RAM;
explicit sparse-wave and single-arm diagnostic scope.

The finalizer requires every registered cell and both arms, checks artifact/input
hashes, and refuses missing or duplicate completions. The small fallback rejects
nonfinite JSON and duplicate root/stage metadata; duplicates inside ignored change
records are ignored consistently by both parsers. Installed parser sources and
bundled `ijson.libs/libyajl-*.so.*` are bound; final replay used yajl2_c.

Not adopted: requiring the C backend in the general verifier. A slower supported
backend is not itself a correctness failure or proof of unbounded retention.
Backend identity is recorded; the actual benchmark uses the C backend. The unit
memory test measures Python accumulation only; the full 2.68GB remote replay
measures total RSS, so a separate large synthetic allocation test adds little.

A `finally` cannot reliably write a terminal marker after SIGKILL or node loss.
Such attempts remain incomplete, and the finalizer refuses them. Scheduler/log
state must establish liveness and termination. Stage0/stage1 prior flags are
checked in finalization. Exception paths preserve partial results.

Initial 240G allocation proposal was reduced to 192G after checking node005's
48G service reservation; heap is 168g. An OOM would demonstrate failure at that
limit, not impossibility on all hardware. Sparse66 covers only 0.896% of claims;
it is useful production execution but cannot stand in for all132 or ensembles.
