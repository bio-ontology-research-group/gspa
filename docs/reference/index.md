# Command-line reference

GSPA is distributed as a single executable jar (`gspa-2.1.0.jar`) plus a small
Bash launcher, `gspa`, at the root of the repository. The launcher locates the
jar, sets JVM options and forwards all arguments to it. One subcommand,
`gspa agent`, is handled by the launcher itself and runs a Python program
instead of the jar.

## Subcommands

| Subcommand | Implemented in | Purpose |
|---|---|---|
| [`gspa annotate`](cli-annotate.md) | jar | Predict protein functions for a genome (FASTA, GFF3 or protein FASTA), optionally compute genome-scale quality metrics and enforce constraints. |
| [`gspa integrate`](cli-integrate.md) | jar | Combine per-protein evidence claims into posterior probabilities with genome-scale priors (Stage 1), and optionally run symbolic enforcement (Stage 2A), flux-certified refinement (Stage 2B), gap suggesters and agentic completion (Stage 3). |
| [`gspa evaluate`](cli-evaluate.md) | jar | Compute completeness, coherence, consistency and information content for an existing annotation (GFF3 + GAF). |
| [`gspa compare`](cli-compare.md) | jar | Evaluate several annotation sets of the same genome and write a side-by-side report. |
| [`gspa report`](cli-report.md) | jar | Render HTML from one or more JSON quality reports. |
| [`gspa visualize`](cli-visualize.md) | jar (calls Python) | Build a self-contained HTML browser for a GSPA workspace. |
| [`gspa agent`](cli-agent.md) | launcher (Python) | Run the Stage 3 agent sidecar directly on an anomalies file. |

Every jar subcommand accepts `-h`/`--help` and `-V`/`--version`. Running `gspa`
without a subcommand prints the top-level usage.

```bash
./gspa --version          # prints: gspa 2.1.0
./gspa integrate --help
```

The version string is read from `version.properties`, which Gradle fills in at
build time.

## How the launcher works

The launcher is a Bash script (`gspa`) that performs these steps in order.

1. **Agent dispatch.** If the first argument is `agent`, the launcher drops it
   and executes the Stage 3 sidecar with the remaining arguments. It looks for
   `benchmark/agent/run_gspa_agent.py` next to the launcher, then
   `run_gspa_agent.py` next to the launcher, and exits with status 1 if neither
   exists. The script is run with `uv run python` when `uv` is on `PATH`,
   otherwise with `python3`. The jar is not involved. See
   [`gspa agent`](cli-agent.md).
2. **Jar lookup.** For every other subcommand the launcher picks the first jar
   that exists from this list:
    1. the file named by the environment variable `GSPA_JAR`;
    2. `gspa-2.1.0.jar` next to the launcher;
    3. `gspa-cli/build/libs/gspa-2.1.0.jar` (the Gradle build output);
    4. if `./gradlew` exists and is executable, the launcher builds the jar with
       `./gradlew :gspa-cli:shadowJar --quiet` and uses the result.

    If none of these applies it prints an error and exits with status 1.
3. **Java and memory.** The Java binary is `$JAVA_HOME/bin/java` when
   `JAVA_HOME` is set, otherwise `java` from `PATH`. JVM options come from
   `GSPA_JVM_OPTS`; the default is `-Xmx8g` (8 GB maximum heap).
4. **Execution.** The launcher `exec`s
   `java $GSPA_JVM_OPTS -jar <jar> <arguments>`.

```bash
# Use a downloaded release jar and give the JVM 32 GB of heap
export GSPA_JAR=/opt/gspa/gspa-2.1.0.jar
export GSPA_JVM_OPTS="-Xmx32g"
./gspa annotate -i genome.fna -k bacteria -o out/
```

!!! note "Memory"
    `GSPA_JVM_OPTS` replaces the default entirely, so include `-Xmx` whenever
    you set it. Loading the full GO OWL file with the ELK reasoner (process
    coherence, `--enforce-coherence`) is the most memory-intensive step;
    `--lite` in `integrate` and `evaluate` skips the reasoner.

You can always bypass the launcher and call the jar directly:

```bash
java -Xmx16g -jar gspa-cli/build/libs/gspa-2.1.0.jar integrate --help
```

All environment variables are listed in [Environment variables](environment.md).

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success. Also returned by `compare` and `report` when no positional input files are given (they print an error message but do not fail). |
| 1 | Launcher could not find the jar or the agent script; agent sidecar could not find its `--anomalies` file. An uncaught exception in a jar subcommand (for example a missing required input in `integrate`, an invalid `--promotion-strategy`, or `visualize` failing because its Python script returned non-zero) also ends with status 1 (picocli's execution-error code). |
| 2 | Command-line usage error reported by the argument parser: an unknown option, a missing required option (for example `--go-owl` in `evaluate`), or a value that cannot be converted to the option type. The Python agent sidecar uses the same code for argparse errors. |

Warnings about optional inputs that fail to load (GO ontology, pathway files,
reaction graph, essential-function profiles) are printed to standard error as
`[warn]` lines and do not change the exit code; the corresponding prior or
feature is then skipped.

## Pages in this section

- One page per subcommand: [integrate](cli-integrate.md),
  [annotate](cli-annotate.md), [evaluate](cli-evaluate.md),
  [compare](cli-compare.md), [report](cli-report.md),
  [visualize](cli-visualize.md), [agent](cli-agent.md).
- [Configuration files](configuration.md): the YAML schema used by
  `gspa annotate -c` and the bundled kingdom presets.
- [Hyperparameters (theta)](theta.md): the JSON file passed to
  `gspa integrate --theta`.
- [Input and output formats](file-formats.md): every file GSPA reads or writes.
- [Bundled resources](bundled-resources.md): ontology constraint tables and
  profiles shipped inside the jar.
- [Environment variables](environment.md).
- [Nextflow parameters](nextflow-params.md) for the `gspa-nf` pipeline.
