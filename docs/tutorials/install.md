# Install GSPA from a clean container

*2026-10-24T13:49:24Z by Showboat 0.6.1*
<!-- showboat-id: a7494ad9-c0ee-4dd7-a66a-e3f016d83ba9 -->

This recorded walkthrough builds the current source tree inside an unconfigured Ubuntu 24.04 image. It uses no host Java installation, compiled JAR, reference database or Python environment. Run commands from the GSPA checkout root. Docker and Python 3 are host prerequisites; install uv only if you want to regenerate or verify the Showboat documents. The source snapshot includes the usability changes in this checkout, which may not yet be in a published release. See tutorials/workflows/Dockerfile for every package and build step.

```bash
set -eu
# Change this path to a disk with at least 25 GiB free, if necessary.
mkdir -p .tutorial-work
python3 tutorials/workflows/prepare_context.py .tutorial-work
```

```output
Build inputs copied; no host JAR, databases or Gradle caches included.
Storage check: at least 25 GiB available.
```

The CPU image installs Java 21, Python, Prodigal and DIAMOND, then builds GSPA with a fresh Gradle cache. The full build log is saved outside the document because download speeds and progress messages vary. A failed build still fails this command.

```bash
set -eu
docker build -t gspa-tutorial:cpu -f .tutorial-work/source/tutorials/workflows/Dockerfile .tutorial-work/source > .tutorial-work/build-cpu.log 2>&1 || { tail -60 .tutorial-work/build-cpu.log; exit 1; }
docker run --rm gspa-tutorial:cpu gspa --version
docker run --rm gspa-tutorial:cpu diamond version
```

```output
gspa 2.2.2
diamond version 2.1.9
```

The image is now usable without a GPU. Its launcher limits Java to a 4 GiB heap and two visible CPUs for these small examples. Containers below are also limited to 6 GiB RAM and two CPUs, and run with your user ID so generated files remain editable. This tutorial image includes build tools for transparency; it is not the larger optional production predictor image in docker/Dockerfile.

```bash
set -eu
docker run --rm gspa-tutorial:cpu bash -euc 'python3 --version; prodigal -v 2>&1; command -v gspa'
```

```output
Python 3.12.3

Prodigal V2.6.3: February, 2016

/opt/gspa/gspa
```

Continue with [resources and teaching data](resources.md), then choose a biological workflow. For neural execution, follow [CPU and NVIDIA GPU setup](neural.md). Keep `.tutorial-work/source-manifest.json` with your results: it records the base commit and hashes of the exact working-tree source snapshot. Keep the image for fast replay; it avoids relying on unchanged package repositories. The recorded OS packages are Ubuntu 24.04 amd64 packages; other architectures or future package versions may change the transcript.

Verify with `uvx showboat verify docs/tutorials/install.md`. Verification rebuilds from source if necessary and checks captured outputs; Docker may reuse unchanged layers. To test installation completely from scratch, use a fresh Docker builder/cache and an empty tutorial work directory.
