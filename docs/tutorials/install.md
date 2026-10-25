# Prepare your first biological analysis

The first lesson asks which *E. coli* proteins might participate in threonine
biosynthesis. This page prepares the software and example sequences; you only
need to do it once for all seven biological lessons.

## Before you start

Use a GSPA checkout on Linux or macOS with **Python 3 and a running Docker
installation**. Docker must be usable by your normal account. Allow at least
25 GiB free on the workspace and Docker storage disks, and make at least 6 GiB
RAM available to Docker. No GPU is needed.

If you do not already have the code:

```bash
git clone https://github.com/bio-ontology-research-group/gspa.git
cd gspa
```

Run tutorial commands from this directory. The launcher and helpers must be
present in your checkout; older releases may not include them.

## Set up once

```bash
./tutorial setup
```

The launcher builds GSPA and its sequence tools in a CPU container, downloads
about 40 MB of teaching sequences and annotations, checks their hashes, and builds
a small search database. The first setup needs internet and can take several
minutes or longer, depending on downloads and Docker build speed. Progress and
log locations appear in the terminal.

When it finishes, you should see:

```text
Ready. Start the first lesson: ./tutorial run bacteria
```

Now go to [the bacterial lesson](bacteria.md). You do not need to download any
additional models or run the resource catalogue first.

## If something goes wrong

| Message or symptom | Next step |
|---|---|
| Docker is not available | Start Docker; check that `docker info` works for your account |
| Need at least 25 GiB free | Free space on the workspace disk; also check Docker's storage disk |
| Build or download failed | Read the log path printed by the launcher, fix the reported issue, then repeat `./tutorial setup` |
| Source changed / checksum mismatch | Keep the error and existing files; use the maintainer's [data-refresh procedure](https://github.com/bio-ontology-research-group/gspa/blob/main/tutorials/workflows/README.md#refresh-input-data-deliberately). Do not disable the check |
| Browser does not open on a server | Copy the generated HTML file to your computer and open it locally |

Repeating setup reuses Docker build layers and checks existing downloads. All
inputs, logs and results live under `.tutorial-work/`. Repeating a lesson replaces
that lesson's results, so copy them elsewhere first if you want to compare runs.

??? details "What is installed?"

    The CPU image builds GSPA from the current checkout on Ubuntu 24.04 and
    installs Java 21, Python, Prodigal and DIAMOND. The launcher limits each analysis
    container to two CPUs and 6 GiB RAM, with a 4 GiB Java heap, and writes files as
    your user. It records the source snapshot in `.tutorial-work/source-manifest.json`.
    See `tutorials/workflows/Dockerfile` for the build recipe. For an installation
    outside the lessons, follow the [installation guide](../installation.md).
