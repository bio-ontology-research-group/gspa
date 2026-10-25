# What can sequence embeddings suggest about protein function?

A protein language model represents sequences as numerical vectors, or
*embeddings*. Could similarity between those vectors help us nominate functions
for proteins?

**You will learn what an embedding-based candidate means, and distinguish model
output from biological validation.** This optional lesson follows
[setup](install.md). It downloads an ESM-2 model and builds a separate CPU neural
image; first execution can take substantially longer than a biological lesson.

## 1. Predict the limitation

We use 64 short reference proteins—16 each from bacteria, archaea, yeast and
phage—to build function centroids, the average embeddings for labelled groups.
We then query 16 different proteins, four from each group. What happens if a
query's true function is absent from such a small reference?

## 2. Run the CPU experiment

```bash
./tutorial neural cpu
```

The launcher downloads the real 35M-parameter ESM-2 backbone, prepares the
reference, builds the centroids and runs GSPA's centroid predictor. No GPU is
needed. The recorded run builds **28 centroids** and produces **48 candidate
annotations for 16 query proteins**: three per protein.

Open `.tutorial-work/neural/cpu/*_annotations.tsv` in a text editor or spreadsheet
(the `*` denotes the generated filename). Select a protein and compare its three
labels and scores. The source column identifies `esm2-centroid`.

**Biological insight:** embedding similarity provides a shortlist of functional
hypotheses relative to the functions represented in the reference. The fixed
three-candidate output does not mean every protein has exactly three functions.
An unrepresented function cannot be recovered simply by ranking these centroids.

## 3. Check your interpretation

**Does the highest cosine-similarity score establish the correct function?**

??? success "Worked answer"

    No. Cosine similarity measures proximity in the model's representation, not
    the probability that a function is correct. This tiny custom reference is
    neither a released functional prediction head nor a calibrated model. Reference
    and query IDs differ, but we have not controlled their homology. Assessing
    functional accuracy needs suitable held-out labels and evaluation design.

## Optional: does changing the device change these candidates?

After the CPU run, a host with a compatible NVIDIA GPU, driver and
[NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)
can run:

```bash
./tutorial neural cuda
```

This builds the CUDA image and uses the **same query proteins and centroids**.
The check requires identical top-three labels and score differences no larger
than 0.001. The recorded run met those conditions. Your result is written to
`.tutorial-work/neural/comparison.json`; a missing GPU or failed comparison
reports an error rather than silently using the CPU.

Matching device outputs supports computational reproducibility. It does not
validate the biological annotations. The GPU accelerates the neural stage;
GSPA's integration and the earlier DIAMOND examples remain CPU workloads.

??? details "Environment and production use"

    The recipe uses PyTorch 2.7.1, fair-esm 2.0.0 and NumPy 1.26.4. CPU and CUDA
    environments are separate images; the latter uses CUDA 12.8 wheels. Consult
    [PyTorch's installation guide](https://pytorch.org/get-started/locally/) for
    hardware compatibility. Short sequences and batch size four keep the example
    small. Production use needs an appropriate function reference or trained head,
    calibration and independent evaluation; backbone weights alone do not assign
    functions. See [resources and models](../reference/resources.md).
