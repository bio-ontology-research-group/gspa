# Download resources and models

GSPA has one resource catalog for the standalone jar and source checkout:

```bash
./gspa resources                         # list IDs, types and download methods
./gspa resources clean                   # sources, dependencies, paths and setup notes
./gspa resources --json                  # machine-readable catalog
./gspa download --all --dry-run          # inspect everything before downloading
./gspa download --all --dir /data/gspa    # attempt every catalog entry
```

The manager needs **Python 3.9+**, using only the standard library. Choose an
interpreter with `--python /path/to/python`. Resource-specific downloaders need
their corresponding installed tools; the manager does not install Python
packages, system software, or accept third-party license agreements.

`--all` includes optional databases and models, including very large resources
(such as GTDB-Tk, over 100 GB compressed for the pinned release). It attempts
independent entries even when another fails. A summary lists everything that
failed or requires manual provision, and the command exits nonzero if anything
remains unavailable. A successful file download does not mean a predictor's
software environment has been installed or tested.

## Individual resources and groups

```bash
./gspa download go ec2go --dir /data/gspa
./gspa download diamond --dir /data/gspa
./gspa download esm2-650m --dir /data/gspa
./gspa download deepfri --dir /data/gspa
./gspa download deepgo-plusplus-light --dir /data/gspa
./gspa download deepgo-plusplus-models --dir /data/gspa
./gspa download --group core --dir /data/gspa
./gspa download --group models --dir /data/gspa
```

| Group | Contents |
|---|---|
| `core` | GO OWL, EC2GO and GSPA's bundled taxonomy, anatomy, complex and essential-function tables |
| `annotation` | GO, annotation-bearing Swiss-Prot, DIAMOND index and predictor sidecars |
| `models` | Every cataloged model, including individual DG++ variants and models needing manual provision |
| `databases` | Optional sequence, structure, domain, pathway and taxonomy resources |
| `bundled` | Reference tables, sidecars and frozen integrators shipped inside the jar; works offline |
| `deepgo-plusplus` | All frozen integrators, runtime assets, optional CNN and sidecars; unpublished assets are reported |

Dependencies are selected automatically: `download diamond` first obtains
Swiss-Prot, preserves GO/EC annotations in FASTA headers, then runs
`diamond makedb`. `download mmseqs2` builds the alternative index. The tools must
be installed. This annotated reference is needed because ordinary Swiss-Prot
FASTA titles do not contain the labels GSPA's hit-title parser transfers.

Resources default to `$GSPA_DATA_DIR`, or `~/.local/share/gspa` if unset.
`--dir` (alias `--dest`) overrides it. Each ID gets a separate subdirectory.
Downloading a group does not change any biological preset.

## Model coverage

Every JSON model in `deepgo-plusplus/models/` has its own catalog ID (filename
without `.json`, with underscores replaced by hyphens and lowercased).
`deepgo-plusplus-models` exports all of them into one directory suitable for
`--deepgo-plusplus-light-models` or the full integrator's model file.
These are the actual shipped coefficients, not links to an unverified mirror.

Other entries cover ESM-2 35M and 650M, ESM-1b, ProtT5, ProstT5, ProteInfer,
CLEAN, DeepFRI, metagenomic-DeepFRI, TMbed, DeepSig, TPpred3, Metapredict,
MusiteDeep, ScanNet, DeepARG and PhiSpy. Some download pinned source archives
that include weights; others call the official tool downloader. Runtime
libraries still need installation. Details identify which assets are supplied
and which paths or environment settings their wrappers consume.

For example, CLEAN needs its upstream source, pretrained weights/centroids and
ESM-1b; the catalog retrieves these separately rather than mistaking the source
archive for the complete model. The pretrained package uses `gdown`, with
upstream layout validation. See the [CLEAN setup instructions](https://github.com/tttianhao/CLEAN).
ProteInfer uses the upstream default single model per task (Pfam, EC and GO),
plus metadata; its optional larger ensembles are not included. See the
[ProteInfer installer](https://github.com/google-research/proteinfer/blob/master/install_models.py).

## Use downloaded paths

After an actual download, the destination contains:

- `paths.json`: resource directories and usage hints, including integration-only inputs.
- `gspa-paths.yaml`: annotation configuration paths. It uses JSON syntax, which
  the YAML loader accepts. It sets paths but **does not enable predictors**.
- `<ID>/.gspa-resource.json`: source/version, catalog fingerprint and content
  hashes for the installed files.

For example, after downloading the `annotation` group:

```bash
./gspa annotate --preset bacteria -i genome.fna \
  -c /data/gspa/gspa-paths.yaml \
  --db /data/gspa/diamond/swissprot.dmnd -o annotation/

./gspa integrate --preset bacteria --claims claims.jsonl --out integrated.tsv \
  --go-owl /data/gspa/go/go.owl --ec2go /data/gspa/ec2go/ec2go
```

Sequence annotation still needs an installed gene caller for raw nucleotide
input. The path file does not override explicit CLI options. For neural models,
select `--base-predictor` or the individual predictor flag after preparing its
runtime environment. Exported sidecars do not replace the repository's training
and inference pipeline code used by DeepGO-PlusPlus.

## Unpublished, custom and licensed resources

Some resources do **not** currently have a verified public complete download:
GSPA-specific ESM2 heads and vocabularies, centroid databases, the DG++ Light
runtime training assets and CNN checkpoint, and the complete flux-refiner
matrix bundle. The previous Nextflow manifest used placeholder URLs for several
of these; those URLs are not used by the new manager.

These entries remain individually addressable. `gspa resources ID` gives their
required files, build instructions or upstream distribution page. Provide a
compatible bundle from your own release with:

```bash
./gspa download deepgo-plusplus-light-assets \
  --source-url https://your-server.example/releases/light-assets.tar.gz \
  --sha256 YOUR_64_HEX_DIGIT_SHA256 --dir /data/gspa
```

Tar bundles must put the files at the paths shown in the resource's `expected`
field (for example, `train_db.dmnd` at the archive root). Single checkpoints
such as `deepgo-plusplus-cnn` accept a file URL instead of an archive.
The same `--source-url` option can override the URL of a single-file public
resource, for example to select an archived GO release or your institutional
mirror. Its compression/layout must match the catalog entry. Supplying a new
URL without `--sha256` does not reuse the old URL's expected checksum.

SignalP/TMHMM, PSORTb and DeepEC follow their upstream distribution requirements.
The Stage 3 agent uses a separately served model; the catalog does not guess a
Hugging Face model from an API model alias. Genome-specific claims, operons,
gaps, layouts, component predictions and learned `theta.json` are run inputs,
not universally downloadable reference files.

## Integrity, retries and updates

Completed resources are reused only after checking their recorded file sizes
and SHA-256 hashes. `--force` fetches them again. An existing unmanaged resource
directory is never overwritten. Failed replacements leave the previous
installation intact. Downloads retry transient failures, reject empty/truncated
responses and HTML error pages, and install only after checking expected files.
Archive extraction rejects escaping paths, external links and device entries.

The cache under `.downloads/` retains completed transfers, allowing a rerun to
reuse them after extraction or an external tool fails. Interrupted transfers
restart from the beginning; byte-range resume is not implemented. Budget disk
space for the cache, extracted files and temporary replacement directory.

Upstream SHA-256 checksums are checked where available; shipped models also
have checksums. Otherwise the receipt records the hash of the bytes actually
retrieved. **A locally recorded hash is not an upstream signature.** Moving
upstream releases are identified as such; `--force` is required to refresh an
already verified copy. Mutable download commands select releases according to
the installed upstream tool; inspect the resulting files/version information.

A destination lock prevents concurrent installs from colliding. If the process
is forcibly killed, confirm that it has stopped before removing the stale
`.download.lock` directory. Do not delete it while another download is running.

In a source checkout, the same manager is directly available as:

```bash
python3 scripts/gspa_resources.py list
python3 scripts/gspa_resources.py download --group bundled --dir /data/gspa
```
