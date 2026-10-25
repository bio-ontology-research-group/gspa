# GSPA genome-scale annotation service.
#
# Builds the gspa-cli JVM distribution, then a CPU-only runtime that can run the
# whole `gspa annotate` genome-scale pipeline behind a small FastAPI endpoint:
#   genome FASTA + GFF3 -> CDS translation -> DeepGO-PlusPlus-Light (CPU)
#   -> per-contig genome-scale metrics -> optional SAT taxon-consistency /
#      completeness / coherence enforcement -> provenance.
#
# Build context = the gspa repo ROOT (needs gspa-*/ + benchmark/ + deepgo-plusplus/):
#   docker build -f service/Dockerfile -t gspa-service .
#
# The heavy DG++Light asset bundle is NOT baked in; mount it at /opt/dgpp_assets
# (the same volume DeepGOWeb's download_assets.sh populates).

# --- Stage 1: build the gspa-cli distribution -------------------------------
FROM gspa-analysis-tools:2026-10-25 AS analysis_tools
FROM eclipse-temurin:21-jdk AS builder
WORKDIR /build
COPY gradlew settings.gradle.kts build.gradle.kts gradle.properties ./
COPY gradle/ gradle/
COPY gspa-core/ gspa-core/
COPY gspa-predictors/ gspa-predictors/
COPY gspa-cli/ gspa-cli/
COPY scripts/ scripts/
COPY benchmark/neural/ benchmark/neural/
COPY deepgo-plusplus/models/ deepgo-plusplus/models/
RUN ./gradlew :gspa-cli:installDist --no-daemon

# --- Stage 2: CPU runtime ----------------------------------------------------
FROM eclipse-temurin:21-jre

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    TORCH_HOME=/opt/torch \
    DGPP_ASSETS=/opt/dgpp_assets \
    NEURAL_SIDECAR=/app/benchmark/neural/run_neural_predictors.py \
    DGPP_LIGHT_SERVER=/tmp/dgpp_light.sock \
    DGPP_LIGHT_THREADS=8

# System deps: python (sidecar), diamond (BLAST-KNN), libgomp (torch).
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip curl ca-certificates libgomp1 prodigal \
    && rm -rf /var/lib/apt/lists/*

# DIAMOND (matches the version used to build train_db.dmnd).
ARG DIAMOND_VERSION=2.1.11
RUN curl -fL "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz" \
      | tar -xz -C /usr/local/bin diamond && diamond --version

COPY --from=analysis_tools /opt/bio /opt/bio
ENV PATH="/opt/bio/bin:${PATH}"
WORKDIR /app

# Python deps: the service itself + numpy for the sidecar's integrator math.
COPY service/requirements.txt /app/service/requirements.txt
RUN pip install --break-system-packages -r /app/service/requirements.txt "numpy>=1.26,<3"
# Compile antiSMASH's packaged HMM/SVM resources during build: the service root is read-only.
RUN python -c 'from antismash.config import build_config; from antismash.main import get_all_modules; build_config([], isolated=True, modules=get_all_modules()); [m.prepare_data(True) for m in get_all_modules() if hasattr(m, "prepare_data") and m.__name__ not in ("antismash.modules.cluster_compare", "antismash.modules.clusterblast")]'
# Keep the validated HiGHS/SciPy wheels separate from antiSMASH's conda pins.
RUN python -m venv /opt/fba && /opt/fba/bin/pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3 python-libsbml==5.21.1
# Optional torch + fair-esm: ONLY the esm2-kNN / CNN components need them; the
# default DG++Light path (DIAMOND BLAST-KNN + STRING Net-KNN) is torch-free.
# Enable optional neural components with --build-arg INSTALL_TORCH=true.
#
# TORCH_INDEX selects the wheel build: default CPU; pass
#   --build-arg INSTALL_TORCH=true --build-arg TORCH_INDEX=https://download.pytorch.org/whl/cu128
# for a GPU image (CUDA 12.8 / Blackwell sm_120) that runs the esm2-kNN + CNN
# stream on the GPU under DGppLight.predict_full (DGPP_LIGHT_DEVICE=cuda),
# overlapping the CPU DIAMOND/Net-KNN stream. Run such a container with `--gpus all`.
ARG TORCH_INDEX=https://download.pytorch.org/whl/cpu
ARG INSTALL_TORCH=false
RUN if [ "$INSTALL_TORCH" = true ]; then \
      pip install --break-system-packages --index-url ${TORCH_INDEX} torch \
      && pip install --break-system-packages fair-esm "scikit-learn~=1.5"; fi

# gspa-cli distribution on PATH.
COPY --from=builder /build/gspa-cli/build/install/gspa-cli /opt/gspa
ENV PATH="/opt/gspa/bin:${PATH}"

# The sidecar + DG++Light models the light predictor needs.
COPY benchmark/ /app/benchmark/
COPY deepgo-plusplus/ /app/deepgo-plusplus/
COPY scripts/stoichiometric_locus_refiner.py /app/scripts/stoichiometric_locus_refiner.py
COPY service/ /app/service/
RUN chmod +x /app/service/entrypoint.sh /app/deepgo-plusplus/service/dgpp_server.py

EXPOSE 8000
# entrypoint launches the warm dgpp_light server, then uvicorn (see entrypoint.sh).
CMD ["/app/service/entrypoint.sh"]
