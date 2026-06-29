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
FROM eclipse-temurin:21-jdk AS builder
WORKDIR /build
COPY gradlew settings.gradle.kts build.gradle.kts gradle.properties ./
COPY gradle/ gradle/
COPY gspa-core/ gspa-core/
COPY gspa-predictors/ gspa-predictors/
COPY gspa-cli/ gspa-cli/
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
        python3 python3-pip curl ca-certificates libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# DIAMOND (matches the version used to build train_db.dmnd).
ARG DIAMOND_VERSION=2.1.11
RUN curl -fL "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz" \
      | tar -xz -C /usr/local/bin diamond && diamond --version

WORKDIR /app

# Python deps: the service itself + numpy for the sidecar's integrator math.
COPY service/requirements.txt /app/service/requirements.txt
RUN pip install --break-system-packages -r /app/service/requirements.txt "numpy>=1.26,<3"
# Optional CPU torch + fair-esm: ONLY the esm2-kNN / CNN light components need them;
# the default DG++Light path (DIAMOND BLAST-KNN + STRING Net-KNN) is torch-free and
# they are lazy-imported. Best-effort + unpinned so the image always builds across
# base-image Python versions (the torch CPU index drops old pins on new Pythons).
RUN pip install --break-system-packages --index-url https://download.pytorch.org/whl/cpu torch \
      && pip install --break-system-packages fair-esm "scikit-learn~=1.5" \
    || echo "[gspa] optional torch/esm components unavailable; default light path unaffected"

# gspa-cli distribution on PATH.
COPY --from=builder /build/gspa-cli/build/install/gspa-cli /opt/gspa
ENV PATH="/opt/gspa/bin:${PATH}"

# The sidecar + DG++Light models the light predictor needs.
COPY benchmark/ /app/benchmark/
COPY deepgo-plusplus/ /app/deepgo-plusplus/
COPY service/ /app/service/
RUN chmod +x /app/service/entrypoint.sh /app/deepgo-plusplus/service/dgpp_server.py

EXPOSE 8000
# entrypoint launches the warm dgpp_light server, then uvicorn (see entrypoint.sh).
CMD ["/app/service/entrypoint.sh"]
